#!/usr/bin/env node
/**
 * Stellome builder
 *
 * Downloads the human STRING network and Reactome pathway mappings, groups every
 * protein into a biological "region" (its dominant Reactome top-level pathway),
 * computes a clustered 3D layout and writes compact files the browser can load:
 *
 *   Stellome.html                  one self-contained file: open it in any browser, no server needed
 *   public/stellome/stellome.json proteins, constellations, metadata   (used by the in-app view)
 *   public/stellome/positions.bin  Float32 x,y,z per protein
 *   public/stellome/edges.bin      Uint16 [a, b, score] per interaction
 *
 * Usage:   npm run build                      (default: STRING score >= 700)
 *          node scripts/build-stellome.mjs --min-score 900   (fewer, stronger links, smaller file)
 *          node scripts/build-stellome.mjs --ticks 400       (longer, smoother layout)
 *          node scripts/build-stellome.mjs --no-uniprot      (STRING proteins only, skip the UniProt additions)
 *
 * UniProt additions: every reviewed (Swiss-Prot) human protein that STRING does not contain is added
 * as a star without links, placed in its Reactome constellation. Unreviewed TrEMBL entries are not
 * added (about 60,000 mostly redundant isoforms/fragments; they would also exceed the 16-bit index limit).
 *
 * Downloads (~150 MB) are cached in .stellome-cache/ and reused on later runs.
 * If a download is blocked, put the files listed below into .stellome-cache/ by hand.
 */
import { createReadStream, createWriteStream, existsSync, mkdirSync, readFileSync, renameSync, writeFileSync, statSync } from 'node:fs';
import { createGunzip } from 'node:zlib';
import { createInterface } from 'node:readline';
import { Readable } from 'node:stream';
import { pipeline } from 'node:stream/promises';
import { join, dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { forceSimulation, forceManyBody, forceLink, forceX, forceY, forceZ } from 'd3-force-3d';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');

/* ───────────── options ───────────── */
const argv = process.argv.slice(2);
const opt = (name, dflt) => {
  const i = argv.indexOf(`--${name}`);
  return i >= 0 && argv[i + 1] ? argv[i + 1] : dflt;
};
const MIN_SCORE = Number(opt('min-score', 700)); // STRING combined score, 0–1000
const LAYOUT_SCORE = Math.max(MIN_SCORE, Number(opt('layout-score', 900))); // links that shape the layout
const TICKS = Number(opt('ticks', 300));
const CACHE = resolve(ROOT, opt('cache', '.stellome-cache'));
const OUT = resolve(ROOT, opt('out', 'public/stellome'));
const STRING_V = opt('string-version', 'v12.0');
const HTML_OUT = resolve(ROOT, opt('html', 'Stellome.html'));
const WITH_UNIPROT = !argv.includes('--no-uniprot');

const SRC = {
  links: { file: `9606.protein.links.${STRING_V}.txt.gz`, url: `https://stringdb-downloads.org/download/protein.links.${STRING_V}/9606.protein.links.${STRING_V}.txt.gz` },
  info: { file: `9606.protein.info.${STRING_V}.txt.gz`, url: `https://stringdb-downloads.org/download/protein.info.${STRING_V}/9606.protein.info.${STRING_V}.txt.gz` },
  aliases: { file: `9606.protein.aliases.${STRING_V}.txt.gz`, url: `https://stringdb-downloads.org/download/protein.aliases.${STRING_V}/9606.protein.aliases.${STRING_V}.txt.gz` },
  pathways: { file: 'ReactomePathways.txt', url: 'https://reactome.org/download/current/ReactomePathways.txt' },
  relations: { file: 'ReactomePathwaysRelation.txt', url: 'https://reactome.org/download/current/ReactomePathwaysRelation.txt' },
  ensembl2r: { file: 'Ensembl2Reactome_All_Levels.txt', url: 'https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt' },
  uniprot: { file: 'uniprot_human_reviewed.tsv', url: 'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=tsv&fields=accession%2Cid%2Cgene_primary%2Cprotein_name%2Cxref_string&query=%28proteome%3AUP000005640%29+AND+%28reviewed%3Atrue%29' },
  uniprot2r: { file: 'UniProt2Reactome_All_Levels.txt', url: 'https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt' },
};

const t0 = Date.now();
const log = (...a) => console.log(`[${((Date.now() - t0) / 1000).toFixed(1).padStart(6)}s]`, ...a);

/* ───────────── download + read helpers ───────────── */
async function ensure(key) {
  const { file, url } = SRC[key];
  const path = join(CACHE, file);
  if (existsSync(path) && statSync(path).size > 0) return path;
  mkdirSync(CACHE, { recursive: true });
  log(`Downloading ${file} …`);
  const res = await fetch(url, { headers: { 'User-Agent': 'stellome (open-source research tool)' } }).catch((e) => ({ ok: false, statusText: e.message }));
  if (!res.ok) {
    throw new Error(
      `Could not download ${url} (${res.status ?? ''} ${res.statusText}).\n` +
        `Download it in your browser and save it as:\n  ${path}\nthen run this command again.`,
    );
  }
  const total = Number(res.headers.get('content-length')) || 0;
  let got = 0, lastPct = -10;
  const body = Readable.fromWeb(res.body);
  body.on('data', (c) => {
    got += c.length;
    const pct = total ? Math.floor((got / total) * 100) : -1;
    if (pct >= lastPct + 10) { lastPct = pct; log(`  ${file}: ${pct}% of ${(total / 1e6).toFixed(0)} MB`); }
  });
  await pipeline(body, createWriteStream(path + '.part'));
  renameSync(path + '.part', path);
  return path;
}

async function* lines(path) {
  const input = path.endsWith('.gz') ? createReadStream(path).pipe(createGunzip()) : createReadStream(path);
  const rl = createInterface({ input, crlfDelay: Infinity });
  for await (const line of rl) if (line) yield line;
}

/* deterministic PRNG so the map is identical on every build */
function mulberry32(a) {
  return () => {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const rand = mulberry32(9606);

/* ───────────── 1 · proteins ───────────── */
async function main() {
  const infoPath = await ensure('info');
  const proteins = []; // { ensp, sym, ann }
  const idx = new Map();
  for await (const l of lines(infoPath)) {
    if (l.startsWith('#') || l.startsWith('string_protein_id')) continue;
    const [id, sym, , ann = ''] = l.split('\t');
    const ensp = id.replace(/^9606\./, '');
    idx.set(ensp, proteins.length);
    const clean = ann.replace(/\s*\[[^\]]*\]/g, '').replace(/\s+/g, ' ').trim();
    proteins.push({ ensp, sym, ann: clean.length > 220 ? clean.slice(0, 217).replace(/\s\S*$/, '') + '…' : clean });
  }
  /* UniProt-only proteins: reviewed human entries that STRING lacks (no links, placed by Reactome). */
  const extraIdx = new Map(); // UniProt accession -> protein index
  if (WITH_UNIPROT) {
    try {
      const have = new Set(proteins.map((p) => p.sym.toUpperCase()));
      let seen = 0;
      for await (const l of lines(await ensure('uniprot'))) {
        const c = l.split('\t');
        if (c[0] === 'Entry' || c.length < 4) continue;
        seen++;
        const [acc, , gene = '', pname = '', str = ''] = c;
        if (str.split(';').some((x) => x && idx.has(x.replace(/^9606\./, '')))) continue; // already in STRING
        const sym = gene.split(/[ ;]/)[0].trim();
        if (sym && have.has(sym.toUpperCase())) continue; // same gene under another identifier
        const name = pname.replace(/\s*\(.*$/, '').trim();
        extraIdx.set(acc, proteins.length);
        proteins.push({ ensp: '', sym: sym || acc, ann: name, extra: 1, acc });
        if (sym) have.add(sym.toUpperCase());
      }
      log(`UniProt: ${seen.toLocaleString()} reviewed human entries, ${extraIdx.size.toLocaleString()} not in STRING were added`);
    } catch (e) {
      log(`UniProt step skipped: ${e.message}`);
    }
  }
  if (proteins.length > 65535) throw new Error('More than 65,535 proteins: edges.bin uses 16-bit indices.');
  log(`${proteins.length.toLocaleString()} proteins`);

  /* ───────────── 2 · interactions ───────────── */
  const linksPath = await ensure('links');
  const ea = [], eb = [], es = [];
  for await (const l of lines(linksPath)) {
    if (l.startsWith('protein1')) continue;
    const [p1, p2, s] = l.split(' ');
    const score = +s;
    if (score < MIN_SCORE) continue;
    const a = idx.get(p1.slice(5)), b = idx.get(p2.slice(5));
    if (a === undefined || b === undefined || a >= b) continue; // file lists both directions; keep one
    ea.push(a); eb.push(b); es.push(score);
  }
  const E = ea.length, N = proteins.length;
  const deg = new Uint32Array(N);
  for (let i = 0; i < E; i++) { deg[ea[i]]++; deg[eb[i]]++; }
  log(`${E.toLocaleString()} interactions with score ≥ ${MIN_SCORE}`);

  const adj = Array.from({ length: N }, () => []);
  for (let i = 0; i < E; i++) { adj[ea[i]].push(eb[i]); adj[eb[i]].push(ea[i]); }

  /* ───────────── 3 · Reactome regions ───────────── */
  const pwName = new Map();
  for await (const l of lines(await ensure('pathways'))) {
    const [st, name, sp] = l.split('\t');
    if (sp === 'Homo sapiens') pwName.set(st, name);
  }
  const hasParent = new Set();
  for await (const l of lines(await ensure('relations'))) hasParent.add(l.split('\t')[1]);
  const topLevel = new Set([...pwName.keys()].filter((st) => !hasParent.has(st)));
  log(`${topLevel.size} Reactome top-level pathways`);

  const tops = Array.from({ length: N }, () => new Set());
  let mapped = 0;
  const addTop = (i, st) => { if (topLevel.has(st)) tops[i].add(st); };

  for await (const l of lines(await ensure('ensembl2r'))) {
    if (!l.startsWith('ENSP')) continue;
    const c = l.split('\t');
    if (c[5] !== 'Homo sapiens') continue;
    const i = idx.get(c[0].split('.')[0]);
    if (i !== undefined) addTop(i, c[1]);
  }
  mapped = tops.filter((s) => s.size).length;
  log(`Reactome via Ensembl protein IDs: ${mapped.toLocaleString()} proteins`);

  if (extraIdx.size) {
    for await (const l of lines(await ensure('uniprot2r'))) {
      const c = l.split('\t');
      if (c[5] !== 'Homo sapiens') continue;
      const i = extraIdx.get(c[0]);
      if (i !== undefined) addTop(i, c[1]);
    }
    mapped = tops.filter((s) => s.size).length;
    log(`Reactome after adding UniProt-only proteins: ${mapped.toLocaleString()} proteins`);
  }

  if (mapped < N * 0.3) {
    log('Low coverage, adding UniProt-based mapping …');
    const uni2i = new Map();
    for await (const l of lines(await ensure('aliases'))) {
      const c = l.split('\t');
      if (!c[2] || !c[2].includes('UniProt_AC')) continue;
      const i = idx.get(c[0].replace(/^9606\./, ''));
      if (i !== undefined && !uni2i.has(c[1])) uni2i.set(c[1], i);
    }
    for await (const l of lines(await ensure('uniprot2r'))) {
      const c = l.split('\t');
      if (c[5] !== 'Homo sapiens') continue;
      const i = uni2i.get(c[0]);
      if (i !== undefined) addTop(i, c[1]);
    }
    mapped = tops.filter((s) => s.size).length;
    log(`Reactome after UniProt mapping: ${mapped.toLocaleString()} proteins`);
  }

  // "Disease" re-annotates normal biology; use it only when nothing else applies.
  const DISEASE = [...topLevel].find((st) => pwName.get(st) === 'Disease');
  for (const s of tops) if (s.size > 1 && DISEASE) s.delete(DISEASE);

  // Pick each protein's region: the top-level pathway it shares with most of its partners.
  const freq = new Map();
  for (const s of tops) for (const st of s) freq.set(st, (freq.get(st) ?? 0) + 1);
  const region = new Array(N).fill(null);
  for (let i = 0; i < N; i++) {
    if (!tops[i].size) continue;
    let best = null, bestScore = -1;
    for (const st of tops[i]) {
      let sc = 0;
      for (const j of adj[i]) if (tops[j].has(st)) sc++;
      sc += 0.5 / Math.log2(2 + (freq.get(st) ?? 0)); // tie-break toward the more specific region
      if (sc > bestScore) { bestScore = sc; best = st; }
    }
    region[i] = best;
  }
  // Unannotated proteins inherit the most common region among their partners (a few rounds).
  for (let round = 0; round < 4; round++) {
    const next = region.slice();
    for (let i = 0; i < N; i++) {
      if (region[i] || !adj[i].length) continue;
      const votes = new Map();
      for (const j of adj[i]) if (region[j]) votes.set(region[j], (votes.get(region[j]) ?? 0) + 1);
      let best = null, bv = 0;
      for (const [st, v] of votes) if (v > bv) { bv = v; best = st; }
      if (best) next[i] = best;
    }
    region.splice(0, N, ...next);
  }
  const inferred = region.filter((r, i) => r && !tops[i].size).length;
  const UNASSIGNED = 'NONE';
  for (let i = 0; i < N; i++) if (!region[i]) region[i] = UNASSIGNED;

  // Region table, largest first.
  const counts = new Map();
  for (const r of region) counts.set(r, (counts.get(r) ?? 0) + 1);
  const regions = [...counts.entries()]
    .sort((a, b) => (a[0] === UNASSIGNED) - (b[0] === UNASSIGNED) || b[1] - a[1])
    .map(([id, size]) => ({ id, name: id === UNASSIGNED ? 'No pathway annotation' : pwName.get(id), size }));
  const rIndex = new Map(regions.map((r, k) => [r.id, k]));
  const cat = Uint8Array.from(region, (r) => rIndex.get(r));
  log(`${regions.length} regions; ${inferred.toLocaleString()} proteins placed by their partners`);

  /* ───────────── 4 · region centres (small force layout of regions) ───────────── */
  const R = regions.length;
  const w = Array.from({ length: R }, () => new Float64Array(R));
  for (let i = 0; i < E; i++) { const a = cat[ea[i]], b = cat[eb[i]]; if (a !== b) { w[a][b]++; w[b][a]++; } }
  const radius = regions.map((r) => 38 * Math.cbrt(r.size));
  const C = regions.map((_, k) => {
    const y = 1 - ((k + 0.5) / R) * 2, rr = Math.sqrt(1 - y * y), th = k * 2.399963;
    return [Math.cos(th) * rr * 900, y * 900, Math.sin(th) * rr * 900];
  });
  for (let it = 0; it < 600; it++) {
    const F = C.map(() => [0, 0, 0]);
    for (let a = 0; a < R; a++) for (let b = a + 1; b < R; b++) {
      const d = [C[b][0] - C[a][0], C[b][1] - C[a][1], C[b][2] - C[a][2]];
      const dist = Math.hypot(...d) || 1, gap = radius[a] + radius[b] + 60;
      const wn = w[a][b] / Math.sqrt(regions[a].size * regions[b].size);
      let f = dist < gap ? (dist - gap) * 0.25 : 0; // never overlap
      f += Math.min(1.5, wn * 3) * (dist - gap) * 0.01; // linked regions sit closer
      f -= 4e4 / (dist * dist); // general spread
      for (let k = 0; k < 3; k++) { F[a][k] += (d[k] / dist) * f; F[b][k] -= (d[k] / dist) * f; }
    }
    for (let a = 0; a < R; a++) for (let k = 0; k < 3; k++) C[a][k] += Math.max(-40, Math.min(40, F[a][k] - C[a][k] * 0.002));
  }

  /* ───────────── 5 · protein layout ───────────── */
  const nodes = proteins.map((_, i) => {
    const c = C[cat[i]], r = radius[cat[i]] * Math.cbrt(rand());
    const u = rand() * 2 - 1, th = rand() * Math.PI * 2, s = Math.sqrt(1 - u * u);
    return { x: c[0] + r * s * Math.cos(th), y: c[1] + r * u, z: c[2] + r * s * Math.sin(th) };
  });
  const links = [];
  for (let i = 0; i < E; i++) if (es[i] >= LAYOUT_SCORE) links.push({ source: ea[i], target: eb[i], same: cat[ea[i]] === cat[eb[i]] });
  log(`Laying out ${N.toLocaleString()} proteins with ${links.length.toLocaleString()} links (${TICKS} steps) …`);

  const sim = forceSimulation(nodes, 3)
    .force('charge', forceManyBody().strength(-6).theta(0.9).distanceMax(400))
    .force('link', forceLink(links).distance((l) => (l.same ? 14 : 60)).strength((l) => (l.same ? 0.25 : 0.01)).iterations(1))
    .force('x', forceX((d, i) => C[cat[i]][0]).strength(0.06))
    .force('y', forceY((d, i) => C[cat[i]][1]).strength(0.06))
    .force('z', forceZ((d, i) => C[cat[i]][2]).strength(0.06))
    .alphaDecay(1 - Math.pow(0.001, 1 / TICKS))
    .stop();
  for (let t = 1; t <= TICKS; t++) {
    sim.tick();
    if (t % 25 === 0) log(`  layout ${Math.round((t / TICKS) * 100)}%`);
  }

  // Final region centroids (for labels) and bounds.
  const cent = regions.map(() => [0, 0, 0]);
  nodes.forEach((n, i) => { const c = cent[cat[i]]; c[0] += n.x; c[1] += n.y; c[2] += n.z; });
  cent.forEach((c, k) => { for (let j = 0; j < 3; j++) c[j] = +(c[j] / regions[k].size).toFixed(1); });

  /* ───────────── 6 · write ───────────── */
  mkdirSync(OUT, { recursive: true });
  const pos = new Float32Array(N * 3);
  nodes.forEach((n, i) => { pos[i * 3] = n.x; pos[i * 3 + 1] = n.y; pos[i * 3 + 2] = n.z; });
  writeFileSync(join(OUT, 'positions.bin'), Buffer.from(pos.buffer));

  // Sort edges by score descending so the viewer can take a prefix for any threshold.
  const order = Array.from({ length: E }, (_, i) => i).sort((p, q) => es[q] - es[p]);
  const eb16 = new Uint16Array(E * 3);
  order.forEach((i, k) => { eb16[k * 3] = ea[i]; eb16[k * 3 + 1] = eb[i]; eb16[k * 3 + 2] = es[i]; });
  writeFileSync(join(OUT, 'edges.bin'), Buffer.from(eb16.buffer));

  const meta = {
    version: 1,
    built: new Date().toISOString(),
    sources: { string: `STRING ${STRING_V} (Homo sapiens, 9606)`, reactome: 'Reactome current release', ...(extraIdx.size ? { uniprot: 'UniProt reviewed human proteome (UP000005640)' } : {}) },
    minScore: MIN_SCORE,
    proteins: N,
    interactions: E,
    regions: regions.map((r, k) => ({ id: r.id, name: r.name, size: r.size, centroid: cent[k] })),
    symbol: proteins.map((p) => p.sym),
    ensp: proteins.map((p) => p.ensp),
    annotation: proteins.map((p) => p.ann),
    region: Array.from(cat),
    degree: Array.from(deg),
    annotated: proteins.map((_, i) => (tops[i].size ? 1 : 0)),
    ...(extraIdx.size ? { extra: proteins.map((p) => (p.extra ? 1 : 0)), acc: proteins.map((p) => p.acc ?? '') } : {}),
  };
  writeFileSync(join(OUT, 'stellome.json'), JSON.stringify(meta));
  // One self-contained page with everything embedded.
  const tpl = readFileSync(join(ROOT, 'scripts', 'stellome-template.html'), 'utf8');
  const html = tpl
    .replace('__META__', () => JSON.stringify(meta).replace(/</g, '\\u003c'))
    .replace('__POS__', () => Buffer.from(pos.buffer).toString('base64'))
    .replace('__EDGES__', () => Buffer.from(eb16.buffer).toString('base64'));
  writeFileSync(HTML_OUT, html);

  const mb = (f) => (statSync(f).size / 1e6).toFixed(1);
  log(`Done → ${HTML_OUT} (${mb(HTML_OUT)} MB). Double-click it to open the map in your browser.`);
}

main().catch((e) => {
  console.error('\n✖ ' + (e?.message ?? e));
  process.exit(1);
});
