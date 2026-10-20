# Contributing to Stellome

Thank you for helping. Biologists, bioinformaticians, designers and developers are all welcome.

## Ways to help

- **Report a problem**: a wrong annotation, a protein in the wrong pathway, a broken button. Use the issue templates.
- **Suggest a feature**: describe the scientific question it would help answer.
- **Propose a dataset**: it must be openly licensed (CC0, CC BY or similar) and cover the human proteome.
- **Send a pull request**: start with an issue for anything larger than a small fix, so we can agree on the approach first.

## How the project is built

| File | What it does |
|---|---|
| `scripts/build-stellome.mjs` | Downloads STRING, Reactome and UniProt, assigns proteins to pathways, computes the layout and writes `Stellome.html`. |
| `scripts/stellome-template.html` | The whole app (HTML, CSS, WebGL) with three placeholders, `__META__`, `__POS__` and `__EDGES__`, that the build fills with data. |
| `.github/workflows/site.yml` | Builds the map and publishes the website; runs on every push and monthly. |
| `.github/workflows/release.yml` | Attaches `Stellome.html` to each GitHub release. |

The app is plain JavaScript with raw WebGL and no framework, so it stays a single file that opens anywhere.

## Development

```bash
npm ci
npm run build:quick    # faster build for testing (links ≥ 0.90)
```

Open `Stellome.html` and check your change in both the circular and globe layouts, in light and dark themes, and at a narrow window width. Browser console errors count as bugs.

When you change the template, keep it working offline: no new network requests at start-up, and any external library must be optional and loaded from a pinned version.

## Pull requests

- One topic per pull request, with a short description of what changed and why.
- Include a screenshot for visual changes.
- For scientific changes (layout rules, statistics, data filters), explain the method and cite a reference.
- By contributing you agree that your code is released under the MIT licence.
