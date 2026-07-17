# Notebook documentation site (fixes #355)

Turn https://aarhuscosmology.github.io/CLASSpp/ from a bare Doxygen dump into a
site that is useful to end-users: a human landing page, executed-and-rendered
tutorial notebooks, a cosmology-model showcase, and the existing C++ API
reference — all in one GitHub Pages deployment.

## Decisions (made with Thomas, 2026-07-17)

- **Full site restructure**: landing page at `/`, Doxygen moves to `/api/`.
- **Notebooks are executed in CI** on the self-hosted runner during the
  nightly docs deploy — never rendered from committed outputs, so plots always
  reflect current master and a broken tutorial breaks the build visibly.
- **Curated content**: classic tutorials + best model showcases.
  PR-validation notebooks stay in the repo but off the site.
- **No Binder / JupyterLite**: classy is a compiled C++ extension, so
  in-browser execution is out; Binder needs a maintained Docker image and is
  slow. Download-`.ipynb` buttons instead. Revisit only if missed.

## Stack

Sphinx + **MyST-NB** (notebook execution + rendering) +
**pydata-sphinx-theme** (the NumPy/SciPy/pandas theme) + sphinx-design (card
grids) + sphinx-copybutton. Pure-pip, pinned in `doc/site/requirements.txt`.

Rejected: Jupyter Book 2 (Node toolchain on self-hosted runners, still-moving
ecosystem), nbconvert + hand-rolled index (hand-built nav for a worse site).

## Site layout

| URL | Content |
| --- | --- |
| `/` | Landing: what CLASSpp is, install quickstart, card grid → Tutorials / Models / C++ API / GitHub |
| `/tutorials/<slug>.html` | Classic notebooks (warmup, distances, thermo, …) |
| `/models/<slug>.html` | Model showcases (axion species, dcdm→wdm decay, …) |
| `/api/` | Doxygen HTML, byte-identical to today's output, moved down one level |

Existing deep links into the Doxygen pages at the root break; accepted.

## Source layout

New `doc/site/`:

- `conf.py` — Sphinx config **and the curation list**: a dict mapping section
  (`tutorials` / `models`) → list of `(slug, source-filename)` pairs. At build
  start, conf.py copies each curated notebook from `notebooks/` into
  `doc/site/<section>/<slug>.ipynb`. Adding a notebook to the site later =
  one line here.
- `index.md` — landing page (MyST markdown, sphinx-design cards).
- `tutorials/index.md`, `models/index.md` — section intros + toctrees.
- `requirements.txt` — pinned doc dependencies.
- `.gitignore` additions: `doc/site/_build/`, copied `doc/site/*/**.ipynb`.

Notebooks stay in `notebooks/` (single source of truth). Each curated
notebook must start with an H1 markdown cell (MyST-NB page title); the
verification pass adds one where missing.

## Execution model

- `nb_execution_mode` is `"off"` by default (fast local/PR builds) and
  switched to `"force"` via env var `CLASSPP_DOCS_EXECUTE=force` in the
  nightly job. No jupyter-cache: it keys on notebook source, not on the CLASS
  build, so caching would silently serve stale plots after code changes.
- `nb_execution_timeout = 900` seconds per notebook,
  `nb_execution_raise_on_error = True` — an erroring notebook fails the
  Sphinx build, which fails the deploy job; Pages keeps the last good site.

## CI changes

**`test_nightly.yml` / `deploy_documentation`** (self-hosted, linux, light):

1. Build Doxygen HTML as today (`doc/manual/html`).
2. Create venv (existing repo pattern), `pip install .` (builds classy),
   `pip install -r doc/site/requirements.txt`.
3. `CLASSPP_DOCS_EXECUTE=force sphinx-build -b html doc/site doc/site/_build/html`
4. `cp -r doc/manual/html doc/site/_build/html/api`
5. Upload/deploy Pages artifact from `doc/site/_build/html`.

**`docs.yml`** (PR check): keep the Doxygen build; add a Sphinx build with
execution off (venv + requirements only, no CLASS build) using `-W` so broken
toctrees/markup fail PRs. Relax `-W` only if it proves noisy.

## Cross-linking

- `doc/input/mainpage.md` gets a prominent link at the top: tutorials and
  examples live at `../index.html` (relative, works from `/api/`).
- Landing-page card and theme header link to `api/index.html`.

## Curation candidates (final list = what survives verification)

- Tutorials: warmup, distances, thermo, cltt_terms, one_k, one_time,
  many_times, varying_neff, neutrinohierarchy, varying_pann, Growth_with_w,
  cl_ST.
- Models: 370-axion-species, dcdm-wdm-massive-decay-products,
  173-Interacting_NCDM, 185-Hot_NEDE, Thomas_PPF — renamed to clean slugs
  (e.g. `axion_species`) by the copy step.

Verification pass: build CLASS + classy locally, execute every candidate,
fix cheap breakages in `notebooks/`, drop what cannot run against master.
Outcomes recorded in the PR description.

## Out of scope (future work)

- Python `classy` API reference via autodoc.
- Binder or a wasm classy build.
- Moving/retiring `doc/manual/CLASS_MANUAL.pdf`.
