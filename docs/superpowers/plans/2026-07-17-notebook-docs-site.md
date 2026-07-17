# Notebook Documentation Site Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn https://aarhuscosmology.github.io/CLASSpp/ into a real documentation site: landing page at `/`, CI-executed notebooks at `/tutorials/` and `/models/`, Doxygen moved to `/api/` (fixes #355).

**Architecture:** A Sphinx project in `doc/site/` renders a MyST landing page plus notebooks copied at build time from `notebooks/` (single source of truth). Notebook execution is off by default and forced via `CLASSPP_DOCS_EXECUTE=force` in the nightly deploy job, which also copies the Doxygen HTML into `<site>/api/` before uploading the Pages artifact.

**Tech Stack:** Sphinx, MyST-NB, pydata-sphinx-theme, sphinx-design, sphinx-copybutton. Pure pip; no Node.

**Spec:** `docs/superpowers/specs/2026-07-17-notebook-docs-site-design.md`

## Global Constraints

- Branch: `docs/notebook-site` (already created off `origin/master`).
- NEVER `git add -A` or `git add .` in this repo — stage explicit paths only (in-source CMake/Xcode artifacts get swept in otherwise).
- Commit messages end with:
  `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_01LMSkqYYndXZdTB6sy8D7co`
- Notebooks stay in `notebooks/`; `doc/site/` only ever holds build-time copies of them (gitignored).
- The site must build with plain `python -m sphinx` and **no CLASS/classy installed** when `CLASSPP_DOCS_EXECUTE` is unset — that is the PR-check invariant.
- Env var name (exact): `CLASSPP_DOCS_EXECUTE`, active value `force`.
- Local venv name: `venv-docs` at repo root (`venv-docs/` is not gitignored — never stage it; the CI jobs create their own).
- PyPI package name is `classy-community` (import name `classy`).

---

### Task 1: Sphinx site skeleton

**Files:**
- Create: `doc/site/requirements.txt`
- Create: `doc/site/conf.py`
- Create: `doc/site/index.md`
- Create: `doc/site/tutorials/index.md`
- Create: `doc/site/models/index.md`
- Modify: `.gitignore` (append at end)

**Interfaces:**
- Produces: a Sphinx project that builds clean with execution off; `conf.py` reads `CLASSPP_DOCS_EXECUTE`; root toctree references `tutorials/index` and `models/index` (section toctrees still empty — filled by Task 3).

- [ ] **Step 1: Write `doc/site/requirements.txt`**

```
sphinx>=8.0
myst-nb>=1.2
pydata-sphinx-theme>=0.16
sphinx-design>=0.6
sphinx-copybutton>=0.5
ipykernel>=6.29
matplotlib>=3.8
scipy>=1.11
```

(`matplotlib`/`scipy`/`ipykernel` are notebook runtime deps for the executing build; `classy` itself is installed separately with `pip install .`.)

- [ ] **Step 2: Write `doc/site/conf.py`**

```python
"""Sphinx configuration for the CLASSpp documentation site.

Builds the human-facing site: landing page plus example notebooks.
The Doxygen C++ reference is built separately (doc/input/doxyconf) and
copied to <site>/api/ by the deploy workflow.
"""

import os

project = "CLASSpp"
author = "Aarhus Cosmology"
copyright = "2026, Aarhus Cosmology"  # noqa: A001

extensions = [
    "myst_nb",
    "sphinx_design",
    "sphinx_copybutton",
]

exclude_patterns = ["_build"]

myst_enable_extensions = ["colon_fence", "dollarmath", "amsmath"]

# -- Notebook execution ------------------------------------------------------
# Off by default so PR checks and local builds are fast and need no CLASS
# build.  The nightly deploy exports CLASSPP_DOCS_EXECUTE=force to run every
# notebook against current master.  Deliberately no jupyter-cache: it keys on
# notebook source, not the CLASS build, so it would serve stale plots after
# code changes.
nb_execution_mode = "force" if os.environ.get("CLASSPP_DOCS_EXECUTE") == "force" else "off"
nb_execution_timeout = 900
nb_execution_raise_on_error = True
nb_merge_streams = True

# -- HTML --------------------------------------------------------------------
html_theme = "pydata_sphinx_theme"
html_title = "CLASSpp"
html_theme_options = {
    "logo": {"text": "CLASSpp"},
    "github_url": "https://github.com/AarhusCosmology/CLASSpp_public",
    "external_links": [
        {"name": "C++ API", "url": "https://aarhuscosmology.github.io/CLASSpp/api/"},
    ],
    "navbar_align": "left",
    # sourcelink = "Show Source" → the .ipynb, i.e. the download link
    "secondary_sidebar_items": ["page-toc", "sourcelink"],
}
```

- [ ] **Step 3: Write `doc/site/index.md`**

`````markdown
# CLASSpp

**CLASSpp** is the C++ evolution of [CLASS](http://class-code.net), the Cosmic
Linear Anisotropy Solving System — an Einstein–Boltzmann solver computing CMB
and large-scale-structure observables for standard and exotic cosmologies.

The notebooks on this site are **executed nightly against CLASSpp master**, so
every plot reflects the current code.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} 🚀 Tutorials
:link: tutorials/index
:link-type: doc
Start here: compute your first spectra, then work through distances,
thermodynamics, perturbations and power spectra.
:::

:::{grid-item-card} 🌌 Cosmology models
:link: models/index
:link-type: doc
Showcases of extended models: axions, dark matter decaying to massive
daughters, self-interacting neutrinos, and more.
:::

:::{grid-item-card} 🛠️ C++ API reference
:link: api/index.html
Doxygen reference for the C++ code base: modules, the species-plugin
hierarchy, and source maps.
:::

:::{grid-item-card} 💻 Source code
:link: https://github.com/AarhusCosmology/CLASSpp_public
CLASSpp on GitHub: issues, releases and contributions welcome.
:::
::::

## Quickstart

Install the Python wrapper from PyPI:

```bash
pip install classy-community
```

Then compute your first CMB spectrum:

```python
from classy import Class

cosmo = Class()
cosmo.set({"output": "tCl", "l_max_scalars": 2500})
cosmo.compute()
cls = cosmo.raw_cl(2500)
```

Head to the [tutorials](tutorials/index) to go further, or build the C++
command-line code with `make class` (see the
[README](https://github.com/AarhusCosmology/CLASSpp_public#readme)).

```{toctree}
:hidden:

tutorials/index
models/index
```
`````

- [ ] **Step 4: Write `doc/site/tutorials/index.md`**

`````markdown
# Tutorials

Classic CLASS example notebooks, updated for CLASSpp and executed nightly
against master. Download any notebook via the "Show Source" link in its
right-hand sidebar.

Recommended order:

```{toctree}
:maxdepth: 1
```
`````

- [ ] **Step 5: Write `doc/site/models/index.md`**

`````markdown
# Cosmology models

Worked examples of extended cosmologies implemented in CLASSpp, executed
nightly against master.

```{toctree}
:maxdepth: 1
```
`````

- [ ] **Step 6: Append to `.gitignore`**

```
# Documentation site (doc/site): build output + build-time notebook copies
doc/site/_build/
doc/site/tutorials/*.ipynb
doc/site/models/*.ipynb
doc/site/models/dcdm_wdm_fig3_9models.py
```

- [ ] **Step 7: Create venv and install doc requirements**

Run:
```bash
python3 -m venv venv-docs
source venv-docs/bin/activate
pip install -r doc/site/requirements.txt
```
Expected: installs succeed (pure-python wheels).

- [ ] **Step 8: Build with execution off — verify it fails first without the files? No — infra task: verify the build passes**

Run:
```bash
source venv-docs/bin/activate
python -m sphinx -W --keep-going -b html doc/site doc/site/_build/html
```
Expected: exit 0, `build succeeded`. If `-W` trips on a benign theme warning, record the warning, fix it (do not blanket-suppress); `-W` must survive because the PR check uses it.

Then: open `doc/site/_build/html/index.html` and check the four cards render and the two section pages exist.

- [ ] **Step 9: Commit**

```bash
git add doc/site/requirements.txt doc/site/conf.py doc/site/index.md doc/site/tutorials/index.md doc/site/models/index.md .gitignore
git commit -m "docs-site: Sphinx skeleton (MyST-NB + pydata theme) for #355"
```

---

### Task 2: Add H1 title cells to curated notebooks

MyST-NB takes the page title from the notebook's first heading; 15 of the 17 curated notebooks have none.

**Files:**
- Modify: 15 files under `notebooks/` (list in the script below)
- Create (temporary, scratchpad — not committed): `add_titles.py`

**Interfaces:**
- Consumes: nothing from Task 1.
- Produces: every curated notebook starts with a markdown cell whose first line is `# <title>`. Task 3's toctrees rely on these titles rendering.

- [ ] **Step 1: Write the check (fails before, passes after)**

Save to scratchpad as `check_titles.py`:

```python
import json, sys

CURATED = [
    "warmup", "distances", "thermo", "cltt_terms", "one_k", "one_time",
    "many_times", "varying_neff", "neutrinohierarchy", "varying_pann",
    "Growth_with_w", "cl_ST", "370-axion-species",
    "dcdm-wdm-massive-decay-products", "173-Interacting_NCDM",
    "185-Hot_NEDE", "Thomas_PPF",
]
bad = []
for name in CURATED:
    nb = json.load(open(f"notebooks/{name}.ipynb"))
    first = nb["cells"][0]
    src = "".join(first["source"]).lstrip()
    if not (first["cell_type"] == "markdown" and src.startswith("# ")):
        bad.append(name)
print("missing H1:", bad)
sys.exit(1 if bad else 0)
```

Run: `python3 check_titles.py` — Expected: exit 1, lists 15 notebooks (all but `370-axion-species` and `dcdm-wdm-massive-decay-products`).

- [ ] **Step 2: Insert title cells**

Save to scratchpad as `add_titles.py` and run with the venv-docs python (has `nbformat` via myst-nb):

```python
import nbformat

TITLES = {
    "warmup": "Warmup: your first CLASS computation",
    "distances": "Cosmological distances",
    "thermo": "Thermodynamics history",
    "cltt_terms": "Anatomy of the CMB temperature spectrum",
    "one_k": "Perturbation evolution at one wavenumber",
    "one_time": "Perturbations at one moment in time",
    "many_times": "Perturbations across cosmic time",
    "varying_neff": "Varying the effective number of neutrinos",
    "neutrinohierarchy": "Neutrino mass hierarchies",
    "varying_pann": "Dark-matter annihilation and the CMB",
    "Growth_with_w": "Growth of structure with dynamical dark energy",
    "cl_ST": "Scalar and tensor CMB spectra",
    "173-Interacting_NCDM": "Self-interacting neutrinos",
    "185-Hot_NEDE": "Hot New Early Dark Energy",
    "Thomas_PPF": "Dark energy with the PPF formalism",
}

for name, title in TITLES.items():
    path = f"notebooks/{name}.ipynb"
    nb = nbformat.read(path, as_version=4)
    first = nb.cells[0]
    if first.cell_type == "markdown" and first.source.lstrip().startswith("# "):
        continue
    nb.cells.insert(0, nbformat.v4.new_markdown_cell(f"# {title}"))
    nbformat.write(nb, path)
    print("titled", name)
```

Run: `source venv-docs/bin/activate && python add_titles.py`
Expected: prints `titled <name>` for all 15.

- [ ] **Step 3: Re-run the check**

Run: `python3 check_titles.py` — Expected: `missing H1: []`, exit 0.

- [ ] **Step 4: Sanity-diff one notebook**

Run: `git diff notebooks/warmup.ipynb | head -40`
Expected: a single inserted markdown cell at the top; no other churn. If `nbformat.write` rewrote unrelated cell ids massively, note it but proceed (ids are cosmetic) — do NOT hand-edit JSON.

- [ ] **Step 5: Commit (explicit paths)**

```bash
git add notebooks/warmup.ipynb notebooks/distances.ipynb notebooks/thermo.ipynb notebooks/cltt_terms.ipynb notebooks/one_k.ipynb notebooks/one_time.ipynb notebooks/many_times.ipynb notebooks/varying_neff.ipynb notebooks/neutrinohierarchy.ipynb notebooks/varying_pann.ipynb notebooks/Growth_with_w.ipynb notebooks/cl_ST.ipynb notebooks/173-Interacting_NCDM.ipynb notebooks/185-Hot_NEDE.ipynb notebooks/Thomas_PPF.ipynb
git commit -m "notebooks: add H1 title cells for the documentation site"
```

---

### Task 3: Wire curated notebooks into the site

**Files:**
- Modify: `doc/site/conf.py` (add curation + copy step)
- Modify: `doc/site/tutorials/index.md` (fill toctree)
- Modify: `doc/site/models/index.md` (fill toctree)

**Interfaces:**
- Consumes: Task 1 conf.py; Task 2 title cells.
- Produces: `CURATED` dict in `conf.py` — shape `{"tutorials": [(slug, src_filename, [extra_files]), ...], "models": [...]}`. Task 4 edits this dict when dropping notebooks. Slugs (exact): tutorials `warmup, distances, thermo, cltt_terms, one_k, one_time, many_times, varying_neff, neutrinohierarchy, varying_pann, growth_with_w, cl_st`; models `axion_species, dcdm_wdm, interacting_ncdm, hot_nede, ppf_dark_energy`.

- [ ] **Step 1: Add curation block to `doc/site/conf.py`**

Insert between the `myst_enable_extensions` line and the `# -- Notebook execution` section:

```python
# -- Notebook curation -------------------------------------------------------
# Notebooks live in notebooks/ (single source of truth).  Each entry copies
# notebooks/<src> to <section>/<slug>.ipynb at build time; extra files are
# copied next to it (helper scripts the notebook executes).  To publish
# another notebook: add a line here and a slug to the section's index.md.

import shutil
from pathlib import Path

CURATED = {
    "tutorials": [
        ("warmup", "warmup.ipynb", []),
        ("distances", "distances.ipynb", []),
        ("thermo", "thermo.ipynb", []),
        ("cltt_terms", "cltt_terms.ipynb", []),
        ("one_k", "one_k.ipynb", []),
        ("one_time", "one_time.ipynb", []),
        ("many_times", "many_times.ipynb", []),
        ("varying_neff", "varying_neff.ipynb", []),
        ("neutrinohierarchy", "neutrinohierarchy.ipynb", []),
        ("varying_pann", "varying_pann.ipynb", []),
        ("growth_with_w", "Growth_with_w.ipynb", []),
        ("cl_st", "cl_ST.ipynb", []),
    ],
    "models": [
        ("axion_species", "370-axion-species.ipynb", []),
        ("dcdm_wdm", "dcdm-wdm-massive-decay-products.ipynb",
         ["dcdm_wdm_fig3_9models.py"]),
        ("interacting_ncdm", "173-Interacting_NCDM.ipynb", []),
        ("hot_nede", "185-Hot_NEDE.ipynb", []),
        ("ppf_dark_energy", "Thomas_PPF.ipynb", []),
    ],
}

_here = Path(__file__).resolve().parent
_notebooks_dir = _here.parent.parent / "notebooks"


def _copy_if_changed(src: Path, dest: Path) -> None:
    if dest.exists() and dest.read_bytes() == src.read_bytes():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dest)


for _section, _entries in CURATED.items():
    for _slug, _src, _extra in _entries:
        _copy_if_changed(_notebooks_dir / _src, _here / _section / f"{_slug}.ipynb")
        for _name in _extra:
            _copy_if_changed(_notebooks_dir / _name, _here / _section / _name)
```

- [ ] **Step 2: Fill `doc/site/tutorials/index.md` toctree**

Replace the empty toctree with:

```
```{toctree}
:maxdepth: 1

warmup
distances
thermo
cltt_terms
one_k
one_time
many_times
varying_neff
neutrinohierarchy
varying_pann
growth_with_w
cl_st
```
```

- [ ] **Step 3: Fill `doc/site/models/index.md` toctree**

```
```{toctree}
:maxdepth: 1

axion_species
dcdm_wdm
interacting_ncdm
hot_nede
ppf_dark_energy
```
```

- [ ] **Step 4: Build with execution off**

Run:
```bash
source venv-docs/bin/activate
rm -rf doc/site/_build
python -m sphinx -W --keep-going -b html doc/site doc/site/_build/html
```
Expected: exit 0; `doc/site/_build/html/tutorials/warmup.html` … `models/ppf_dark_energy.html` all exist (17 pages). Notebook-internal markdown may raise myst warnings under `-W` (e.g. header-level jumps); fix in `notebooks/` (they are our files) rather than suppressing, and fold those edits into this task's commit.

- [ ] **Step 5: Spot-check in browser**

Open `doc/site/_build/html/tutorials/warmup.html`: title "Warmup: your first CLASS computation", code cells visible (no outputs yet — execution off), "Show Source" link in right sidebar points at the `.ipynb`.

- [ ] **Step 6: Commit**

```bash
git add doc/site/conf.py doc/site/tutorials/index.md doc/site/models/index.md
git commit -m "docs-site: curation list + build-time notebook copy step"
```
(plus any `notebooks/*.ipynb` warning fixes from Step 4, staged explicitly.)

---

### Task 4: Execute all curated notebooks locally; fix or drop

The heart of the verification promise: every published notebook must run against current master.

**Files:**
- Modify: `notebooks/*.ipynb` (fixes), possibly `doc/site/conf.py` + section `index.md` (drops)

**Interfaces:**
- Consumes: `CURATED` dict (Task 3), `CLASSPP_DOCS_EXECUTE` (Task 1).
- Produces: a curation list where every entry executes clean; a triage table (kept in the commit message and later the PR body) recording fixed/dropped/clean per notebook.

- [ ] **Step 1: Build classy into the venv**

Run:
```bash
source venv-docs/bin/activate
pip install .
python -c "import classy; print(classy.__version__)"
```
Expected: build succeeds (scikit-build-core), version prints.

- [ ] **Step 2: Refresh copies, then execute notebooks one by one (survey pass)**

A per-notebook loop gives a full failure survey in one pass (sphinx would stop reporting usefully at the first hard failure):

```bash
source venv-docs/bin/activate
python -m sphinx -b html doc/site doc/site/_build/html -Q || true   # refresh copies
for nb in doc/site/tutorials/*.ipynb doc/site/models/*.ipynb; do
  echo "=== $nb"
  (cd "$(dirname "$nb")" && jupyter execute --timeout=900 "$(basename "$nb")") \
    && echo "OK $nb" || echo "FAIL $nb"
done 2>&1 | tee execute-survey.log
grep -E "^(OK|FAIL)" execute-survey.log
```

Expected: a mix of OK/FAIL. Known issue to expect: `interacting_ncdm` fails on `import classy_pivot`.

- [ ] **Step 3: Triage each FAIL — fix in `notebooks/`, or drop**

Rules:
- **Fix** when mechanical: dead imports (`import classy_pivot as classy` → `import classy`), renamed input parameters, removed classy methods with obvious modern equivalents, matplotlib API drift. Edit the source notebook under `notebooks/` (with nbformat or NotebookEdit — never the copies under `doc/site/`).
- **Drop** when the failure needs physics judgment, depends on non-repo data, or resists 20 minutes of effort: delete its `CURATED` line and toctree slug, keep the notebook in the repo, record the reason.
- Re-run Step 2's loop after each round until every remaining notebook prints OK.

- [ ] **Step 4: Full executing site build (the real pipeline)**

Run:
```bash
source venv-docs/bin/activate
rm -rf doc/site/_build
CLASSPP_DOCS_EXECUTE=force python -m sphinx -b html doc/site doc/site/_build/html
```
Expected: exit 0 (this is the exact nightly command; `-W` is deliberately absent here, matching the deploy job). Open `tutorials/warmup.html` and 2–3 others — plots present.

- [ ] **Step 5: Commit fixes and drops**

```bash
git add <each modified notebooks/*.ipynb> doc/site/conf.py doc/site/tutorials/index.md doc/site/models/index.md
git commit -m "notebooks: make curated set execute against master
<triage table: clean/fixed/dropped + reasons>"
```

---

### Task 5: Doxygen cross-link

**Files:**
- Modify: `doc/input/mainpage.md`

**Interfaces:**
- Consumes: site layout (`/api/` one level below `/`).
- Produces: Doxygen mainpage links back to the site root.

- [ ] **Step 1: Add banner link to `doc/input/mainpage.md`**

Insert after the first paragraph (after the line ending "…the exact class, function, or source-level contract."):

```markdown
> **Looking for tutorials and examples?** This page is the C++ API
> reference. Executed example notebooks and the getting-started guide live
> on the [CLASSpp documentation site](../index.html).
```

(The relative link works because Doxygen HTML is served from `/api/`.)

- [ ] **Step 2: Verify Doxygen still builds — if `doxygen` is installed locally**

Run: `command -v doxygen && (cd doc/input && doxygen doxyconf) | tail -2 || echo "doxygen not installed — CI's docs.yml covers this"`
Expected: either a successful doxygen run (check `doc/manual/html/index.html` contains the new link) or the skip message.

- [ ] **Step 3: Commit**

```bash
git add doc/input/mainpage.md
git commit -m "docs: link Doxygen mainpage back to the documentation site"
```

---

### Task 6: CI workflows

**Files:**
- Modify: `.github/workflows/docs.yml` (add site job)
- Modify: `.github/workflows/test_nightly.yml` (rewrite `deploy_documentation` job only)

**Interfaces:**
- Consumes: `doc/site/` project, `CLASSPP_DOCS_EXECUTE=force`, Doxygen output at `doc/manual/html`.
- Produces: Pages artifact rooted at `doc/site/_build/html` with `api/` inside it.

- [ ] **Step 1: Add `site` job to `docs.yml`**

Append to `.github/workflows/docs.yml` (keep the existing `doxygen` job untouched):

```yaml
  site:
    runs-on: [self-hosted, linux, light]

    steps:
      - uses: actions/checkout@v4

      - name: Create virtual Python environment
        run: |
          rm -rf venv-docs
          DEB_PYTHON_INSTALL_LAYOUT='deb' python -m virtualenv venv-docs --system-site-packages
          source venv-docs/bin/activate
          pip install -r doc/site/requirements.txt
          deactivate

      - name: Build site (no notebook execution)
        run: |
          source venv-docs/bin/activate
          python -m sphinx -W --keep-going -b html doc/site doc/site/_build/html
          deactivate

      - name: Remove virtual Python environment
        if: success() || failure()
        run: rm -rf venv-docs
```

- [ ] **Step 2: Rewrite `deploy_documentation` in `test_nightly.yml`**

Replace the existing job body (keep `name:`, `runs-on:`, `needs:`, `permissions:`, `concurrency:`, `environment:` exactly as they are) so the steps read:

```yaml
    steps:
      - uses: actions/checkout@v4

      - name: Build Doxygen HTML
        run: cd doc/input && doxygen doxyconf

      - name: Create virtual Python environment
        run: |
          rm -rf venv-docs
          DEB_PYTHON_INSTALL_LAYOUT='deb' python -m virtualenv venv-docs --system-site-packages
          source venv-docs/bin/activate
          pip install .
          pip install -r doc/site/requirements.txt
          deactivate

      - name: Build site with executed notebooks
        run: |
          source venv-docs/bin/activate
          CLASSPP_DOCS_EXECUTE=force python -m sphinx -b html doc/site doc/site/_build/html
          deactivate

      - name: Add C++ API reference
        run: cp -r doc/manual/html doc/site/_build/html/api

      - name: Configure GitHub Pages
        uses: actions/configure-pages@v5

      - name: Upload Pages artifact
        uses: actions/upload-pages-artifact@v4
        with:
          path: doc/site/_build/html

      - name: Deploy to GitHub Pages
        id: deployment
        uses: actions/deploy-pages@v4

      - name: Remove virtual Python environment
        if: success() || failure()
        run: rm -rf venv-docs
```

Note in the PR body: `pip install .` builds classy — if the `light` runner lacks a C++ toolchain, flip this job's `runs-on` to `[self-hosted, linux, heavy]` (Thomas owns the runners).

- [ ] **Step 3: Validate YAML**

Run:
```bash
python3 -c "import yaml; yaml.safe_load(open('.github/workflows/docs.yml')); yaml.safe_load(open('.github/workflows/test_nightly.yml')); print('yaml ok')"
command -v actionlint && actionlint .github/workflows/docs.yml .github/workflows/test_nightly.yml || true
```
Expected: `yaml ok`; actionlint clean if installed.

- [ ] **Step 4: Commit**

```bash
git add .github/workflows/docs.yml .github/workflows/test_nightly.yml
git commit -m "ci: build + deploy the documentation site with executed notebooks"
```

---

### Task 7: Final verification and PR

**Files:** none new.

- [ ] **Step 1: Clean full build replicating the nightly job**

Run:
```bash
source venv-docs/bin/activate
rm -rf doc/site/_build
CLASSPP_DOCS_EXECUTE=force python -m sphinx -b html doc/site doc/site/_build/html
mkdir -p doc/site/_build/html/api
[ -d doc/manual/html ] && cp -r doc/manual/html/. doc/site/_build/html/api/ || cp doc/site/_build/html/index.html doc/site/_build/html/api/index.html
```
Expected: exit 0. Browse the site root: cards → tutorials (with plots) → models → C++ API link resolves.

- [ ] **Step 2: Invoke superpowers:verification-before-completion, then superpowers:finishing-a-development-branch**

PR: base `master`, title `docs: notebook documentation site (landing + executed tutorials + /api/) (#355)`. Body: what moved where (Doxygen root → `/api/` — old deep links break), the triage table from Task 4, the `light`-runner toolchain caveat from Task 6, and `Fixes #355`. End the body with the standard generation footer from the global constraints.

---

## Self-Review (done at planning time)

- **Spec coverage:** layout (T1/T3), execution model (T1/T4), curation + copy step (T3), verification pass (T4), cross-links (T1 landing card + T5), CI (T6), gitignore (T1), download links (T1 sourcelink), out-of-scope items untouched. ✓
- **Placeholders:** none — all file contents and commands are concrete; Task 4 is inherently discovery-shaped but has explicit rules and expected failures. ✓
- **Type consistency:** `CURATED` 3-tuple shape and slugs identical in T3/T4; `CLASSPP_DOCS_EXECUTE=force` everywhere; venv name `venv-docs` everywhere. ✓
