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

exclude_patterns = ["_build", "**.ipynb_checkpoints"]

myst_enable_extensions = ["colon_fence", "dollarmath", "amsmath"]

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
        # many_times: dropped until the array_hunt_growing_closeby failure at
        # z_max_pk=46000 is fixed (#380)
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
    _expected = set()
    for _slug, _src, _extra in _entries:
        _copy_if_changed(_notebooks_dir / _src, _here / _section / f"{_slug}.ipynb")
        _expected.add(f"{_slug}.ipynb")
        for _name in _extra:
            _copy_if_changed(_notebooks_dir / _name, _here / _section / _name)
            _expected.add(_name)
    # Remove everything else: stale copies of formerly-curated notebooks are
    # orphan pages that fail -W builds, and execution byproducts (savefig
    # PDFs, .npz caches) would go stale on reused self-hosted runners.
    for _stray in (_here / _section).iterdir():
        if _stray.is_file() and _stray.name not in _expected | {"index.md"}:
            _stray.unlink()

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
# Keep copied notebook sources as real .ipynb files (no ".txt" suffix) so the
# sidebar "Show Source" link doubles as a notebook download.
html_sourcelink_suffix = ""

html_theme = "pydata_sphinx_theme"
html_title = "CLASSpp"
html_theme_options = {
    "logo": {"text": "CLASSpp"},
    "github_url": "https://github.com/AarhusCosmology/CLASSpp_public",
    "external_links": [
        {"name": "C++ API", "url": "https://aarhuscosmology.github.io/CLASSpp/api/"},
    ],
    "navbar_align": "left",
    # sourcelink = "Show Source" -> the .ipynb, i.e. the download link
    "secondary_sidebar_items": ["page-toc", "sourcelink"],
}
