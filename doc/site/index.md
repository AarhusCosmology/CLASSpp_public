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
