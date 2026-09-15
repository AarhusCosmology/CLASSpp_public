"""Drive CLASS++'s classy through Cobaya's own `classy` theory, end to end.

The probe likelihood below requests what the CMB, BAO and supernova likelihoods
request (lensed Cls up to l=2500, H(z), D_A(z), D_M(z), r_drag), so a failure here
is a failure of the Cobaya<->classy interface and needs no likelihood data.

Cobaya is an optional dependency, so this module skips without it -- unless
CLASSPP_REQUIRE_COBAYA=1, which CI sets so that a missing install fails loudly
instead of passing green without having run anything.
"""
import os

import numpy as np
import pytest

if os.environ.get("CLASSPP_REQUIRE_COBAYA") == "1":
    import cobaya  # noqa: F401  (ImportError here is the point)
else:
    pytest.importorskip("cobaya")

from cobaya.likelihood import Likelihood  # noqa: E402
from cobaya.model import get_model  # noqa: E402

Z = np.array([0.3, 0.51, 0.93, 1.32, 2.33])
L_MAX = 2500


class Probe(Likelihood):
    def get_requirements(self):
        return {
            "Cl": {"tt": L_MAX, "te": L_MAX, "ee": L_MAX, "pp": L_MAX},
            "Hubble": {"z": Z},
            "angular_diameter_distance": {"z": Z},
            "comoving_radial_distance": {"z": Z},
            "rdrag": None,
        }

    def logp(self, **params):
        cl = self.provider.get_Cl(ell_factor=True, units="muK2")
        return -0.5 * ((cl["tt"][220] - 5700.0) / 100.0) ** 2


def _model():
    return get_model({
        "likelihood": {"probe": Probe},
        "theory": {"classy": {"extra_args": {
            "N_ur": 2.0308,
            "nu1.type": "ncdm_standard",
        }}},
        "params": {
            # Cobaya's documented pattern for a CLASS name that is not a Python
            # identifier: sample a valid name, pass the CLASS one as an input
            # function, and drop the valid one so it is not sent to CLASS too.
            # CLASS++'s dot-notation keys (nu1.m, dncdm1.Gamma) need exactly this.
            "theta_s_100": {"prior": {"min": 1.03, "max": 1.05}, "ref": 1.04110,
                            "drop": True},
            "100*theta_s": {"value": "lambda theta_s_100: theta_s_100",
                            "derived": False},
            "m_nu": {"prior": {"min": 0.0, "max": 0.5}, "ref": 0.06, "drop": True},
            "nu1.m": {"value": "lambda m_nu: m_nu", "derived": False},
            "omega_b": 0.02237,
            "omega_cdm": 0.1200,
            "A_s": 2.1e-9,
            "n_s": 0.9649,
            "tau_reio": 0.0544,
            "H0": None,
            "sigma8": None,
            "Omega_nu": None,
        },
    })


@pytest.fixture(scope="module")
def model():
    return _model()


def _evaluate(model, m_nu):
    result = model.logposterior({"theta_s_100": 1.04110, "m_nu": m_nu})
    derived = dict(zip(model.parameterization.derived_params(), result.derived))
    return result, derived


def test_logposterior_is_finite(model):
    result, derived = _evaluate(model, 0.06)
    assert np.isfinite(result.logpost)
    # theta_s shooting with a massive neutrino lands on a Planck-like H0
    assert 66.0 < derived["H0"] < 69.0
    assert 0.7 < derived["sigma8"] < 0.9


def test_requested_products_reach_the_likelihood(model):
    _evaluate(model, 0.06)
    provider = model.provider
    cl = provider.get_Cl(ell_factor=True, units="muK2")
    for spectrum in ("tt", "te", "ee", "pp"):
        assert len(cl[spectrum]) == L_MAX + 1
        assert np.all(np.isfinite(cl[spectrum]))
    assert 5000.0 < cl["tt"][220] < 6500.0
    for getter in (provider.get_Hubble, provider.get_angular_diameter_distance,
                   provider.get_comoving_radial_distance):
        values = getter(Z)
        assert np.all(np.isfinite(values)) and np.all(values > 0)
    assert 140.0 < provider.get_param("rdrag") < 155.0


def test_dotted_input_parameter_reaches_class(model):
    # Omega_nu h^2 = m_nu / 93.14 eV: proves nu1.m was really set from m_nu.
    for m_nu in (0.06, 0.12):
        _, derived = _evaluate(model, m_nu)
        omega_nu = derived["Omega_nu"] * (derived["H0"] / 100.0) ** 2
        assert omega_nu == pytest.approx(m_nu / 93.14, rel=0.02)
