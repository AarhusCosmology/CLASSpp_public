"""Regression test for issue #391 (composites must name their children's columns).

d_tot is the rho-weighted sum of the per-species density contrasts, so it is
reconstructible from the d_i columns -- but only if every species that carries
density actually emits one. A composite that allocates and fills a child's
transfer slot without naming it makes the sum come up short by exactly that
child's weight in rho_tot, which is what this test catches.
"""
import numpy as np
import pytest

from gen_transfer_golden import COMMON, CASES
from classy import Class

# Species whose background rho column is not named `(.)rho_<suffix>` after the
# `d_<suffix>` transfer column. DR channels are the only irregular ones.
RHO_COLUMN_OVERRIDES = {"DR": "(.)rho_dr"}


def rho_and_transfer(extra):
    cosmo = Class()
    cosmo.set({**COMMON, **extra})
    try:
        cosmo.compute()
        return cosmo.get_transfer(), cosmo.get_background()
    finally:
        cosmo.struct_cleanup()
        cosmo.empty()


def test_dcdm_dr_delta_tot_is_reconstructible():
    """Every density-carrying species in the dcdm_dr sector emits a d_ column.

    Before #391 the DCDM_DR composite named only its DR daughter, so d_dcdm was
    computed every step and thrown away and this sum was short by the dcdm
    weight in rho_tot (7.8e-3).
    """
    tk, bg = rho_and_transfer(CASES["dcdm_dr"])

    z = bg["z"][::-1]

    def rho_today(column):
        return np.interp(0.0, z, bg[column][::-1])

    # CLASS follows the CMBFAST/CAMB convention of excluding Lambda from rho_tot.
    rho_tot = rho_today("(.)rho_tot") - rho_today("(.)rho_lambda")

    suffixes = [key[2:] for key in tk if key.startswith("d_") and key != "d_tot"]
    assert "d_dcdm" in tk, (
        "d_dcdm is missing: the DCDM_DR composite allocated and filled its DCDM "
        f"child's transfer slot but never named it. Columns: {sorted(tk)}"
    )

    recon = np.zeros_like(tk["d_tot"])
    for suffix in suffixes:
        column = RHO_COLUMN_OVERRIDES.get(suffix, f"(.)rho_{suffix}")
        assert column in bg, f"no background rho column {column!r} for d_{suffix}"
        recon += rho_today(column) * tk["d_" + suffix]
    recon /= rho_tot

    # Tolerance is set by interpolating the background rho to z=0, not by the
    # perturbation solve; the residual is ~1e-9 when no species is missing and
    # 7.8e-3 when dcdm is dropped, so the two regimes are orders apart.
    np.testing.assert_allclose(recon, tk["d_tot"], rtol=1e-6)
