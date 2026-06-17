"""Characterization test for issue #319 (clean-loop background_indices()).

Locks the get_background() output (column names, order, values) against goldens
captured from the pre-refactor build. The refactor must keep every assertion
green: it relocates index registration and per-step writes, never the physics.

Exception: '(.)p_tot_prime' is checked at a tolerance, not exactly. The PPrime
refactor distributes the a*H factor across the per-species sum (Sum aH*d_i instead
of aH*Sum d_i) and folds the scalar field into the sum rather than adding it last;
both are ULP-level reorderings. Per project policy bit-identical is the wrong bar.
"""
import os
import numpy as np
import pytest

from gen_background_golden import CASES, background_for

GOLDEN_DIR = os.path.join(os.path.dirname(__file__), "background_golden")
TOL_COLUMNS = {"(.)p_tot_prime"}   # checked at rtol below, not exactly
PPRIME_RTOL = 1e-9


@pytest.mark.parametrize("name", sorted(CASES))
def test_background_columns_match_golden(name):
    golden = np.load(os.path.join(GOLDEN_DIR, name + ".npz"), allow_pickle=False)
    expected_order = list(golden["__order__"])

    bg = background_for(CASES[name])

    assert list(bg.keys()) == expected_order, (
        f"{name}: background column set/order changed.\n"
        f"  expected: {expected_order}\n  got:      {list(bg.keys())}"
    )

    for key in expected_order:
        if key in TOL_COLUMNS:
            np.testing.assert_allclose(
                bg[key], golden[key], rtol=PPRIME_RTOL, atol=0.0,
                err_msg=f"{name}: column {key!r} drifted beyond tolerance",
            )
        else:
            np.testing.assert_array_equal(
                bg[key], golden[key], err_msg=f"{name}: column {key!r} changed"
            )
