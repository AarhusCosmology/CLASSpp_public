"""Characterization test for issue #309 (species-registered transfer sources).

Locks the Tk/vTk transfer-function output (column names, column order, and
values) against goldens captured from the pre-refactor build. The #309 refactor
must keep every assertion green: it changes only where source-slot integers are
stored, never the physics.
"""
import os
import numpy as np
import pytest

from gen_transfer_golden import CASES, transfer_for

GOLDEN_DIR = os.path.join(os.path.dirname(__file__), "transfer_golden")


@pytest.mark.parametrize("name", sorted(CASES))
def test_transfer_columns_match_golden(name):
    golden = np.load(os.path.join(GOLDEN_DIR, name + ".npz"), allow_pickle=False)
    expected_order = list(golden["__order__"])

    flat = transfer_for(CASES[name])

    assert list(flat.keys()) == expected_order, (
        f"{name}: transfer column set/order changed.\n"
        f"  expected: {expected_order}\n  got:      {list(flat.keys())}"
    )

    for key in expected_order:
        np.testing.assert_array_equal(
            flat[key], golden[key], err_msg=f"{name}: column {key!r} changed"
        )
