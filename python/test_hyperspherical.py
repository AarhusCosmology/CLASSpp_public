import numpy as np
import pytest

# scipy provides the ground-truth spherical Bessel; skip the whole module if it
# is unavailable rather than failing at collection (CI installs only pytest).
spherical_jn = pytest.importorskip("scipy.special").spherical_jn
from classy import hyperspherical_bessel_direct


def test_direct_flat_matches_spherical_jn():
    """Flat case: Phi_l^beta(x) = j_l(beta * x)."""
    beta = 30.0
    x = np.linspace(0.05, 2.0, 50)
    for l in (0, 1, 2, 5, 10, 20):
        got = hyperspherical_bessel_direct(0, beta, l, x)
        ref = spherical_jn(l, beta * x)
        assert got.shape == (x.size,)
        assert np.allclose(got, ref, rtol=1e-6, atol=1e-9), f"l={l}"


def test_direct_array_l_shape():
    beta = 30.0
    x = np.linspace(0.05, 2.0, 40)
    l = np.array([0, 3, 7])
    got = hyperspherical_bessel_direct(0, beta, l, x)
    assert got.shape == (3, 40)
    for i, ll in enumerate(l):
        assert np.allclose(got[i], spherical_jn(ll, beta * x), rtol=1e-6, atol=1e-9)


from classy import hyperspherical_bessel_interpolate


def test_interpolate_flat_matches_spherical_jn():
    beta = 50.0
    x = np.linspace(0.1, 2.0, 60)
    for l in (0, 2, 8, 25):
        got = hyperspherical_bessel_interpolate(0, beta, l, x)
        ref = spherical_jn(l, beta * x)
        assert np.allclose(got, ref, rtol=1e-4, atol=1e-7), f"l={l}"


def test_interpolate_matches_direct_flat():
    beta = 50.0
    x = np.linspace(0.1, 2.0, 60)
    l = np.array([1, 4, 12])
    interp = hyperspherical_bessel_interpolate(0, beta, l, x)
    direct = hyperspherical_bessel_direct(0, beta, l, x)
    assert np.allclose(interp, direct, rtol=1e-4, atol=1e-7)


def test_interpolate_derivatives_shapes():
    beta = 50.0
    x = np.linspace(0.1, 2.0, 30)
    l = np.array([2, 5])
    out = hyperspherical_bessel_interpolate(0, beta, l, x, derivatives=True)
    assert isinstance(out, tuple) and len(out) == 3
    for arr in out:
        assert arr.shape == (2, 30)
    # First derivative of j_l(beta x) is beta * j_l'(beta x).
    Phi, dPhi, d2Phi = out
    ref0 = beta * spherical_jn(2, beta * x, derivative=True)
    assert np.allclose(dPhi[0], ref0, rtol=1e-3, atol=1e-6)


from classy import hyperspherical_bessel_wkb


def test_wkb_matches_direct_open_large_l():
    """Open case: WKB approximates the recurrence well at large l, away from
    the turning point."""
    K, beta = -1, 200.0
    l = 80
    # Turning point asinh(sqrt(l(l+1))/beta); sample in the oscillatory region.
    x = np.linspace(1.0, 3.0, 40)
    wkb = hyperspherical_bessel_wkb(K, beta, l, x)
    direct = hyperspherical_bessel_direct(K, beta, l, x)
    # WKB is asymptotic; compare on a robust normalized error.
    scale = np.maximum(np.abs(direct), 1e-3 * np.abs(direct).max())
    assert np.median(np.abs(wkb - direct) / scale) < 0.05


def test_wkb_rejects_flat():
    with pytest.raises(ValueError):
        hyperspherical_bessel_wkb(0, 100.0, 5, np.array([0.5, 1.0]))


def test_wkb_rejects_l_zero():
    with pytest.raises(ValueError):
        hyperspherical_bessel_wkb(-1, 100.0, 0, np.array([0.5, 1.0]))


def test_closed_case_symmetry():
    """Closed case: direct evaluation beyond pi/2 equals the symmetry-folded
    value. Phi_l(pi - y) = (-1)^(beta-l-1) Phi_l(y)."""
    K = 1
    beta = 40  # integer
    l = 5
    y = np.linspace(0.05, 0.5 * np.pi - 0.05, 30)
    base = hyperspherical_bessel_direct(K, beta, l, y)
    reflected = hyperspherical_bessel_direct(K, beta, l, np.pi - y)
    sign = (-1.0) ** (beta - l - 1)
    assert np.allclose(reflected, sign * base, rtol=1e-6, atol=1e-9)


def test_closed_case_interpolate_matches_direct_beyond_halfpi():
    K, beta, l = 1, 40, 4
    x = np.linspace(0.6 * np.pi, 1.3 * np.pi, 25)  # outside [0, pi/2]
    interp = hyperspherical_bessel_interpolate(K, beta, l, x)
    direct = hyperspherical_bessel_direct(K, beta, l, x)
    assert np.allclose(interp, direct, rtol=1e-4, atol=1e-7)


def test_direct_continuous_across_turning_point():
    """Open case: direct is continuous across xfwd (the bwd->fwd switch)."""
    K, beta, l = -1, 60.0, 30
    xfwd = np.arcsinh(np.sqrt(l * (l + 1.0)) / beta)
    x = np.linspace(xfwd - 0.02, xfwd + 0.02, 41)
    got = hyperspherical_bessel_direct(K, beta, l, x)
    # No jump: max adjacent difference is small relative to local scale.
    jumps = np.abs(np.diff(got))
    assert jumps.max() < 5e-3, jumps.max()
    # And it matches interpolation through the seam.
    interp = hyperspherical_bessel_interpolate(K, beta, l, x)
    assert np.allclose(got, interp, rtol=1e-3, atol=1e-6)


def test_scalar_l_returns_1d():
    out = hyperspherical_bessel_direct(0, 30.0, 3, np.linspace(0.1, 1.0, 10))
    assert out.ndim == 1 and out.shape == (10,)


def test_closed_case_rejects_noninteger_beta():
    with pytest.raises(ValueError):
        hyperspherical_bessel_direct(1, 40.5, 3, np.array([0.5, 1.0]))


def test_closed_case_rejects_l_ge_beta():
    with pytest.raises(ValueError):
        hyperspherical_bessel_direct(1, 10, 10, np.array([0.5, 1.0]))


def test_invalid_K_rejected():
    with pytest.raises(ValueError):
        hyperspherical_bessel_direct(2, 30.0, 3, np.array([0.5]))


def test_derivatives_only_on_interpolate():
    with pytest.raises(TypeError):
        hyperspherical_bessel_direct(0, 30.0, 3, np.array([0.5]), derivatives=True)


def test_interpolate_preserves_unordered_duplicate_l():
    """interpolate sorts l internally (unique) but must return rows in the
    requested order, including duplicates."""
    beta = 50.0
    x = np.linspace(0.1, 2.0, 30)
    l = np.array([12, 4, 1, 4])  # unsorted, with a duplicate
    out = hyperspherical_bessel_interpolate(0, beta, l, x)
    assert out.shape == (4, 30)
    for i, ll in enumerate(l):
        assert np.allclose(out[i], spherical_jn(ll, beta * x), rtol=1e-4, atol=1e-7), f"row {i} (l={ll})"
    # the two l=4 rows must be identical
    assert np.array_equal(out[1], out[3])


def test_interpolate_and_wkb_scalar_l_return_1d():
    x = np.linspace(0.1, 1.0, 12)
    assert hyperspherical_bessel_interpolate(0, 50.0, 3, x).shape == (12,)
    assert hyperspherical_bessel_wkb(-1, 100.0, 5, x).shape == (12,)


def test_empty_l_rejected():
    with pytest.raises(ValueError):
        hyperspherical_bessel_direct(0, 30.0, np.array([], dtype=int), np.array([0.5]))


def test_empty_x_rejected():
    with pytest.raises(ValueError):
        hyperspherical_bessel_direct(0, 30.0, 3, np.array([]))


def test_interpolate_small_l_no_corruption():
    """Regression: small l makes the classical turning point xfwd fall below
    xmin, which used to drive xfwdidx negative inside HIS_create and read/write
    pHIS->x/phi/dphi out of bounds (heap corruption -> intermittent NaN). Run
    interpolate for the smallest l repeatedly and confirm finite, correct
    values every time."""
    beta = 50.0
    x = np.linspace(0.1, 2.0, 60)
    for _ in range(50):
        for l in (0, 1, 2):
            got = hyperspherical_bessel_interpolate(0, beta, l, x)
            assert np.all(np.isfinite(got)), f"non-finite at l={l}"
            assert np.allclose(got, spherical_jn(l, beta * x), rtol=1e-4, atol=1e-7), f"l={l}"
