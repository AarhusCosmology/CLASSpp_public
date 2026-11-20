"""Hyperspherical Bessel functions Phi_l^beta(x), as the transfer module computes them.

The input checks and the choice of grid are here; the compiled core makes the C++ calls.
"""
import numpy as np

from . import _core


def _normalize_bessel_inputs(K, beta, l, x, method):
    """Validate and coerce (K, beta, l, x). Returns (lvec_int32, xvec_float64,
    scalar_l_bool). Raises ValueError on invalid input. method in
    {'interpolate','direct','wkb'} for method-specific checks."""
    if K not in (-1, 0, 1):
        raise ValueError("K must be -1, 0, or 1 (got %r)" % (K,))
    beta = float(beta)
    if beta <= 0.0:
        raise ValueError("beta must be positive (got %r)" % (beta,))
    scalar_l = np.isscalar(l)
    lvec = np.ascontiguousarray(np.atleast_1d(l), dtype=np.intc)
    if lvec.size == 0:
        raise ValueError("l must contain at least one value")
    if np.any(lvec < 0):
        raise ValueError("all l must be non-negative")
    xvec = np.ascontiguousarray(np.atleast_1d(x), dtype=np.double)
    if xvec.size == 0:
        raise ValueError("x must contain at least one value")
    if np.any(xvec <= 0.0):
        raise ValueError("all x must be positive")
    if K == 1:
        if abs(beta - round(beta)) > 1e-6:
            raise ValueError("closed case (K=1) requires integer beta (got %r)" % (beta,))
        if int(lvec.max()) >= int(round(beta)):
            raise ValueError("closed case requires l < beta")
    if method == 'wkb':
        if K == 0:
            raise ValueError("wkb is only defined for curved space (K = +/-1)")
        if int(lvec.min()) < 1:
            raise ValueError("wkb requires l >= 1")
    return lvec, xvec, scalar_l


def hyperspherical_bessel_direct(K, beta, l, x):
    """Hyperspherical Bessel function Phi_l^beta(x) by direct recurrence.

    Parameters
    ----------
    K : int          curvature sign, -1 (open), 0 (flat), or 1 (closed)
    beta : float     wavenumber nu (> 0; integer for K=1)
    l : int or array of int   multipole(s) (>= 0; < beta for K=1)
    x : array of float        evaluation points (> 0)

    Returns a float64 array: shape (n_x,) for scalar l, else (n_l, n_x).
    Closed-case symmetry and turning-point stability are handled internally.
    """
    lvec, xvec, scalar_l = _normalize_bessel_inputs(K, beta, l, x, 'direct')
    Phi = _core._bessel_direct(K, float(beta), lvec, xvec)
    if scalar_l:
        return np.ascontiguousarray(Phi[0])
    return Phi


def hyperspherical_bessel_interpolate(K, beta, l, x, sampling=None, derivatives=False):
    """Hyperspherical Bessel function Phi_l^beta(x) by grid + Hermite interpolation
    (the path used internally by the transfer module).

    Parameters as in hyperspherical_bessel_direct, plus:
    sampling : float or None    grid points per wavelength (None -> CLASS default:
                                8.0 flat; 7.0 for nu<1000, 3.0 otherwise).
    derivatives : bool          if True return (Phi, dPhi, d2Phi).

    Returns a float64 array (or a 3-tuple of them if derivatives=True): shape
    (n_x,) for scalar l, else (n_l, n_x).

    Note: for very narrow x ranges (span < ~1 wavelength) accuracy near the grid
    edges degrades; use hyperspherical_bessel_direct for single-point evaluation.
    """
    beta = float(beta)
    lvec_arr, xvec_arr, scalar_l = _normalize_bessel_inputs(K, beta, l, x, 'interpolate')

    # Sorted unique l for the interpolation structure; remember each request's index.
    uniq = np.unique(lvec_arr)
    luniq = np.ascontiguousarray(uniq, dtype=np.intc)

    if sampling is None:
        samp = 8.0 if K == 0 else (7.0 if beta < 1000.0 else 3.0)
    else:
        samp = float(sampling)

    pi = np.pi
    if K == 1:
        # Closed case: requested x are folded by ClosedModY into the fundamental
        # domain, so the grid must always span all of [eps, pi/2 - eps]
        # regardless of the requested x range.
        xmin = 1.0e-5
        xmax = 0.5 * pi - 1.0e-5
    else:
        # The Hermite interpolation returns 0 for any point that lands strictly
        # outside [grid_xmin, grid_xmax]; because the grid endpoints suffer
        # floating-point rounding (xmin + (nx-1)*deltax) a requested point at the
        # exact min/max can fall just outside and read back as 0. Also, if all
        # requested x are equal, xmin == xmax yields a zero-width grid (deltax=0)
        # and undefined behaviour. Pad the grid outward so every requested point
        # is strictly interior, with a minimum span tied to the wavelength.
        req_min = float(xvec_arr.min())
        req_max = float(xvec_arr.max())
        # Grid step the C routine will use is lambda/sampling = 2*pi/(beta*samp).
        # The interpolation is only reliable strictly inside the grid: the first
        # and last cells couple to boundary derivatives that are noisy (and the
        # closed/forward recurrence leaves the very last cell unstable). Pad the
        # grid by a few cells on each side so that every requested point lands in
        # an interior cell, and guarantee a non-degenerate span (handles the
        # all-x-equal / single-point case where req_min == req_max).
        grid_step = 2.0 * pi / (beta * samp)
        # second term: fallback for pathologically large beta where grid_step
        # underflows to 0 (req_max is always > 0; _normalize rejects x <= 0).
        pad = max(3.0 * grid_step, 1.0e-6 * req_max)
        xmin = max(req_min - pad, 1.0e-5)
        xmax = req_max + pad
        # Guarantee enough grid cells that the requested span is comfortably
        # interior even for a single point / very narrow request (otherwise the
        # grid collapses to a handful of boundary-dominated cells). Grow the
        # high side until the grid spans at least ~16 cells.
        min_span = 16.0 * grid_step
        if xmax - xmin < min_span:
            xmax = xmin + min_span
    l_WKB = int(luniq[-1]) + 1
    phiminabs = 1.0e-10

    Phi, dPhi, d2Phi = _core._bessel_interpolate(K, beta, luniq, xvec_arr, xmin, xmax, samp,
                                                 l_WKB, phiminabs, derivatives)

    # Map sorted-unique rows back to the requested l order.
    idx = np.searchsorted(uniq, lvec_arr)

    def _shape(M):
        out = M[idx]
        return np.ascontiguousarray(out[0]) if scalar_l else out

    if derivatives:
        return _shape(Phi), _shape(dPhi), _shape(d2Phi)
    return _shape(Phi)


def hyperspherical_bessel_wkb(K, beta, l, x):
    """Hyperspherical Bessel function Phi_l^beta(x) by WKB (Airy) approximation.

    Curved space only (K = +/-1) and l >= 1. Parameters as in
    hyperspherical_bessel_direct. Returns a float64 array: shape (n_x,) for
    scalar l, else (n_l, n_x). Closed-case symmetry is handled internally.
    """
    lvec, xvec, scalar_l = _normalize_bessel_inputs(K, beta, l, x, 'wkb')
    Phi = _core._bessel_wkb(K, float(beta), lvec, xvec)
    if scalar_l:
        return np.ascontiguousarray(Phi[0])
    return Phi
