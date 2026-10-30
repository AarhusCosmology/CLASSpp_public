# distutils: language = c++
"""
.. module:: classy
    :synopsis: Python wrapper around CLASS
.. moduleauthor:: Karim Benabed <benabed@iap.fr>
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
.. moduleauthor:: Julien Lesgourgues <lesgourg@cern.ch>
.. moduleauthor:: Thomas Tram <thomas.tram@phys.au.dk>

This module defines a class called Class. It is used with Monte Python to
extract cosmological parameters.

"""
from libc.math cimport exp, log
from libc.stdlib cimport malloc, calloc, free
from libc.stdio cimport snprintf

from libcpp cimport bool
from libcpp.memory cimport nullptr
from libcpp.memory cimport unique_ptr
from libcpp.memory cimport shared_ptr
from libcpp.memory cimport make_shared
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

from numpy.math cimport EULER, LOGE2

import cython
cimport cython
from cython.operator cimport dereference as deref

import numpy as np
cimport numpy as np
import os
import sys

# Implement a specific Exception (this might not be optimally designed, nor
# even acceptable for python standards. It, however, does the job).
# The idea is to raise either an AttributeError if the problem happened while
# reading the parameters (in the normal Class, this would just return a line in
# the unused_parameters file), or a NameError in other cases. This allows
# MontePython to handle things differently.
class CosmoError(Exception):
    def __init__(self, message=None):
        if message is None:
            message = ''
        elif isinstance(message, bytes):
            try:
                # If message has garbage inside, decode may fail.
                message = message.decode()
            except:
                pass

        self.message = message

    def __str__(self):
        return '\n\nError in Class: ' + self.message


class CosmoSevereError(CosmoError):
    """
    Raised when Class failed to understand one or more input parameters.

    This case would not raise any problem in Class default behaviour. However,
    for parameter extraction, one has to be sure that all input parameters were
    understood, otherwise the wrong cosmological model would be selected.
    """
    pass


class CosmoComputationError(CosmoError):
    """
    Raised when Class could not compute the cosmology at this point.

    This will be caught by the parameter extraction code to give an extremely
    unlikely value to this point
    """

cdef int raise_my_py_error() except *:
    cdef:
        pair[string, string] cpp_exception
        str cpp_exception_type
    cpp_exception = get_my_py_error_message()
    cpp_exception_type = cpp_exception.first.decode()
    if "invalid_argument" in cpp_exception_type:
        raise CosmoSevereError(cpp_exception.second)
    elif "runtime_error" in cpp_exception_type:
        raise CosmoComputationError(cpp_exception.second)
    else:
        # logic_error (programmer-error invariants) and anything unrecognized:
        # abort loudly rather than surfacing a misleading NotImplementedError.
        raise CosmoSevereError(cpp_exception.second)

# The declarations generate_wrapper.py harvests from the C++ headers. Every
# function among them carries `except +raise_my_py_error`, and Cython looks an
# exception handler up where the function is declared, so they are included
# into this module after the handler, not cimported from a .pxd of their own.
include "cclassy.pxi"

cdef extern from "version.h":
    const char* _VERSION_

# The version CMake stamps into version.h from pyproject.toml, so the module and
# the wheel metadata cannot disagree. Cobaya refuses a classy without it.
__version__ = _VERSION_.decode("utf-8")

cdef ClassConstants constvals

# Derived parameters of the form "<species key>.<field>", served by
# BackgroundModule::GetSpeciesParam (the species answers through its GetParam
# override, or generically from the background table). GetSpeciesParam reports
# an unanswered name as 0.0, which is indistinguishable from a value, so a
# request is refused unless all three of these hold: the field is listed below,
# HasSpecies() knows the species key, and HasSpeciesParam() says that species
# answers that field.
#   f_peak, a_peak, z_peak  -- any species: peak of rho/rho_crit(a) and where it
#                              sits. For an early-dark-energy species (the
#                              pheno_axion fluid, or an axion scalar field) this
#                              is f_EDE and its peak location.
#   a_c, n_axion, m_fld, alpha_fld, omega_axion -- pheno-axion fluid scales.
#   alpha, q0, x, M_2, M_3, M_4                 -- greybody NCDM shape/moments.
#   z_star                                      -- redshift of the NEDE transition.
#   G_eff_today, gamma_PPN_minus_1              -- scalar-tensor: Cavendish Newton constant
#                                                  over the action's, and gamma_PPN - 1,
#                                                  today, for a massless field.
# Derived distances also offered in units of Mpc/h, under the name with '_h' appended.
_H_SCALED_DISTANCES = frozenset(('rs_rec_h', 'ds_rec_h', 'ra_rec_h', 'da_rec_h', 'ds_d_h', 'rs_d_h'))

_SPECIES_DERIVED_FIELDS = frozenset((
    'f_peak', 'a_peak', 'z_peak',
    'a_c', 'n_axion', 'm_fld', 'alpha_fld', 'omega_axion',
    'alpha', 'q0', 'x', 'M_2', 'M_3', 'M_4',
    'z_star',
    'G_eff_today', 'gamma_PPN_minus_1',
))

# Nils : Added for python 3.x and python 2.x compatibility
cdef viewdictitems(dict d):
    if sys.version_info >= (3,0):
        return d.items()
    else:
        return d.viewitems()

cdef extern from "bessel_table_cache.h":
    cdef cppclass BesselTableCache:
        BesselTableCache()

cdef extern from "cosmology.h":
    ctypedef shared_ptr[const InputModule] InputModulePtr
    ctypedef shared_ptr[const BackgroundModule] BackgroundModulePtr
    ctypedef shared_ptr[const ThermodynamicsModule] ThermodynamicsModulePtr
    ctypedef shared_ptr[const PerturbationsModule] PerturbationsModulePtr
    ctypedef shared_ptr[const PrimordialModule] PrimordialModulePtr
    ctypedef shared_ptr[const NonlinearModule] NonlinearModulePtr
    ctypedef shared_ptr[const TransferModule] TransferModulePtr
    ctypedef shared_ptr[const SpectraModule] SpectraModulePtr
    ctypedef shared_ptr[const LensingModule] LensingModulePtr

    cdef cppclass Cosmology:
        Cosmology(FileContent& fc, shared_ptr[BesselTableCache] bessel_cache) except +raise_my_py_error
        InputModulePtr& GetInputModule()
        BackgroundModulePtr& GetBackgroundModule() except +raise_my_py_error
        ThermodynamicsModulePtr& GetThermodynamicsModule() except +raise_my_py_error
        PerturbationsModulePtr& GetPerturbationsModule() except +raise_my_py_error
        PrimordialModulePtr& GetPrimordialModule() except +raise_my_py_error
        NonlinearModulePtr& GetNonlinearModule() except +raise_my_py_error
        TransferModulePtr& GetTransferModule() except +raise_my_py_error
        SpectraModulePtr& GetSpectraModule() except +raise_my_py_error
        LensingModulePtr& GetLensingModule() except +raise_my_py_error

cdef extern from "hyperspherical.h":
    cdef cppclass HyperInterpStruct:
        HyperInterpStruct(int K, double beta, int nl, int* lvec, double xmin,
                          double xmax, double sampling, int l_WKB, double phiminabs) except +
    void hyperspherical_Hermite_interpolation_all "hyperspherical_Hermite_interpolation<6, true, true, true>"(
        const HyperInterpStruct* pHIS, int nxi, int lnum, const double* xinterp,
        double* Phi, double* dPhi, double* d2Phi) except +
    void hyperspherical_Hermite_interpolation_Phi "hyperspherical_Hermite_interpolation<6, true, false, false>"(
        const HyperInterpStruct* pHIS, int nxi, int lnum, const double* xinterp,
        double* Phi, double* dPhi, double* d2Phi) except +
    void hyperspherical_bessel_direct_vector(int K, double beta, int* lvec, int nl,
                                             double* xvec, int nx, double* Phi) except +
    bint hyperspherical_WKB(int K, int l, double beta, double y, double* Phi) except +


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
    cdef int Kc = K
    cdef double betac = float(beta)
    cdef int[::1] lvec
    cdef double[::1] xvec
    lvec_arr, xvec_arr, scalar_l = _normalize_bessel_inputs(K, beta, l, x, 'direct')
    lvec = lvec_arr
    xvec = xvec_arr
    cdef int nl = lvec_arr.shape[0]
    cdef int nx = xvec_arr.shape[0]
    cdef np.ndarray[double, ndim=2] Phi = np.empty((nl, nx), dtype=np.double)
    hyperspherical_bessel_direct_vector(Kc, betac, &lvec[0], nl, &xvec[0], nx, &Phi[0, 0])
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
    cdef int Kc = K
    cdef double betac = float(beta)
    lvec_arr, xvec_arr, scalar_l = _normalize_bessel_inputs(K, beta, l, x, 'interpolate')
    cdef int nx = xvec_arr.shape[0]

    # Sorted unique l for the interpolation structure; remember each request's index.
    uniq = np.unique(lvec_arr)
    cdef int[::1] luniq = np.ascontiguousarray(uniq, dtype=np.intc)
    cdef int nl = luniq.shape[0]
    cdef double[::1] xv = xvec_arr

    if sampling is None:
        samp = 8.0 if K == 0 else (7.0 if betac < 1000.0 else 3.0)
    else:
        samp = float(sampling)

    cdef double pi = np.pi
    cdef double xmin, xmax
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
        grid_step = 2.0 * pi / (betac * samp)
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
    cdef int l_WKB = int(luniq[nl - 1]) + 1
    cdef double phiminabs = 1.0e-10

    # Heap-allocate with `new` so the C++ constructor runs and fills the
    # struct. A stack `cdef HyperInterpStruct his` would bypass the constructor
    # (Cython zero-initializes the raw bytes), leaving the std::vector members
    # in an invalid state and producing sporadic NaNs/heap corruption.
    cdef HyperInterpStruct* his = new HyperInterpStruct(Kc, betac, nl, &luniq[0], xmin, xmax,
                                                        samp, l_WKB, phiminabs)
    cdef np.ndarray[double, ndim=2] Phi = np.empty((nl, nx), dtype=np.double)
    cdef np.ndarray[double, ndim=2] dPhi
    cdef np.ndarray[double, ndim=2] d2Phi
    cdef int i
    try:
        if derivatives:
            dPhi = np.empty((nl, nx), dtype=np.double)
            d2Phi = np.empty((nl, nx), dtype=np.double)
            for i in range(nl):
                hyperspherical_Hermite_interpolation_all(his, nx, i, &xv[0],
                                                         &Phi[i, 0], &dPhi[i, 0], &d2Phi[i, 0])
        else:
            for i in range(nl):
                hyperspherical_Hermite_interpolation_Phi(his, nx, i, &xv[0],
                                                         &Phi[i, 0], NULL, NULL)
    finally:
        del his

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
    cdef int Kc = K
    cdef double betac = float(beta)
    lvec_arr, xvec_arr, scalar_l = _normalize_bessel_inputs(K, beta, l, x, 'wkb')
    cdef int[::1] lvec = lvec_arr
    cdef double[::1] xvec = xvec_arr
    cdef int nl = lvec_arr.shape[0]
    cdef int nx = xvec_arr.shape[0]
    cdef np.ndarray[double, ndim=2] Phi = np.empty((nl, nx), dtype=np.double)
    cdef int il, ix
    cdef double val
    for il in range(nl):
        for ix in range(nx):
            if not hyperspherical_WKB(Kc, lvec[il], betac, xvec[ix], &val):
                raise CosmoSevereError("hyperspherical_WKB failed")
            Phi[il, ix] = val
    if scalar_l:
        return np.ascontiguousarray(Phi[0])
    return Phi


# To support the legacy name Class for the cosmology class.
cdef class Class(PyCosmology):
    pass

cdef class PyCosmology:

    """
    Cython wrapper class for C++ class Cosmology
    """

    cdef unique_ptr[Cosmology] _thisptr
    # The transfer module's Bessel table, kept across the Cosmology objects this wrapper builds
    # one after another (include/bessel_table_cache.h). A cached table never changes a result.
    cdef shared_ptr[BesselTableCache] _bessel_cache
    cdef dict _pars
    cdef bool parameters_changed
    cdef FileContent _fc

    property pars:
        def __get__(self):
            # Snapshot, not the live dict: a caller mutating the returned
            # dict would otherwise bypass the parameters_changed gate and
            # stale-serve. set() is the mutation API.
            return dict(self._pars)
    property state:
        def __get__(self):
            return True
    property Omega_nu:
        def __get__(self):
            return self._background().GetOmega0NcdmTot()
    property nonlinear_method:
        def __get__(self):
            return self.nl().method

    def __init__(self, input_parameters=None):
        if input_parameters is None:
            input_parameters = {}
        # Copy: never alias the caller's dict, so external mutation cannot
        # bypass the parameters_changed gate (and the class_dir default does
        # not leak into the caller's dict).
        self._pars = dict(input_parameters)
        if 'class_dir' not in self._pars:
            self._pars['class_dir'] = os.path.join(os.path.dirname(__file__), 'classy')
        self.reset()

    cdef reset(self):
        self._update_fc_from_pars()
        # Cosmology(FileContent&) rejects parameters that nothing read
        # (InputModule::RejectUnreadParameters), as a CosmoSevereError, unless
        # the input sets allow_unread_parameters = yes.
        if self._bessel_cache.get() == NULL:
            self._bessel_cache = make_shared[BesselTableCache]()
        self._thisptr.reset(new Cosmology(self._fc, self._bessel_cache))
        # Only clear the dirty flag after all throwing operations have
        # succeeded. If Cosmology construction throws, parameters_changed
        # stays True so the next _cosmo() call re-runs reset() and re-raises
        # instead of serving stale results.
        self.parameters_changed = False
        return self

    # Single gate to the C++ object: rebuilds the Cosmology if the parameter
    # dict changed since it was last built. Everything must reach the C++
    # object through _cosmo() (or the struct accessors below), never through
    # _thisptr directly, so results always reflect the current parameters.
    cdef Cosmology* _cosmo(self) except NULL:
        if self.parameters_changed or self._thisptr.get() == NULL:
            self.reset()
        return self._thisptr.get()

    # The input structs live inside the InputModule, and DoShooting replaces
    # the InputModule instance on first use; fetching through GetInputModule()
    # on every access is the only way a cached pointer cannot go stale.
    cdef const precision* pr(self) except NULL:
        return &deref(deref(self._cosmo()).GetInputModule()).precision_

    cdef const background* ba(self) except NULL:
        return &deref(deref(self._cosmo()).GetInputModule()).background_

    cdef const thermo* th(self) except NULL:
        return &deref(deref(self._cosmo()).GetInputModule()).thermodynamics_

    cdef const perturbs* pt(self) except NULL:
        return &deref(deref(self._cosmo()).GetInputModule()).perturbations_

    cdef const primordial* pm(self) except NULL:
        return &deref(deref(self._cosmo()).GetInputModule()).primordial_

    cdef const nonlinear* nl(self) except NULL:
        return &deref(deref(self._cosmo()).GetInputModule()).nonlinear_

    cdef const transfers* tr(self) except NULL:
        return &deref(deref(self._cosmo()).GetInputModule()).transfers_

    cdef const spectra* sp(self) except NULL:
        return &deref(deref(self._cosmo()).GetInputModule()).spectra_

    cdef const lensing* le(self) except NULL:
        return &deref(deref(self._cosmo()).GetInputModule()).lensing_

    cdef const output* op(self) except NULL:
        return &deref(deref(self._cosmo()).GetInputModule()).output_

    # The computed modules, each built on first use: a caller asking for one
    # thermodynamics number must not pay for the perturbations.
    cdef const BackgroundModule* _background(self) except NULL:
        return deref(self._cosmo()).GetBackgroundModule().get()

    cdef const ThermodynamicsModule* _thermo(self) except NULL:
        return deref(self._cosmo()).GetThermodynamicsModule().get()

    cdef const PrimordialModule* _primordial(self) except NULL:
        return deref(self._cosmo()).GetPrimordialModule().get()

    cdef const NonlinearModule* _nonlinear(self) except NULL:
        return deref(self._cosmo()).GetNonlinearModule().get()

    cdef _update_fc_from_pars(self):
        cdef:
            FileContent new_file_content
            char* dumc
            char* dumv

        # _fc will be cleaned up by its destructor
        self._fc = new_file_content

        for name, value in self._pars.items():
            dumcp = name.encode()
            dumc = dumcp
            dumvp = str(value).encode()
            dumv = dumvp
            self._fc.set(dumc, dumv)
        return self

    # The functions struct_cleanup(), empty(), set() and compute() are not neccessary, but they are here to support
    # legacy code and MontePython.
    cpdef struct_cleanup(self):
        return

    cpdef empty(self):
        self._pars = {}
        self.parameters_changed = True
        return self

    def set(self, input_parameters=None, **kwargs):
        # Keywords as well as a dict: Cobaya calls set(**args). Keys that are not
        # Python identifiers ('100*theta_s', 'nu1.m') can only travel in the dict.
        input_parameters = dict(input_parameters or {}, **kwargs)
        if viewdictitems(input_parameters) <= viewdictitems(self._pars):
            return
        self._pars.update(input_parameters)
        self.parameters_changed = True
        return self

    cpdef compute(self, level=None):
        # Legacy shim (MontePython): modules are built lazily on first access,
        # so compute() only pre-triggers the requested level. Its remaining
        # value is that computation errors surface here, where callers expect
        # to catch them.
        if level is None:
            level = ['lensing']
        final_level = level[0].lower()
        if final_level == 'background':
            deref(self._cosmo()).GetBackgroundModule()
        elif final_level == 'thermodynamics':
            deref(self._cosmo()).GetThermodynamicsModule()
        elif final_level == 'perturb':
            deref(self._cosmo()).GetPerturbationsModule()
        elif final_level == 'primordial':
            deref(self._cosmo()).GetPrimordialModule()
        elif final_level == 'nonlinear':
            deref(self._cosmo()).GetNonlinearModule()
        elif final_level == 'transfer':
            deref(self._cosmo()).GetTransferModule()
        elif final_level == 'spectra':
            deref(self._cosmo()).GetSpectraModule()
        elif final_level == 'lensing':
            deref(self._cosmo()).GetLensingModule()
        return self

    cpdef get_input_precision(self):
        return deref(self.pr())

    cpdef get_input_background(self):
        return deref(self.ba())

    cpdef get_input_thermodynamics(self):
        return deref(self.th())

    cpdef get_input_perturbations(self):
        return deref(self.pt())

    cpdef get_input_transfers(self):
        return deref(self.tr())

    cpdef get_input_primordial(self):
        return deref(self.pm())

    cpdef get_input_spectra(self):
        return deref(self.sp())

    cpdef get_input_nonlinear(self):
        return deref(self.nl())

    cpdef get_input_lensing(self):
        return deref(self.le())

    cpdef get_input_output(self):
        return deref(self.op())

    cdef general_cl(self, int lmax, bool is_lensed):
        cdef:
            dict out_dict
            int lmaxpp
            map[string, vector[double]] cl_data

        if is_lensed:
            le = deref(self._cosmo()).GetLensingModule()
            if lmax == -1:
                lmax = deref(le).l_lensed_max_
            cl_data = deref(le).cl_output(lmax)
        else:
            sp = deref(self._cosmo()).GetSpectraModule()
            if lmax == -1:
                lmax = deref(sp).l_max_tot_
            cl_data = deref(sp).cl_output(lmax)
        lmaxpp = lmax + 1

        out_dict = {}
        for element in cl_data:
            key = <bytes> element.first
            key = str(key.decode())
            out_dict[key] = np.asarray(<double[:lmaxpp]> &element.second[0]).copy()
        out_dict['ell'] = np.arange(lmax + 1)
        return out_dict

    cpdef raw_cl_no_copy(self, int lmax=-1):
        cdef:
            dict out_dict
            double[::1] double_view
            int ct_size
            map[string, int] index_map
            vector[double*] output_pointers

        sp = deref(self._cosmo()).GetSpectraModule()
        if lmax == -1:
            lmax = deref(sp).l_max_tot_

        ct_size = deref(sp).ct_size_
        index_map = deref(sp).cl_output_index_map()
        output_pointers.resize(ct_size)
        out_dict = {}
        for element in index_map:
            key = <bytes> element.first
            key = str(key.decode())
            out_dict[key] = np.empty(lmax + 1, dtype=np.double)
            double_view = out_dict[key]
            output_pointers[element.second] = &double_view[0]

        deref(sp).cl_output_no_copy(lmax, output_pointers)
        out_dict['ell'] = np.arange(lmax + 1)
        return out_dict

    cpdef raw_cl(self, int lmax=-1):
        return self.general_cl(lmax, False)

    cpdef lensed_cl(self, int lmax=-1):
        return self.general_cl(lmax, True)

    cpdef lensed_cl_computed(self):
        cdef:
            dict out_dict
            map[string, vector[double]] cl_data

        le = deref(self._cosmo()).GetLensingModule()
        cl_data = deref(le).cl_output_computed()
        out_dict = {}
        for element in cl_data:
            key = <bytes> element.first
            key = str(key.decode())
            out_dict[key] = np.asarray(<double[:element.second.size()]> &element.second[0]).copy()
        out_dict['ell'] = out_dict['ell'].astype(int)
        return out_dict

    cpdef z_of_r (self, double[::1] z_array):
        cdef:
            double tau
            int last_index
            int bg_size
            short long_info
            short inter_normal
            int index_bg_conf_distance
            int index_bg_H
            Py_ssize_t i
            Py_ssize_t z_array_size
            double z
            vector[double] pvecback
            double[::1] r
            double[::1] dzdr

        z_array_size = z_array.shape[0]
        r_arr = np.empty(z_array_size, np.double)
        r = r_arr
        dzdr_arr = np.empty(z_array_size, np.double)
        dzdr = dzdr_arr

        long_info = self.ba().long_info
        inter_normal = self.ba().inter_normal

        background_module = deref(self._cosmo()).GetBackgroundModule()
        bg_size = deref(background_module).bg_size_
        index_bg_conf_distance = deref(background_module).index_bg_conf_distance_
        index_bg_H = deref(background_module).index_bg_H_

        pvecback.resize(bg_size)

        for i in range(z_array_size):
            z = z_array[i]
            deref(background_module).background_tau_of_z(z, &tau)
            deref(background_module).background_at_tau(tau, long_info, inter_normal, &last_index, &pvecback[0])

            # store r
            r[i] = pvecback[index_bg_conf_distance]
            # store dz/dr = H
            dzdr[i] = pvecback[index_bg_H]

        return r_arr, dzdr_arr

    cpdef luminosity_distance(self, double z):
        """
        luminosity_distance(z)
        """
        cdef int index_bg_lum_distance
        background_module = deref(self._cosmo()).GetBackgroundModule()
        index_bg_lum_distance = deref(background_module).index_bg_lum_distance_
        return self.get_background_value_at_z(z, index_bg_lum_distance)

    cdef pk_outputs _pk_output_default(self) except *:
        """Which P(k) a bare pk()/get_pk() call means: non-linear if one was asked for."""
        if self.nl().method == nl_none:
            return pk_linear
        return pk_nonlinear

    cdef pk_general(self, double k, double z, int index_pk, pk_outputs linear_or_nonlinear):
        cdef:
            double pk

        nonlinear_module = deref(self._cosmo()).GetNonlinearModule()
        deref(nonlinear_module).nonlinear_pk_at_k_and_z(linear_or_nonlinear, k, z, index_pk, &pk, NULL)
        return pk

    cpdef pk(self, double k, double z):
        """
        Gives the total matter pk (in Mpc**3) for a given k (in 1/Mpc) and z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        return self.pk_general(k, z, self._nonlinear().index_pk_m_, self._pk_output_default())

    # Gives the cdm+b pk for a given (k,z)
    cpdef pk_cb(self, double k, double z):
        """
        Gives the cdm+b pk (in Mpc**3) for a given k (in 1/Mpc) and z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        return self.pk_general(k, z, self._nonlinear().index_pk_cluster_, self._pk_output_default())

    # Gives the total matter pk for a given (k,z)
    cpdef pk_lin(self, double k, double z):
        """
        Gives the linear total matter pk (in Mpc**3) for a given k (in 1/Mpc) and z

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        return self.pk_general(k, z, self._nonlinear().index_pk_m_, pk_linear)

    # Gives the cdm+b pk for a given (k,z)
    cpdef pk_cb_lin(self,double k,double z):
        """
        Gives the linear cdm+b pk (in Mpc**3) for a given k (in 1/Mpc) and z

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        return self.pk_general(k, z, self._nonlinear().index_pk_cluster_, pk_linear)

    cdef get_pk_general(self, double[:,:,::1] k, double[::1] z, Py_ssize_t k_size, Py_ssize_t z_size, Py_ssize_t mu_size, Py_ssize_t index_pk, pk_outputs linear_or_nonlinear):
        cdef:
            Py_ssize_t index_k
            Py_ssize_t index_z
            Py_ssize_t index_mu
            double pk_val
            double[:,:,::1] pk

        pk_arr = np.empty((k_size, z_size, mu_size), np.double)
        pk = pk_arr
        nonlinear_module = deref(self._cosmo()).GetNonlinearModule()
        # TODO: Consider prange for outer loop
        for index_k in range(k_size):
            for index_z in range(z_size):
                for index_mu in range(mu_size):
                    pk[index_k, index_z, index_mu] = self.pk_general(k[index_k, index_z, index_mu], z[index_z], index_pk, linear_or_nonlinear)

        return pk_arr

    cpdef get_pk(self, double[:,:,::1] k, double[::1] z, Py_ssize_t k_size, Py_ssize_t z_size, Py_ssize_t mu_size):
        """ Fast function to get the power spectrum on a k and z array """
        return self.get_pk_general(k, z, k_size, z_size, mu_size, self._nonlinear().index_pk_m_, self._pk_output_default())

    cpdef get_pk_cb(self, double[:,:,::1] k, double[::1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the power spectrum on a k and z array """
        return self.get_pk_general(k, z, k_size, z_size, mu_size, self._nonlinear().index_pk_cluster_, self._pk_output_default())

    cpdef get_pk_lin(self, double[:,:,::1] k, double[::1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the linear power spectrum on a k and z array """
        return self.get_pk_general(k, z, k_size, z_size, mu_size, self._nonlinear().index_pk_m_, pk_linear)

    cpdef get_pk_cb_lin(self, double[:,:,::1] k, double[::1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the linear power spectrum on a k and z array """
        return self.get_pk_general(k, z, k_size, z_size, mu_size, self._nonlinear().index_pk_cluster_, pk_linear)

    def _z_grid_for_output(self):
        """The redshifts P(k,z) and T_i(k,z) are tabulated at, decreasing to 0.

        One owner for this axis. Both tables live on the perturbation module's
        ln_tau grid -- the Fourier module holds a copy -- so anything that
        combines them, as get_Weyl_pk_and_k_and_z does, is aligned by
        construction rather than by coincidence.
        """
        cdef Py_ssize_t index_tau
        cdef Py_ssize_t ln_tau_size

        perturbations_module = deref(self._cosmo()).GetPerturbationsModule()
        ln_tau_size = deref(perturbations_module).ln_tau_size_

        z = np.empty(ln_tau_size, np.double)
        for index_tau in range(ln_tau_size - 1):
            z[index_tau] = self.z_of_tau(exp(deref(perturbations_module).ln_tau_[index_tau]))
        # The grid ends at today; say so exactly rather than to interpolation accuracy.
        z[ln_tau_size - 1] = 0.

        return z

    def get_pk_and_k_and_z(self, nonlinear=True, only_clustering_species=False, h_units=False):
        """
        Return the stored grid of matter power spectrum values with the k and z it
        was computed on. Useful for building interpolators.

        Parameters
        ----------
        nonlinear : bool
                Return the non-linear spectrum (default) or the linear one
        only_clustering_species : bool
                Exclude massive neutrinos, as galaxy clustering does, rather than
                including everything (default)
        h_units : bool
                Return k in h/Mpc rather than 1/Mpc (default)

        Returns
        -------
        pk : grid of power spectrum values, pk[index_k, index_z], in Mpc**3
        k : vector of k values, k[index_k]
        z : vector of z values, z[index_z]
        """
        cdef vector[double] k_vector
        cdef vector[double] pk_vector
        cdef pk_outputs pk_output
        cdef int index_pk

        pk_output = pk_nonlinear if nonlinear else pk_linear
        nonlinear_module = deref(self._cosmo()).GetNonlinearModule()
        if only_clustering_species:
            index_pk = deref(nonlinear_module).index_pk_cluster_
        else:
            index_pk = deref(nonlinear_module).index_pk_m_
        deref(nonlinear_module).GetPkGrid(pk_output, index_pk, k_vector, pk_vector)

        z = self._z_grid_for_output()
        k = np.asarray(<double[:k_vector.size()]> k_vector.data()).copy()
        # GetPkGrid hands back the table in its stored order, z-major.
        pk = np.asarray(<double[:pk_vector.size()]> pk_vector.data()).reshape(z.size, k.size).T.copy()

        if h_units:
            k = k/self.ba().h

        return pk, k, z

    def get_transfer_and_k_and_z(self, output_format='class', h_units=False):
        """
        Return a dictionary of density and/or velocity transfer function grids with
        the k and z they were computed on. Useful for building interpolators.

        Include at least one of 'dTk' or 'vTk' in 'output'. With the default
        output_format='class' the transfer functions are normalised to curvature
        R=1 at initial time, and the newtonian metric fluctuations phi and psi are
        included; output_format='camb' switches to the CAMB convention, which
        outputs density transfer functions only.

        Parameters
        ----------
        output_format : 'class' or 'camb'
        h_units : bool
                Return k in h/Mpc rather than 1/Mpc (default)

        Returns
        -------
        tk : dict of grids, tk[name][index_k, index_z]
        k : vector of k values, k[index_k]
        z : vector of z values, z[index_z]
        """
        cdef Py_ssize_t index_k
        cdef Py_ssize_t k_size
        cdef int index_md

        if (not self.pt().has_density_transfers) and (not self.pt().has_velocity_transfers):
            raise CosmoSevereError(
                "No transfer functions were computed. Add 'dTk' and/or 'vTk' to the list of outputs.")

        perturbations_module = deref(self._cosmo()).GetPerturbationsModule()
        index_md = deref(perturbations_module).index_md_scalars_

        if deref(perturbations_module).ic_size_[index_md] > 1:
            raise CosmoSevereError(
                "get_transfer_and_k_and_z() assumes adiabatic initial conditions: the grids it "
                "returns have no axis to hold several of them.")

        z = self._z_grid_for_output()

        k_size = deref(perturbations_module).k_size_[index_md]
        k = np.empty(k_size, np.double)
        for index_k in range(k_size):
            k[index_k] = deref(perturbations_module).k_[index_md][index_k]
        if h_units:
            k = k/self.ba().h

        tk = {}
        for index_z in range(z.size):
            # By node, not by redshift: the stored table is what these grids are,
            # and its topmost node does not survive the round trip through z.
            transfers = self.transfer_general(0., index_z, output_format)
            # get_transfer leads with the k column; k is returned separately here.
            transfers.pop('k (h/Mpc)', None)
            for name, values in transfers.items():
                if name not in tk:
                    tk[name] = np.empty((k_size, z.size), np.double)
                tk[name][:, index_z] = values

        return tk, k, z

    def get_Weyl_pk_and_k_and_z(self, nonlinear=False, h_units=False):
        """
        Return the grid of k**4 times the power spectrum of the Weyl potential
        (phi+psi)/2, with the k and z it was computed on.

        This is the matter power spectrum rescaled by [(phi+psi)/2/d_m]**2. The
        factor k**4 is a convention: the Poisson equation carries a k**2, so the
        rescaled spectrum has a shape similar to the matter one.

        Parameters
        ----------
        nonlinear : bool
                Rescale the non-linear matter spectrum; the default is to rescale
                the linear one
        h_units : bool
                Return k in h/Mpc rather than 1/Mpc (default)

        Returns
        -------
        Weyl_pk : grid of values, Weyl_pk[index_k, index_z]
        k : vector of k values, k[index_k]
        z : vector of z values, z[index_z]
        """
        pk, k, z = self.get_pk_and_k_and_z(nonlinear=nonlinear,
                                           only_clustering_species=False,
                                           h_units=h_units)
        tk, _, _ = self.get_transfer_and_k_and_z(output_format='class', h_units=h_units)

        for name in ('phi', 'psi', 'd_m'):
            if name not in tk:
                raise CosmoSevereError(
                    "The Weyl spectrum is built from the '%s' transfer function, which was not "
                    "computed. It needs 'dTk' in 'output' (and, for d_m, 'mPk')." % name)

        # k enters in whichever units the caller asked for, so that the convention
        # is the same on both sides of the rescaling.
        return pk*((tk['phi'] + tk['psi'])/2./tk['d_m'])**2*k[:, None]**4, k, z

    cdef double sigma_general(self, double R, double z, int index_pk, bint h_units) except? -1.:
        cdef double sigma

        # Not NonlinearModule::nonlinear_sigmas_at_z's to refuse: the module also
        # calls it for sigma8 in runs without mPk (halofit for lensing) and
        # with any P_k_max.
        if self.pt().has_pk_matter == constvals.sFALSE:
            raise CosmoSevereError("Power spectrum not computed. In order to get sigma(R, z) you must add mPk to the list of outputs.")
        if self.pt().k_max_for_pk < self.ba().h:
            raise CosmoSevereError("In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc.")

        if h_units:
            R = R/self.ba().h

        nonlinear_module = deref(self._cosmo()).GetNonlinearModule()
        deref(nonlinear_module).nonlinear_sigmas_at_z(R, z, index_pk, out_sigma, &sigma)

        return sigma

    # Gives sigma(R,z) for a given (R,z)
    cpdef sigma(self, double R, double z, bint h_units=False):
        """
        Gives sigma (total matter) for a given R and z

        R is the radius in units of Mpc, so R=8/h gives the usual sigma8(z) --
        unless h_units is set, in which case R is in Mpc/h and R=8 gives it.

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        return self.sigma_general(R, z, self._nonlinear().index_pk_m_, h_units)

    # Gives sigma(R, z) for a given (R, z)
    cpdef sigma_cb(self, double R, double z, bint h_units=False):
        """
        Gives sigma (cdm+b) for a given R and z

        R is the radius in units of Mpc, so R=8/h gives the usual sigma8(z) --
        unless h_units is set, in which case R is in Mpc/h and R=8 gives it.

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        return self.sigma_general(R, z, self._nonlinear().index_pk_cluster_, h_units)

    def effective_f_sigma8(self, z, z_step=0.1):
        """
        Estimate f(z)*sigma8(z) as (d sigma8 / d ln a)(z), by finite difference.

        Parameters
        ----------
        z : float
                Desired redshift
        z_step : float
                Step of the two-sided derivative. For z below z_step the step
                shrinks to z, so the derivative stays two-sided; below z_step/10
                it goes one-sided instead.

        Returns
        -------
        (d sigma8/d ln a)(z) (dimensionless)
        """
        # d sigma8/d ln a = -(d sigma8/dz)*(1+z)
        if z < z_step/10.:
            z_step = z_step/10.
            return (self.sigma(8, z, h_units=True)
                    - self.sigma(8, z + z_step, h_units=True))/z_step*(1 + z)

        if z < z_step:
            z_step = z

        return (self.sigma(8, z - z_step, h_units=True)
                - self.sigma(8, z + z_step, h_units=True))/(2.*z_step)*(1 + z)

    # Gives effective logarithmic slope of P_L(k,z) (total matter) for a given (k,z)
    cpdef pk_tilt(self, double k, double z):
        """
        Gives effective logarithmic slope of P_L(k,z) (total matter) for a given k and z
        (k is the wavenumber in units of 1/Mpc, z is the redshift, the output is dimensionless)

        .. note::

            there is an additional check to verify whether output contains `mPk` and whether k is in the right range

        """
        cdef double pk_tilt

        nonlinear_module = deref(self._cosmo()).GetNonlinearModule()
        deref(nonlinear_module).nonlinear_pk_tilt_at_k_and_z(pk_linear, k, z, deref(nonlinear_module).index_pk_m_, &pk_tilt)
        return pk_tilt


    cpdef age(self):
        return self._background().age_

    cpdef h(self):
        return self.ba().h

    cpdef n_s(self):
        return self._primordial().n_s_

    cpdef A_s(self):
        return self._primordial().A_s_

    cpdef tau_reio(self):
        return self._thermo().tau_reionization_

    cpdef Omega_m(self):
        return self._background().Omega0_m_

    cpdef Omega_r(self):
        return self._background().Omega0_r_

    cpdef theta_s_100(self):
        return 100.*self._thermo().rs_rec_/self._thermo().ra_rec_

    cpdef theta_star_100(self):
        return 100.*self._thermo().rs_star_/self._thermo().ra_star_

    cpdef theta_d_100(self):
        if (self.th().compute_damping_scale == constvals.sFALSE):
            raise CosmoSevereError(r"Photon damping scale has not been computed - you must add 'compute damping scale':'yes'.")
        return 100.*self._thermo().rd_rec_/self._thermo().ra_rec_

    cpdef Omega_Lambda(self):
        return self._background().GetOmega0Species(b"Lambda")

    cpdef Omega_g(self):
        return self.ba().Omega0_g

    cpdef Omega_b(self):
        return self.ba().Omega0_b

    cpdef omega_b(self):
        return self.ba().Omega0_b*self.ba().h*self.ba().h

    cpdef Neff(self):
        return self._background().Neff_

    cpdef f_idr_drmd(self):
        return self._background().f_idr_drmd()

    cpdef z_dec_drmd(self):
        return self._background().z_dec_drmd()

    cpdef k_eq(self):
        return self._background().a_eq_*self._background().H_eq_

    cpdef sigma8(self):
        return self._nonlinear().GetSigma8(self._nonlinear().index_pk_m_)

    cpdef sigma8_cb(self):
        return self._nonlinear().GetSigma8(self._nonlinear().index_pk_cluster_)

    cpdef rs_drag(self):
        return self._thermo().rs_d_

    cpdef z_reio(self):
        return self._thermo().z_reionization_

    cdef get_background_value_at_z(self, double z, int index_bg):
        cdef:
            double output_value
            double tau
            int bg_size
            int last_index
            short inter_normal
            short long_info
            vector[double] pvecback

        long_info = self.ba().long_info
        inter_normal = self.ba().inter_normal

        background_module = deref(self._cosmo()).GetBackgroundModule()
        bg_size = deref(background_module).bg_size_

        pvecback.resize(bg_size)
        deref(background_module).background_tau_of_z(z, &tau)
        deref(background_module).background_at_tau(tau, long_info, inter_normal, &last_index, &pvecback[0])
        output_value = pvecback[index_bg]
        return output_value

    cdef get_thermodynamics_value_at_z(self, double z, int index_th):
        cdef:
            double tau
            int last_index
            vector[double] pvecback
            vector[double] pvecthermo
            short long_info
            short inter_normal
            int bg_size
            int th_size
            short inter_mode
            double output_value

        long_info = self.ba().long_info
        inter_normal = self.ba().inter_normal

        background_module = deref(self._cosmo()).GetBackgroundModule()
        thermodynamics_module = deref(self._cosmo()).GetThermodynamicsModule()
        bg_size = deref(background_module).bg_size_
        th_size = deref(thermodynamics_module).th_size_

        pvecback.resize(bg_size)
        pvecthermo.resize(th_size)

        inter_mode = deref(thermodynamics_module).inter_normal_

        deref(background_module).background_tau_of_z(z, &tau)
        deref(background_module).background_at_tau(tau, long_info, inter_normal, &last_index, &pvecback[0])
        deref(thermodynamics_module).thermodynamics_at_z(z, inter_mode, &last_index, &pvecback[0], &pvecthermo[0])

        output_value = pvecthermo[index_th]
        return output_value

    cpdef angular_distance(self, double z):
        """
        angular_distance(z)

        Return the angular diameter distance (exactly, the quantity defined by Class
        as index_bg_ang_distance in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int index_bg_ang_distance
        background_module = deref(self._cosmo()).GetBackgroundModule()
        index_bg_ang_distance = deref(background_module).index_bg_ang_distance_
        return self.get_background_value_at_z(z, index_bg_ang_distance)

    def angular_distance_from_to(self, z1, z2):
        """
        Angular diameter distance of an object at z2 as seen by an observer at z1,
        that is, sin_K((chi2 - chi1)*sqrt(|K|))/sqrt(|K|)/(1 + z2). Zero if z1 >= z2.

        Parameters
        ----------
        z1 : float
                Observer redshift
        z2 : float
                Source redshift

        Returns
        -------
        d_A(z1, z2) in Mpc
        """
        cdef int index_bg_conf_distance
        cdef double K

        if z1 >= z2:
            return 0.

        background_module = deref(self._cosmo()).GetBackgroundModule()
        index_bg_conf_distance = deref(background_module).index_bg_conf_distance_
        chi1 = self.get_background_value_at_z(z1, index_bg_conf_distance)
        chi2 = self.get_background_value_at_z(z2, index_bg_conf_distance)

        K = self.ba().K
        if K == 0.:
            return (chi2 - chi1)/(1 + z2)
        if K > 0.:
            return np.sin(np.sqrt(K)*(chi2 - chi1))/np.sqrt(K)/(1 + z2)
        return np.sinh(np.sqrt(-K)*(chi2 - chi1))/np.sqrt(-K)/(1 + z2)

    cpdef scale_independent_growth_factor(self, z):
        """
        scale_independent_growth_factor(z)

        Return the scale invariant growth factor D(a) for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_D in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int index_bg_D
        background_module = deref(self._cosmo()).GetBackgroundModule()
        index_bg_D = deref(background_module).index_bg_D_
        return self.get_background_value_at_z(z, index_bg_D)

    cpdef scale_independent_growth_factor_f(self, z):
        """
        scale_independent_growth_factor_f(z)

        Return the scale invariant growth factor f(z)=d ln D / d ln a for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_f in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int index_bg_f
        background_module = deref(self._cosmo()).GetBackgroundModule()
        index_bg_f = deref(background_module).index_bg_f_
        return self.get_background_value_at_z(z, index_bg_f)

    @cython.cdivision(True)
    cpdef z_of_tau(self, double tau):
        """
        Redshift corresponding to a given conformal time.

        Parameters
        ----------
        tau : float
                Conformal time
        """
        cdef:
            double output_value
            double z
            int bg_size_short
            int index_bg_a
            int last_index
            short inter_normal
            short short_info
            vector[double] pvecback

        short_info = self.ba().short_info
        inter_normal = self.ba().inter_normal

        background_module = deref(self._cosmo()).GetBackgroundModule()
        bg_size_short = deref(background_module).bg_size_short_
        index_bg_a = deref(background_module).index_bg_a_

        pvecback.resize(bg_size_short)
        deref(background_module).background_at_tau(tau, short_info, inter_normal, &last_index, &pvecback[0])
        z = 1./pvecback[index_bg_a] - 1.
        return z

    cpdef Hubble(self, double z):
        """
        Hubble(z)

        Return the Hubble rate (exactly, the quantity defined by Class as index_bg_H
        in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int index_bg_H
        background_module = deref(self._cosmo()).GetBackgroundModule()
        index_bg_H = deref(background_module).index_bg_H_
        return self.get_background_value_at_z(z, index_bg_H)

    cpdef Om_m(self, double z):
        """
        Omega_m(z)

        Return the matter density fraction (exactly, the quantity defined by Class as index_bg_Omega_m
        in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int index_bg_Omega_m
        background_module = deref(self._cosmo()).GetBackgroundModule()
        index_bg_Omega_m = deref(background_module).index_bg_Omega_m_
        return self.get_background_value_at_z(z, index_bg_Omega_m)


    cpdef ionization_fraction(self, double z):
        """
        ionization_fraction(z)

        Return the ionization fraction for a given redshift z

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int index_th_xe
        thermodynamics_module = deref(self._cosmo()).GetThermodynamicsModule()
        index_th_xe = deref(thermodynamics_module).index_th_xe_
        return self.get_thermodynamics_value_at_z(z, index_th_xe)

    cpdef baryon_temperature(self, double z):
        """
        baryon_temperature(z)

        Give the baryon temperature for a given redshift z

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int index_th_Tb
        thermodynamics_module = deref(self._cosmo()).GetThermodynamicsModule()
        index_th_Tb = deref(thermodynamics_module).index_th_Tb_
        return self.get_thermodynamics_value_at_z(z, index_th_Tb)

    cpdef T_cmb(self):
        """
        Return the CMB temperature
        """
        return self.ba().T_cmb

    # redundent with a previous Omega_m() funciton,
    # but we leave it not to break compatibility
    cpdef Omega0_m(self):
        """
        Return the sum of Omega0 for all non-relativistic components
        """
        return self.Omega_m()

    cpdef get_background(self):
        """
        Return an array of the background quantities at all times.

        Parameters
        ----------

        Returns
        -------
        background : dictionary containing background.
        """
        cdef:
            dict background
            double[:,::1] view
            int i
            int number_of_titles
            int timesteps
            string titles
            vector[double] data

        background_module = deref(self._cosmo()).GetBackgroundModule()
        deref(background_module).background_output_titles(titles)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.strip("\t").split("\t")
        number_of_titles = len(names)
        timesteps = deref(background_module).bt_size_

        data.resize(timesteps*number_of_titles)
        deref(background_module).background_output_data(number_of_titles, &data[0])

        background = {}
        view = <double[:timesteps,:number_of_titles]> &data[0]
        for i in range(number_of_titles):
            background[names[i]] = np.asarray(view[:, i]).copy()

        return background

    cpdef get_thermodynamics(self):
        """
        Return the thermodynamics quantities.

        Returns
        -------
        thermodynamics : dictionary containing thermodynamics.
        """
        cdef:
            dict thermodynamics
            double[:,::1] view
            int i
            int number_of_titles
            int timesteps
            string titles
            vector[double] data

        thermodynamics_module = deref(self._cosmo()).GetThermodynamicsModule()

        deref(thermodynamics_module).thermodynamics_output_titles(titles)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.strip("\t").split("\t")
        number_of_titles = len(names)
        timesteps = deref(thermodynamics_module).tt_size_

        data.resize(timesteps*number_of_titles)

        deref(thermodynamics_module).thermodynamics_output_data(number_of_titles, &data[0])

        thermodynamics = {}
        view = <double[:timesteps,:number_of_titles]> &data[0]
        for i in range(number_of_titles):
            thermodynamics[names[i]] = np.asarray(view[:, i]).copy()

        return thermodynamics

    cpdef get_primordial(self):
        """
        Return the primordial scalar and/or tensor spectrum depending on 'modes'.
        'output' must be set to something, e.g. 'tCl'.

        Returns
        -------
        primordial : dictionary containing k-vector and primordial scalar and tensor P(k).
        """
        cdef:
            string titles
            vector[double] data
            int timesteps
            int number_of_titles
            double[:,::1] view
            dict primordial
            int i

        primordial_module = deref(self._cosmo()).GetPrimordialModule()

        deref(primordial_module).primordial_output_titles(titles)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.strip("\t").split("\t")
        number_of_titles = len(names)
        timesteps = deref(primordial_module).lnk_size_

        data.resize(timesteps*number_of_titles)
        deref(primordial_module).primordial_output_data(number_of_titles, &data[0])

        primordial = {}
        view = <double[:timesteps,:number_of_titles]> &data[0]
        for i in range(number_of_titles):
            primordial[names[i]] = np.asarray(view[:, i]).copy()

        return primordial


    @cython.returns(dict)
    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.cdivision(True)
    @cython.ccall
    cpdef get_perturbations(self):
        """
        Return scalar, vector and/or tensor perturbations as arrays for requested
        k-values.

        .. note::

            you need to specify both 'k_output_values', and have some
            perturbations computed, for instance by setting 'output' to 'tCl'.

        Returns
        -------
        perturbations : dict of array of dicts
                perturbations['scalar'] is an array of length 'k_output_values' of
                dictionary containing scalar perturbations.
                Similar for perturbations['vector'] and perturbations['tensor'].
        """

        perturbations = {}

        if self.pt().k_output_values_num<1:
            return perturbations

        cdef:
            Py_ssize_t i
            Py_ssize_t j
            Py_ssize_t k_size
            Py_ssize_t number_of_titles
            Py_ssize_t timesteps
            dict tmpdict
            double[:,::1] view
            double* thedata_ptr
            int* thesizes
            list names
            list tmparray

        perturbation_module = deref(self._cosmo()).GetPerturbationsModule()
        k_size = self.pt().k_output_values_num
        for mode in ['scalar','vector','tensor']:
            if mode=='scalar' and self.pt().has_scalars:
                thetitles = <bytes> deref(perturbation_module).scalar_titles_
                thesizes = <int*> deref(perturbation_module).size_scalar_perturbation_data_
            elif mode=='vector' and self.pt().has_vectors:
                thetitles = <bytes> deref(perturbation_module).vector_titles_
                thesizes = <int*> deref(perturbation_module).size_vector_perturbation_data_
            elif mode=='tensor' and self.pt().has_tensors:
                thetitles = <bytes> deref(perturbation_module).tensor_titles_
                thesizes = <int*> deref(perturbation_module).size_tensor_perturbation_data_
            else:
                continue
            thetitles = str(thetitles.decode())
            names = thetitles.strip("\t").split("\t")
            number_of_titles = len(names)
            tmparray = []
            if number_of_titles == 0:
                continue
            for j in range(k_size):
                timesteps = thesizes[j]//number_of_titles
                tmpdict = {}
                if mode == 'scalar':
                    thedata_ptr = <double*>deref(perturbation_module).scalar_perturbations_data_[j].data()
                elif mode == 'vector':
                    thedata_ptr = <double*>deref(perturbation_module).vector_perturbations_data_[j].data()
                else:
                    thedata_ptr = <double*>deref(perturbation_module).tensor_perturbations_data_[j].data()
                view = <double[:timesteps,:number_of_titles]> thedata_ptr
                for i in range(number_of_titles):
                    tmpdict[names[i]] = np.asarray(view[:,i])
                tmparray.append(tmpdict)
            perturbations[mode] = tmparray

        return perturbations

    cpdef get_transfer(self, double z=0., output_format='class'):
        """
        Return the density and/or velocity transfer functions for all initial
        conditions today. You must include 'mTk' and/or 'vCTk' in the list of
        'output'. The transfer functions can also be computed at higher redshift z
        provided that 'z_pk' has been set and that 0<z<z_pk.

        Parameters
        ----------
        z  : redshift (default = 0)
        output_format  : ('class' or 'camb') Format transfer functions according to
                         CLASS convention (default) or CAMB convention.

        Returns
        -------
        tk : dictionary containing transfer functions.
        """
        return self.transfer_general(z, -1, output_format)

    cdef transfer_general(self, double z, int index_tau, output_format):
        """Transfer functions at a redshift, or -- when index_tau is not negative
        -- at that node of the stored late-time sampling, which is exact where
        converting the node to a redshift and back is not."""
        cdef:
            string titles
            vector[double] data
            string ic_info
            string ic_suffix
            file_format outf
            Py_ssize_t number_of_titles
            Py_ssize_t timesteps
            Py_ssize_t size_ic_data
            Py_ssize_t ic_num
            dict tmpdict
            dict transfers
            double[:,:,::1] view
            int index_ic
            int i

        if (not self.pt().has_density_transfers) and (not self.pt().has_velocity_transfers):
            return {}

        if output_format == 'camb':
            outf = file_format.camb_format
        else:
            outf = file_format.class_format

        perturbations_module = deref(self._cosmo()).GetPerturbationsModule()
        index_md = deref(perturbations_module).index_md_scalars_

        deref(perturbations_module).perturb_output_titles(outf, titles)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.strip("\t").split("\t")
        number_of_titles = len(names)
        timesteps = deref(perturbations_module).k_size_[index_md]

        size_ic_data = timesteps*number_of_titles
        ic_num = deref(perturbations_module).ic_size_[index_md]

        data.resize(size_ic_data*ic_num)

        if index_tau < 0:
            deref(perturbations_module).perturb_output_data(outf, z, number_of_titles, &data[0])
        else:
            deref(perturbations_module).perturb_output_data_at_index_tau(outf, index_tau, number_of_titles, &data[0])

        transfers = {}
        view = <double[:ic_num,:timesteps,:number_of_titles]> &data[0]

        for index_ic in range(ic_num):
            deref(perturbations_module).perturb_output_firstline_and_ic_suffix(index_ic, ic_info, ic_suffix)
            ic_key = <bytes> ic_suffix
            ic_key = str(ic_key.decode())

            tmpdict = {}
            for i in range(number_of_titles):
                tmpdict[names[i]] = np.asarray(view[index_ic, :, i]).copy()

            if ic_num == 1:
                transfers = tmpdict
            else:
                transfers[ic_key] = tmpdict

        return transfers

    @cython.cdivision(True)
    cdef get_slowroll_parameters(self):
        cdef:
            (double, double, double) epsilons
            double C
            double alpha_s
            double eps1
            double eps2
            double eps23
            double n_s
            double r

        primordial_module = deref(self._cosmo()).GetPrimordialModule()
        r = deref(primordial_module).r_
        n_s = deref(primordial_module).n_s_
        alpha_s = deref(primordial_module).alpha_s_
        C = EULER + LOGE2 - 2.0

        eps1 = r*(1./16. + C/16.*(r/8. + n_s - 1.))
        eps2 = -n_s + 1. + C*alpha_s - r*(1./8. + 1./8.*(n_s - 1.)*(C - 1.5)) - (r/8.)**2*(C - 1.)
        eps23 = 1./8.*(r**2/8. + (n_s - 1.)*r - 8.*alpha_s)
        epsilons = eps1, eps2, eps23
        return epsilons

    @cython.cdivision(True)
    cpdef get_current_derived_parameters(self, names):
        """
        get_current_derived_parameters(names)

        Return a dictionary containing an entry for all the names defined in the
        input list.

        Parameters
        ----------
        names : list
                Derived parameters that can be asked from Monte Python, or
                elsewhere. Besides the flat names, "<species key>.<field>"
                reads a species' own derived quantities, e.g. "Fluid.f_peak"
                for the peak EDE fraction f_EDE of a pheno-axion fluid,
                "Fluid.z_peak" for the redshift where it peaks, or
                "ScalarField.f_peak" for the same quantity of an axion field
                integrated exactly. See _SPECIES_DERIVED_FIELDS for the fields.

        Returns
        -------
        derived : dict

        .. warning::

            This method used to take as an argument directly the data class from
            Monte Python. To maintain compatibility with this old feature, a
            check is performed to verify that names is indeed a list. If not, it
            returns a TypeError. The old version of this function, when asked
            with the new argument, will raise an AttributeError.

        """
        cdef:
            dict derived
            double value
            double eps1
            double eps2
            double eps23
            (double, double, double) epsilons

        if not isinstance(names, list):
            raise TypeError("Deprecated")

        derived = {}
        for name in names:
            # rs_rec_h = rs_rec * h, and likewise for the other distances below
            key, scale = name, 1.
            if name in _H_SCALED_DISTANCES:
                key, scale = name[:-2], self.ba().h
            if key == 'h':
                value = self.ba().h
            elif key == 'H0':
                value = self.ba().h*100
            elif key == 'Omega0_lambda' or key == 'Omega_Lambda':
                value = self._background().GetOmega0Species(b"Lambda")
            elif key == 'Omega0_fld':
                value = self._background().GetOmega0Species(b"Fluid")
            elif key == 'age':
                value = self.age()
            elif key == 'conformal_age':
                value = self._background().conformal_age_
            #TODO: Need to wrap NCDM class to get the masses
            #elif key == 'm_ncdm_in_eV':
            #    value = self.ba().m_ncdm_in_eV[0]
            elif key == 'm_ncdm_tot':
                # TODO: Is this really what we want??
                value = self._background().GetOmega0NcdmTot()*self.ba().h*self.ba().h*93.14
            elif key == 'Neff':
                value = self.Neff()
            elif key == 'Omega_m':
                value = self.Omega_m()
            elif key == 'omega_m':
                value = self.Omega_m()/self.ba().h/self.ba().h
            elif key == 'Omega_g':
                value = self.ba().Omega0_g
            elif key == 'xi_idr':
                value = self._background().GetSpeciesParam(b"IDM_DR_IDR", b"T_idr")/self.ba().T_cmb
            elif key == 'N_dg':
                value = self._background().GetOmega0Species(b"IDR")/self.ba().Omega0_g*8./7.*pow(11./4.,4./3.)
            elif key == 'Gamma_0_nadm':
                bam = self._background()
                value = bam.GetSpeciesParam(b"IDM_DR_IDR", b"a_idm_dr")*(4./3.)*(self.ba().h*self.ba().h*bam.GetOmega0Species(b"IDR"))
            elif key == 'a_dark':
                value = self._background().GetSpeciesParam(b"IDM_DR_IDR", b"a_idm_dr")
            elif key == 'tau_reio':
                value = self.tau_reio()
            elif key == 'z_reio':
                value = self.z_reio()
            elif key == 'z_rec':
                value = self._thermo().z_rec_
            elif key == 'tau_rec':
                value = self._thermo().tau_rec_
            elif key == 'rs_rec':
                value = self._thermo().rs_rec_
            elif key == 'ds_rec':
                value = self._thermo().ds_rec_
            elif key == 'ra_rec':
                value = self._thermo().ra_rec_
            elif key == 'da_rec':
                value = self._thermo().da_rec_
            elif key == 'z_star':
                value = self._thermo().z_star_
            elif key == 'tau_star':
                value = self._thermo().tau_star_
            elif key == 'rs_star':
                value = self._thermo().rs_star_
            elif key == 'ds_star':
                value = self._thermo().ds_star_
            elif key == 'ra_star':
                value = self._thermo().ra_star_
            elif key == 'da_star':
                value = self._thermo().da_star_
            elif key == 'rd_star':
                value = self._thermo().rd_star_
            elif key == 'z_d':
                value = self._thermo().z_d_
            elif key == 'tau_d':
                value = self._thermo().tau_d_
            elif key == 'ds_d':
                value = self._thermo().ds_d_
            elif key == 'rs_d':
                value = self._thermo().rs_d_
            elif key == '100*theta_s':
                value = self.theta_s_100()
            elif key == '100*theta_star':
                value = self.theta_star_100()
            elif key == '100*theta_d':
                value = self.theta_d_100()
            elif key == 'YHe':
                value = self._thermo().YHe_
            elif key in ('D_H', 'He3_H', 'Li7_H', 'Li6_H', 'eta10'):
                # Primordial abundances from the solved BBN network. The
                # interpolation table yields Y_He alone, so asking for these
                # without 'YHe = network' is a request the run cannot answer;
                # say so rather than returning a zero that looks like a result.
                thm = self._thermo()
                if not thm.bbn_network_ran_:
                    raise CosmoSevereError(
                        "%s is only computed by the solved BBN network. Set "
                        "'YHe' to 'network' to obtain it; the sBBN interpolation "
                        "table provides Y_He only." % key)
                if key == 'D_H':
                    value = thm.DoverH_
                elif key == 'He3_H':
                    value = thm.He3overH_
                elif key == 'Li7_H':
                    value = thm.Li7overH_
                elif key == 'Li6_H':
                    value = thm.Li6overH_
                else:
                    value = thm.bbn_eta10_
            elif key == 'n_e':
                value = self._thermo().n_e_
            elif key == 'A_s':
                value = self.A_s()
            elif key == 'ln10^{10}A_s':
                value = log(1.e10*self.A_s())
            elif key == 'n_s':
                value = self._primordial().n_s_
            elif key == 'alpha_s':
                value = self._primordial().alpha_s_
            elif key == 'beta_s':
                value = self._primordial().beta_s_
            elif key == 'r':
                # This is at the pivot scale
                value = self._primordial().r_
            elif key == 'r_0002':
                # at k_pivot = 0.002/Mpc
                prm = self._primordial()
                value = prm.r_*(0.002/self.pm().k_pivot)**(prm.n_t_ - prm.n_s_ - 1 + 0.5*prm.alpha_s_*log(0.002/self.pm().k_pivot))
            elif key == 'n_t':
                value = self._primordial().n_t_
            elif key == 'alpha_t':
                value = self._primordial().alpha_t_
            elif key == 'V_0':
                value = self.pm().V0
            elif key == 'V_1':
                value = self.pm().V1
            elif key == 'V_2':
                value = self.pm().V2
            elif key == 'V_3':
                value = self.pm().V3
            elif key == 'V_4':
                value = self.pm().V4
            elif key in ('epsilon_V', 'eta_V', 'ksi_V^2'):
                epsilons = self.get_slowroll_parameters()
                eps1, eps2, eps23 = epsilons
                if key == 'epsilon_V':
                    value = eps1*((1. - eps1/3. + eps2/6.)/(1. - eps1/3.))**2
                elif key == 'eta_V':
                    value = (2.*eps1 - eps2/2. - 2./3.*eps1**2 + 5./6.*eps1*eps2 - eps2**2/12. - eps23/6.)/(1. - eps1/3.)
                else:
                    value = 2.*(1. - eps1/3. + eps2/6.)*(2.*eps1**2 - 3./2.*eps1*eps2 + eps23/4.)/(1. - eps1/3.)**2
            elif key == 'exp_m_2_tau_As':
                value = exp(-2*self.tau_reio())*self.A_s()
            elif key == 'phi_min':
                value = self._primordial().phi_min_
            elif key == 'phi_max':
                value = self._primordial().phi_max_
            elif key == 'sigma8':
                value = self.sigma8()
            elif key == 'sigma8_cb':
                value = self.sigma8_cb()
            elif key == 'k_eq':
                value = self.k_eq()
            elif key == 'z_dec_drmd':
                value = self.z_dec_drmd()
            elif key == 'f_idr_drmd':
                value = self.f_idr_drmd()
            elif '.' in key and key.rpartition('.')[2] in _SPECIES_DERIVED_FIELDS:
                instance, _, field = key.rpartition('.')
                bam = self._background()
                if not bam.HasSpecies(instance.encode('utf-8')):
                    raise CosmoSevereError(
                        "%s is not a derived parameter of this cosmology: it has "
                        "no species '%s'" % (key, instance))
                if not bam.HasSpeciesParam(instance.encode('utf-8'), field.encode('utf-8')):
                    raise CosmoSevereError(
                        "%s is not a derived parameter of this cosmology: species "
                        "'%s' does not report '%s'" % (key, instance, field))
                value = bam.GetSpeciesParam(instance.encode('utf-8'), field.encode('utf-8'))
            else:
                raise CosmoSevereError("%s was not recognized as a derived parameter" % name)
            derived[name] = value*scale
        return derived

    cdef _nonlinear_scales(self, double[::1] z, Py_ssize_t z_size, bint cb):
        cdef:
            Py_ssize_t index_z
            double k_nl_val
            double k_nl_cb_val
            double[:] k_nl

        k_nl_arr = np.empty(z_size, np.double)
        k_nl = k_nl_arr
        nonlinear_module = self._nonlinear()
        for index_z in range(z_size):
            nonlinear_module.nonlinear_k_nl_at_z(z[index_z], &k_nl_val, &k_nl_cb_val)
            k_nl[index_z] = k_nl_cb_val if cb else k_nl_val

        return k_nl_arr

    cpdef nonlinear_scale(self, double[::1] z, Py_ssize_t z_size):
        """
        nonlinear_scale(z, z_size)

        Return the nonlinear scale for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        return self._nonlinear_scales(z, z_size, False)

    cpdef nonlinear_scale_cb(self, double[::1] z, Py_ssize_t z_size):
        """
        nonlinear_scale_cb(z, z_size)

        Return the nonlinear scale of the cdm+baryon spectrum for all the
        redshift specified in z, of size z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        return self._nonlinear_scales(z, z_size, True)

    def __call__(self, ctx):
        """
        Function to interface with CosmoHammer

        Parameters
        ----------
        ctx : context
                Contains several dictionaries storing data and cosmological
                information

        """
        data = ctx.get('data')  # recover data from the context

        # Set the module to the current values
        self._pars = dict(data.cosmo_arguments)
        self.parameters_changed = True
        self.compute()

        # Compute the derived paramter value and store them
        params = ctx.getData()
        self.get_current_derived_parameters(
            data.get_mcmc_parameters(['derived']))
        for elem in data.get_mcmc_parameters(['derived']):
            data.mcmc_parameters[elem]['current'] /= \
                data.mcmc_parameters[elem]['scale']
            params[elem] = data.mcmc_parameters[elem]['current']

        ctx.add('boundary', True)
        # Store itself into the context, to be accessed by the likelihoods
        ctx.add('cosmo', self)

    cdef get_pk_array_general(self, double[::1] k,  double[::1] z, int k_size, int z_size, nonlinear):
        """ Fast function to get the power spectrum on a k and z array """
        cdef:
            cdef double[::1] pk
            cdef double[::1] pk_cb
            int status
            pk_outputs linear_or_nonlinear
        if nonlinear == 0:
            linear_or_nonlinear = pk_linear
        else:
            linear_or_nonlinear = pk_nonlinear

        pk_arr = np.empty(k_size*z_size, np.double)
        pk = pk_arr
        pk_cb_arr = np.empty(k_size*z_size, np.double)
        pk_cb = pk_cb_arr

        nonlinear_module = deref(self._cosmo()).GetNonlinearModule()
        deref(nonlinear_module).nonlinear_pks_at_kvec_and_zvec(linear_or_nonlinear, &k[0], k_size, &z[0], z_size, &pk[0], &pk_cb[0])

        if not deref(nonlinear_module).has_pk_cb_:
            # The module stores P_cb apart only when it differs from P_m.
            pk_cb_arr = pk_arr
        return pk_arr, pk_cb_arr

    cpdef get_pk_array(self, double[::1] k,  double[::1] z, int k_size, int z_size, nonlinear):
        return self.get_pk_array_general(k, z, k_size, z_size, nonlinear)[0]

    cpdef get_pk_cb_array(self, double[::1] k,  double[::1] z, int k_size, int z_size, nonlinear):
        return self.get_pk_array_general(k, z, k_size, z_size, nonlinear)[1]

    cpdef Omega0_k(self):
        """ Curvature contribution """
        return self.ba().Omega0_k

    cpdef Omega0_cdm(self):
        return self._background().GetOmega0Species(b"CDM")
