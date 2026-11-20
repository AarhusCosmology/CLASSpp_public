# distutils: language = c++
"""
.. module:: classy._core
    :synopsis: The compiled core of classy: the Cosmology's lifetime, one method per C++
               call, exception translation and loops over input arrays
.. moduleauthor:: Karim Benabed <benabed@iap.fr>
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
.. moduleauthor:: Julien Lesgourgues <lesgourg@cern.ch>
.. moduleauthor:: Thomas Tram <thomas.tram@phys.au.dk>

This module defines PyCosmology, the compiled half of classy.Class: the C++ Cosmology's
lifetime and parameter gate, one method per C++ call, exception translation and loops over
input arrays. classy.Class (_class.py) subclasses it.

"""
from libc.string cimport memcpy

from libcpp cimport bool
from libcpp.memory cimport unique_ptr
from libcpp.memory cimport shared_ptr
from libcpp.memory cimport make_shared
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

import cython
cimport cython
from cython.operator cimport dereference as deref

import numpy as np
cimport numpy as np
import os

# The exceptions are plain Python (_errors.py), defined once for the whole package, so that
# `except classy.CosmoSevereError` catches what this module raises. Imported by the package's
# runtime name, which is not "classy" in a build renamed with CLASS_PYTHON_MODULE_NAME.
import importlib
_errors = importlib.import_module(__name__.rpartition('.')[0] + '._errors')
CosmoError = _errors.CosmoError
CosmoSevereError = _errors.CosmoSevereError
CosmoComputationError = _errors.CosmoComputationError

# The package's data (bbn/, hyrec/), the default class_dir.
_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

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

cdef object _ndarray(const vector[double]& values):
    """A numpy copy of a C++ vector, which the C++ object may free after this call."""
    array = np.empty(values.size(), np.double)
    cdef double[::1] view = array
    if values.size() > 0:
        memcpy(&view[0], values.data(), values.size()*sizeof(double))
    return array

cdef dict _columns_dict(const NamedColumns& columns):
    """A table of named columns from the C++ modules as a dict of arrays, in column order."""
    cdef dict out = {}
    cdef size_t index
    for index in range(columns.size()):
        out[(<bytes> columns[index].first).decode()] = _ndarray(columns[index].second)
    return out

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
        # By value: Cython wraps a reference returned under except + in a proxy deref cannot take.
        InputModulePtr GetInputModule() except +raise_my_py_error
        BackgroundModulePtr& GetBackgroundModule() except +raise_my_py_error
        ThermodynamicsModulePtr& GetThermodynamicsModule() except +raise_my_py_error
        PerturbationsModulePtr& GetPerturbationsModule() except +raise_my_py_error
        PrimordialModulePtr& GetPrimordialModule() except +raise_my_py_error
        NonlinearModulePtr& GetNonlinearModule() except +raise_my_py_error
        TransferModulePtr& GetTransferModule() except +raise_my_py_error
        SpectraModulePtr& GetSpectraModule() except +raise_my_py_error
        LensingModulePtr& GetLensingModule() except +raise_my_py_error
        double GetDerivedParameter(const string& name) except +raise_my_py_error

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


def _bessel_direct(int K, double beta, int[::1] l, double[::1] x):
    """Phi_l^beta(x) by direct recurrence, as an (n_l, n_x) array (see classy._bessel)."""
    Phi = np.empty((l.shape[0], x.shape[0]), dtype=np.double)
    cdef double[:, ::1] out = Phi
    hyperspherical_bessel_direct_vector(K, beta, &l[0], l.shape[0], &x[0], x.shape[0], &out[0, 0])
    return Phi


def _bessel_interpolate(int K, double beta, int[::1] l_unique, double[::1] x, double xmin,
                        double xmax, double sampling, int l_WKB, double phiminabs,
                        bint derivatives):
    """Phi_l^beta(x) on a grid over [xmin, xmax] with Hermite interpolation, one row per l in
    l_unique (sorted, unique); with derivatives, also dPhi and d2Phi, else None for both."""
    cdef int nl = l_unique.shape[0]
    cdef int nx = x.shape[0]
    Phi = np.empty((nl, nx), dtype=np.double)
    cdef double[:, ::1] phi = Phi
    cdef double[:, ::1] dphi
    cdef double[:, ::1] d2phi
    dPhi = d2Phi = None
    if derivatives:
        dPhi = np.empty((nl, nx), dtype=np.double)
        d2Phi = np.empty((nl, nx), dtype=np.double)
        dphi = dPhi
        d2phi = d2Phi
    # Heap-allocate with `new` so the C++ constructor runs and fills the
    # struct. A stack `cdef HyperInterpStruct his` would bypass the constructor
    # (Cython zero-initializes the raw bytes), leaving the std::vector members
    # in an invalid state and producing sporadic NaNs/heap corruption.
    cdef HyperInterpStruct* his = new HyperInterpStruct(K, beta, nl, &l_unique[0], xmin, xmax,
                                                        sampling, l_WKB, phiminabs)
    cdef int i
    try:
        for i in range(nl):
            if derivatives:
                hyperspherical_Hermite_interpolation_all(his, nx, i, &x[0],
                                                         &phi[i, 0], &dphi[i, 0], &d2phi[i, 0])
            else:
                hyperspherical_Hermite_interpolation_Phi(his, nx, i, &x[0],
                                                         &phi[i, 0], NULL, NULL)
    finally:
        del his
    return Phi, dPhi, d2Phi


def _bessel_wkb(int K, double beta, int[::1] l, double[::1] x):
    """Phi_l^beta(x) by the WKB (Airy) approximation, as an (n_l, n_x) array."""
    Phi = np.empty((l.shape[0], x.shape[0]), dtype=np.double)
    cdef double[:, ::1] out = Phi
    cdef int il, ix
    cdef double value
    for il in range(l.shape[0]):
        for ix in range(x.shape[0]):
            if not hyperspherical_WKB(K, l[il], beta, x[ix], &value):
                raise CosmoSevereError("hyperspherical_WKB failed")
            out[il, ix] = value
    return Phi


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
            self._pars['class_dir'] = _DATA_DIR
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
    cdef const perturbs* pt(self) except NULL:
        return &deref(deref(self._cosmo()).GetInputModule()).perturbations_

    cdef const nonlinear* nl(self) except NULL:
        return &deref(deref(self._cosmo()).GetInputModule()).nonlinear_

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

    cdef const PerturbationsModule* _perturbations(self) except NULL:
        return deref(self._cosmo()).GetPerturbationsModule().get()

    cdef _update_fc_from_pars(self):
        cdef:
            FileContent new_file_content
            char* dumc
            char* dumv

        # _fc will be cleaned up by its destructor
        self._fc = new_file_content

        # empty() and set() can leave _pars without class_dir; the C++ default would then be
        # the source tree this module was built from, which a wheel need not have.
        if 'class_dir' not in self._pars:
            self._fc.set(b'class_dir', _DATA_DIR.encode())
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
        if input_parameters.items() <= self._pars.items():
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

    cdef general_cl(self, int lmax, bool is_lensed):
        cdef:
            dict out_dict
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

        out_dict = {}
        for element in cl_data:
            out_dict[(<bytes> element.first).decode()] = _ndarray(element.second)
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
            out_dict[(<bytes> element.first).decode()] = _ndarray(element.second)
        out_dict['ell'] = out_dict['ell'].astype(int)
        return out_dict

    cpdef z_of_r (self, double[::1] z_array):
        cdef:
            Py_ssize_t i
            vector[double] pvecback
            double[::1] r
            double[::1] dzdr
            const BackgroundModule* background_module = self._background()

        r_arr = np.empty(z_array.shape[0], np.double)
        r = r_arr
        dzdr_arr = np.empty(z_array.shape[0], np.double)
        dzdr = dzdr_arr

        for i in range(z_array.shape[0]):
            pvecback = background_module.GetBackgroundAtZ(z_array[i])
            r[i] = pvecback[background_module.index_bg_conf_distance_]
            # dz/dr = H
            dzdr[i] = pvecback[background_module.index_bg_H_]

        return r_arr, dzdr_arr

    cpdef luminosity_distance(self, double z):
        """
        luminosity_distance(z)
        """
        return self.get_background_value_at_z(z, self._background().index_bg_lum_distance_)

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
        return _ndarray(self._perturbations().GetLateTimeRedshifts())

    def _pk_grid(self, nonlinear, only_clustering_species):
        """The stored P(k, z) grid, flat and z-major, and its k in 1/Mpc."""
        cdef vector[double] k_vector
        cdef vector[double] pk_vector
        cdef const NonlinearModule* nonlinear_module = self._nonlinear()
        index_pk = (nonlinear_module.index_pk_cluster_ if only_clustering_species
                    else nonlinear_module.index_pk_m_)
        nonlinear_module.GetPkGrid(pk_nonlinear if nonlinear else pk_linear, index_pk,
                                   k_vector, pk_vector)
        return _ndarray(pk_vector), _ndarray(k_vector)

    def _scalar_transfer_k(self):
        """The k axis of the scalar transfer grids, in 1/Mpc; refuses runs it does not fit."""
        return _ndarray(self._perturbations().GetScalarTransferK())

    def _transfer_at_node(self, int index_tau, output_format):
        """The transfer functions at a node of the stored late-time sampling."""
        return self.transfer_general(0., index_tau, output_format)

    def _sigma(self, double R, double z, bint cb):
        """sigma(R, z) of all matter, or of cdm + baryons with cb; R in Mpc."""
        cdef const NonlinearModule* nonlinear_module = self._nonlinear()
        return nonlinear_module.GetSigma(
            R, z, nonlinear_module.index_pk_cluster_ if cb else nonlinear_module.index_pk_m_)

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


    cdef double get_background_value_at_z(self, double z, int index_bg) except? -1.:
        return self._background().GetBackgroundAtZ(z)[index_bg]

    cdef double get_thermodynamics_value_at_z(self, double z, int index_th) except? -1.:
        return self._thermo().GetThermodynamicsAtZ(z)[index_th]

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
        return self.get_background_value_at_z(z, self._background().index_bg_ang_distance_)

    cpdef angular_distance_from_to(self, double z1, double z2):
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
        return self._background().GetAngularDistanceFromTo(z1, z2)

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
        return self.get_background_value_at_z(z, self._background().index_bg_D_)

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
        return self.get_background_value_at_z(z, self._background().index_bg_f_)

    cpdef z_of_tau(self, double tau):
        """
        Redshift corresponding to a given conformal time.

        Parameters
        ----------
        tau : float
                Conformal time
        """
        return self._background().GetZOfTau(tau)

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
        return self.get_background_value_at_z(z, self._background().index_bg_H_)

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
        return self.get_background_value_at_z(z, self._background().index_bg_Omega_m_)


    cpdef ionization_fraction(self, double z):
        """
        ionization_fraction(z)

        Return the ionization fraction for a given redshift z

        Parameters
        ----------
        z : float
                Desired redshift
        """
        return self.get_thermodynamics_value_at_z(z, self._thermo().index_th_xe_)

    cpdef baryon_temperature(self, double z):
        """
        baryon_temperature(z)

        Give the baryon temperature for a given redshift z

        Parameters
        ----------
        z : float
                Desired redshift
        """
        return self.get_thermodynamics_value_at_z(z, self._thermo().index_th_Tb_)

    cpdef get_background(self):
        """
        Return an array of the background quantities at all times.

        Parameters
        ----------

        Returns
        -------
        background : dictionary containing background.
        """
        return _columns_dict(self._background().GetOutputColumns())

    cpdef get_thermodynamics(self):
        """
        Return the thermodynamics quantities.

        Returns
        -------
        thermodynamics : dictionary containing thermodynamics.
        """
        return _columns_dict(self._thermo().GetOutputColumns())

    cpdef get_primordial(self):
        """
        Return the primordial scalar and/or tensor spectrum depending on 'modes'.
        'output' must be set to something, e.g. 'tCl'.

        Returns
        -------
        primordial : dictionary containing k-vector and primordial scalar and tensor P(k).
        """
        return _columns_dict(self._primordial().GetOutputColumns())


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

        cdef:
            const PerturbationsModule* perturbations_module
            vector[NamedColumns] tables
            size_t index_k

        perturbations = {}

        if self.pt().k_output_values_num<1:
            return perturbations

        perturbations_module = self._perturbations()
        for mode, has_mode, index_md in (
                ('scalar', self.pt().has_scalars, perturbations_module.index_md_scalars_),
                ('vector', self.pt().has_vectors, perturbations_module.index_md_vectors_),
                ('tensor', self.pt().has_tensors, perturbations_module.index_md_tensors_)):
            if not has_mode:
                continue
            tables = perturbations_module.GetPerturbationColumns(index_md)
            perturbations[mode] = [_columns_dict(tables[index_k]) for index_k in range(tables.size())]

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
            string ic_info
            string ic_suffix
            file_format outf
            const PerturbationsModule* perturbations_module
            vector[NamedColumns] tables
            int index_ic

        if (not self.pt().has_density_transfers) and (not self.pt().has_velocity_transfers):
            return {}

        if output_format == 'camb':
            outf = file_format.camb_format
        else:
            outf = file_format.class_format

        perturbations_module = self._perturbations()
        if index_tau < 0:
            tables = perturbations_module.GetTransferColumns(outf, z)
        else:
            tables = perturbations_module.GetTransferColumnsAtIndexTau(outf, index_tau)

        if tables.size() == 1:
            return _columns_dict(tables[0])

        transfers = {}
        for index_ic in range(tables.size()):
            perturbations_module.perturb_output_firstline_and_ic_suffix(index_ic, ic_info, ic_suffix)
            transfers[(<bytes> ic_suffix).decode()] = _columns_dict(tables[index_ic])
        return transfers

    def _derived_parameter(self, name):
        """One derived parameter, by name (Cosmology::GetDerivedParameter)."""
        return deref(self._cosmo()).GetDerivedParameter(name.encode())

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
