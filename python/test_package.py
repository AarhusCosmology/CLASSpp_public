"""classy as a package: the compiled core, the plain-Python Class on top of it, and the derived
parameters, which are one C++ lookup each
(docs/superpowers/specs/2026-11-20-classy-package-design.md)."""
import os

import numpy as np
import pytest

import classy
import classy._core

LCDM = {'output': 'tCl,mPk', 'P_k_max_h/Mpc': 1., 'compute damping scale': 'yes'}

# Each scalar getter of Class, and the derived parameter it returns.
GETTERS = {
    'age': 'age', 'h': 'h', 'n_s': 'n_s', 'A_s': 'A_s', 'tau_reio': 'tau_reio',
    'Omega_m': 'Omega_m', 'Omega0_m': 'Omega_m', 'Omega_r': 'Omega_r',
    'theta_s_100': '100*theta_s', 'theta_star_100': '100*theta_star',
    'theta_d_100': '100*theta_d', 'Omega_Lambda': 'Omega_Lambda', 'Omega_g': 'Omega_g',
    'Omega_b': 'Omega_b', 'omega_b': 'omega_b', 'Neff': 'Neff', 'f_idr_drmd': 'f_idr_drmd',
    'z_dec_drmd': 'z_dec_drmd', 'k_eq': 'k_eq', 'sigma8': 'sigma8', 'sigma8_cb': 'sigma8_cb',
    'rs_drag': 'rs_d', 'z_reio': 'z_reio', 'T_cmb': 'T_cmb', 'Omega0_k': 'Omega0_k',
    'Omega0_cdm': 'Omega0_cdm',
}


@pytest.fixture(scope='module')
def cosmo():
    c = classy.Class(LCDM)
    yield c
    c.struct_cleanup()


# ---------------------------------------------------------------- the package --

def test_class_is_a_subclass_of_the_compiled_core():
    # Inherited, so that a per-redshift call reaches C++ without a Python layer in between.
    assert issubclass(classy.Class, classy._core.PyCosmology)
    assert classy.Class.Hubble is classy._core.PyCosmology.Hubble


def test_exceptions_are_one_class_across_core_and_package():
    assert classy._core.CosmoSevereError is classy.CosmoSevereError
    assert classy._core.CosmoComputationError is classy.CosmoComputationError
    with pytest.raises(classy.CosmoSevereError):
        classy.Class({'no_such_input': 1}).compute(['background'])


def test_data_ships_inside_the_package():
    data = os.path.join(os.path.dirname(classy.__file__), 'data')
    assert os.path.isdir(os.path.join(data, 'hyrec'))
    assert os.path.isdir(os.path.join(data, 'bbn'))


def test_get_input_methods_are_gone():
    assert not hasattr(classy.Class, 'get_input_background')


# ----------------------------------------------------------------- parameters --
# The dict is Python's (_class.py); the core holds (name, value) strings. The gate itself
# (fresh results after set(), errors that stay until the parameters change) is pinned in
# test_class.py, TestLazyLifecycle.

def test_set_returns_self_or_none_when_nothing_changed():
    c = classy.Class({'h': 0.7, 'output': ''})
    assert c.set(h=0.7) is None
    assert c.set({'output': ''}, h=0.7) is None
    assert c.set(h=0.6) is c
    assert c.empty() is c and c.pars == {}


def test_pars_holds_what_was_set_and_no_more():
    # The class_dir default is the core's, applied when the Cosmology is built.
    c = classy.Class({'h': 0.7})
    c.set({'Omega_k': 0.}, output='')
    assert c.pars == {'h': 0.7, 'Omega_k': 0., 'output': ''}
    assert classy.Class().pars == {}


def test_a_name_that_is_not_a_str_changes_nothing():
    c = classy.Class({'h': 0.7, 'output': ''})
    with pytest.raises(TypeError):
        c.set({1: 0.6})
    assert c.pars == {'h': 0.7, 'output': ''}
    assert c.h() == 0.7


def test_set_builds_nothing_until_compute():
    # A sampler may call set() several times, with parameters that fit together only after
    # the last call. set() builds nothing, so nothing is checked or shot before compute() or
    # the first result asked for.
    c = classy.Class({'output': ''})
    assert c.set({'100*theta_s': 3.0}) is c  # a shooting would fail: no h gives this
    assert c.set(h=0.7) is c                 # and h with 100*theta_s is refused at the build
    c.empty()
    c.set({'100*theta_s': 1.0411, 'output': ''})
    assert c.compute(['background']) is c
    assert 0.6 < c.h() < 0.75


def test_compute_with_an_unknown_level_builds_nothing():
    # MontePython asks for class_sz's levels; the parameters below fail only when built.
    c = classy.Class({'h': 0.67, 'output': ''})
    c.set({'100*theta_s': 1.042})
    assert c.compute(['szcount']) is c
    with pytest.raises(classy.CosmoSevereError):
        c.compute(['background'])


# --------------------------------------------------------- derived parameters --

@pytest.mark.parametrize('getter,name', sorted(GETTERS.items()))
def test_every_getter_equals_its_derived_parameter(cosmo, getter, name):
    assert getattr(cosmo, getter)() == cosmo.get_current_derived_parameters([name])[name]


def test_Omega_nu_equals_its_derived_parameter(cosmo):
    assert cosmo.Omega_nu == cosmo.get_current_derived_parameters(['Omega_nu'])['Omega_nu']


def test_omega_m_is_Omega_m_times_h_squared(cosmo):
    # classy returned Omega_m/h^2 until 2026-11-20.
    d = cosmo.get_current_derived_parameters(['omega_m', 'Omega_m', 'h'])
    assert d['omega_m'] == pytest.approx(d['Omega_m']*d['h']**2, rel=1e-15)


def test_unknown_derived_parameter_is_severe(cosmo):
    with pytest.raises(classy.CosmoSevereError,
                       match='was not recognized as a derived parameter'):
        cosmo.get_current_derived_parameters(['no_such_parameter'])


def test_derived_parameter_follows_set():
    c = classy.Class({'h': 0.7, 'output': ''})
    assert c.get_current_derived_parameters(['h'])['h'] == 0.7
    c.set(h=0.6)
    assert c.get_current_derived_parameters(['h'])['h'] == 0.6
    assert c.h() == 0.6


def test_derived_names_may_be_numpy_strings(cosmo):
    d = cosmo.get_current_derived_parameters([np.str_('h')])
    assert list(d) == ['h'] and d['h'] == cosmo.h()


def test_non_list_names_raise_type_error(cosmo):
    with pytest.raises(TypeError):
        cosmo.get_current_derived_parameters(('h',))


# ----------------------------------------------------------------- CosmoHammer --

class _Data:
    def __init__(self, arguments):
        self.cosmo_arguments = arguments
        self.mcmc_parameters = {'h': {'current': arguments['h'], 'scale': 1.}}

    def get_mcmc_parameters(self, kinds):
        return ['h'] if kinds == ['derived'] else []


class _Context:
    def __init__(self, data):
        self.store = {'data': data}
        self.params = {}

    def get(self, key):
        return self.store[key]

    def getData(self):
        return self.params

    def add(self, key, value):
        self.store[key] = value


def test_cosmohammer_call_recomputes():
    c = classy.Class()
    for h in (0.7, 0.6):
        ctx = _Context(_Data({'h': h, 'output': ''}))
        c(ctx)
        assert ctx.store['cosmo'] is c and ctx.store['boundary'] is True
        assert c.h() == h and ctx.params['h'] == h


# --------------------------------------------------------- refusals, not crashes --

def test_nonlinear_scale_needs_a_nonlinear_method():
    # Without one, k_nl was never computed, and reading it crashed the interpreter.
    c = classy.Class({'output': 'mPk'})
    with pytest.raises(classy.CosmoSevereError, match='nonlinear method'):
        c.nonlinear_scale(np.array([0., 1.]), 2)


def test_raw_cl_with_two_initial_conditions():
    # The buffer for the cross-correlated initial conditions was sized only with several
    # modes, so scalars with ad + cdi wrote past it and crashed the interpreter.
    cl = classy.Class({'output': 'tCl,pCl', 'ic': 'ad,cdi', 'f_cdi': 0.5}).raw_cl(200)
    assert sorted(cl) == ['bb', 'ee', 'ell', 'te', 'tt']
    for key in ('tt', 'ee', 'te', 'bb'):
        assert np.all(np.isfinite(cl[key])) and np.all(cl[key][:2] == 0.)
    assert np.all(cl['tt'][2:] > 0.)


def test_a_failed_shooting_raises_instead_of_aborting():
    # The shooting runs on first use of the input; its exception escaped into CPython
    # and aborted the interpreter.
    c = classy.Class({'100*theta_s': 3.0, 'output': 'mPk'})
    with pytest.raises(classy.CosmoComputationError):
        c.nonlinear_method


def test_derived_parameters_obey_their_definitions(cosmo):
    # Relations that hold within one build. A name routed to the wrong field breaks them
    # (the rec, star and drag fields differ only in their redshift), which a comparison with
    # master can miss: ra_rec and ra_star are 1e-4 apart.
    names = ['z_rec', 'ra_rec', 'da_rec', 'rs_rec', 'ds_rec', 'z_star', 'ra_star', 'da_star',
             'rs_star', 'ds_star', 'z_d', 'rs_d', 'ds_d', '100*theta_s', '100*theta_star',
             'h', 'rs_rec_h', 'rs_d_h', 'A_s', 'ln10^{10}A_s', 'tau_reio', 'exp_m_2_tau_As']
    d = cosmo.get_current_derived_parameters(names)
    for when in ('rec', 'star'):
        z = d['z_' + when]
        assert d['ra_' + when] == pytest.approx(d['da_' + when]*(1 + z), rel=1e-12)
        assert d['rs_' + when] == pytest.approx(d['ds_' + when]*(1 + z), rel=1e-12)
    assert d['rs_d'] == pytest.approx(d['ds_d']*(1 + d['z_d']), rel=1e-12)
    assert d['100*theta_s'] == pytest.approx(100*d['rs_rec']/d['ra_rec'], rel=1e-12)
    assert d['100*theta_star'] == pytest.approx(100*d['rs_star']/d['ra_star'], rel=1e-12)
    assert d['rs_rec_h'] == pytest.approx(d['rs_rec']*d['h'], rel=1e-12)
    assert d['rs_d_h'] == pytest.approx(d['rs_d']*d['h'], rel=1e-12)
    assert d['ln10^{10}A_s'] == pytest.approx(np.log(1e10*d['A_s']), rel=1e-12)
    assert d['exp_m_2_tau_As'] == pytest.approx(np.exp(-2*d['tau_reio'])*d['A_s'], rel=1e-12)
    assert cosmo.rs_drag() == d['rs_d'] and cosmo.theta_s_100() == d['100*theta_s']


# ------------------------------------------------------------------ arguments --
# Typed in the core: nanobind converts and refuses, and the public flags are truth values.

def test_integer_arguments_take_numpy_integers_but_not_floats(cosmo):
    assert np.array_equal(cosmo.raw_cl(np.int64(100))['tt'], cosmo.raw_cl(100)['tt'])
    with pytest.raises(TypeError):
        cosmo.raw_cl(100.)
    with pytest.raises(TypeError):
        cosmo.get_pk_array(np.array([0.1]), np.zeros(1), -1, 1, False)


def test_a_size_beyond_its_array_is_a_value_error(cosmo):
    k, z = np.array([0.01, 0.1]), np.zeros(1)
    with pytest.raises(ValueError, match='k_size and z_size'):
        cosmo.get_pk_array(k, z, 3, 1, False)
    with pytest.raises(ValueError, match='mu_size'):
        cosmo.get_pk(np.full((2, 1, 1), 0.1), z, 2, 1, 2)
    assert cosmo.get_pk_array(k[:0], z, 0, 1, False).shape == (0,)


def test_flags_are_truth_values(cosmo):
    k, z = np.array([0.01, 0.1]), np.zeros(1)
    linear = cosmo.get_pk_array(k, z, 2, 1, False)
    for flag in (0, None, np.False_, []):
        assert np.array_equal(cosmo.get_pk_array(k, z, 2, 1, flag), linear)


def test_output_format_is_class_or_camb(cosmo):
    with pytest.raises(ValueError, match="'class' or 'camb', not 'CAMB'"):
        cosmo.get_transfer(0., 'CAMB')


# ------------------------------------- methods no caller in Cobaya or MontePython uses --
# Public classy API (class_public has them too), so scripts may: each against another way
# to get the same number.

def test_Om_m_scales_as_the_matter_density_over_H_squared(cosmo):
    # Omega_m(z) = Omega_m (1+z)^3 (H0/H)^2 with no massive neutrinos.
    H0 = cosmo.Hubble(0.)
    for z in (0., 0.5, 3., 1100.):
        expected = cosmo.Omega0_m()*(1. + z)**3*(H0/cosmo.Hubble(z))**2
        assert cosmo.Om_m(z) == pytest.approx(expected, rel=1e-6)


def test_z_of_tau_inverts_the_comoving_distance(cosmo):
    # Flat: tau(z) = tau_0 - r(z), with r from z_of_r, which is r(z) despite its name.
    tau_0 = cosmo.get_background()['conf. time [Mpc]'][-1]
    z = np.array([0.5, 3., 1100.])
    r, _ = cosmo.z_of_r(z)
    for z_i, r_i in zip(z, r):
        assert cosmo.z_of_tau(tau_0 - r_i) == pytest.approx(z_i, rel=1e-6)


def test_nonlinear_scale_cb():
    z = np.array([0., 1.])
    halofit = {'output': 'mPk', 'non_linear': 'halofit', 'P_k_max_1/Mpc': 3.}
    # Without massive neutrinos, the cdm+baryon spectrum is the matter spectrum.
    c = classy.Class(halofit)
    assert np.array_equal(c.nonlinear_scale_cb(z, 2), c.nonlinear_scale(z, 2))
    # With them, P_cb has more power than P_m, so it turns nonlinear on larger scales.
    c = classy.Class(dict(halofit, N_ur=0., N_ncdm=1, m_ncdm=0.3, deg_ncdm=3.))
    k_nl, k_nl_cb = c.nonlinear_scale(z, 2), c.nonlinear_scale_cb(z, 2)
    assert np.all(k_nl_cb > 0.) and np.all(k_nl_cb < k_nl)


def test_lensed_cl_computed_is_lensed_cl_at_the_computed_multipoles():
    c = classy.Class({'output': 'tCl,pCl,lCl', 'lensing': 'yes', 'l_max_scalars': 500})
    computed, lensed = c.lensed_cl_computed(), c.lensed_cl()
    ell = computed['ell']
    # The lensing module's own l grid, a sparse one, up to l_max_scalars.
    assert ell.dtype.kind == 'i' and ell[0] == 0 and ell[-1] <= lensed['ell'][-1]
    assert np.all(np.diff(ell) > 0) and len(ell) < len(lensed['ell'])/5
    assert sorted(computed) == sorted(lensed)
    for key in ('tt', 'ee', 'te', 'bb', 'pp', 'tp'):
        assert np.array_equal(computed[key], lensed[key][ell]), key
