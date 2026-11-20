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


def test_raw_cl_no_copy_with_two_initial_conditions():
    # Its buffer for the cross-correlated initial conditions was sized only with several
    # modes, so scalars with ad + cdi wrote past it and crashed the interpreter.
    c = classy.Class({'output': 'tCl,pCl', 'ic': 'ad,cdi', 'f_cdi': 0.5})
    without_copy, with_copy = c.raw_cl_no_copy(200), c.raw_cl(200)
    assert sorted(without_copy) == sorted(with_copy)
    for key in with_copy:
        assert np.array_equal(without_copy[key], with_copy[key])


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
