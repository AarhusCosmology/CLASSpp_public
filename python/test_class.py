"""
.. module:: test_class
    :synopsis: python script for testing CLASS using nose
.. moduleauthor:: Benjamin Audren <benjamin.audren@gmail.com>
.. credits:: Benjamin Audren, Thomas Tram
.. version:: 1.0

This is a python script for testing CLASS and its wrapper Classy using nose.
To run the test suite, type
nosetests test_class.py
If you want to extract the problematic input parameters at a later stage,
you should type
nosetests test_class.py 1>stdoutfile 2>stderrfile
and then use the python script extract_errors.py on the stderrfile.

When adding a new input parameter to CLASS (by modifying input.c), you
should also include tests of this new input. You will be in one of the
two cases:
1:  The new input is supposed to be compatible with any existing input.
    This is the standard case when adding a new species for instance.
2:  The new input is incompatible with one of the existing inputs. This
    would be the case if you have added (or just want to test) some other
    value of an already defined parameter. (Maybe you have allowed for
    negative mass neutrinos and you want to test CLASS using a negative mass.)

In case 1, you must add an entry in the CLASS_INPUT dictionary:
CLASS_INPUT['Mnu'] = (
    [{'N_eff': 0.0, 'N_ncdm': 1, 'm_ncdm': 0.06, 'deg_ncdm': 3.0},
     {'N_eff': 1.5, 'N_ncdm': 1, 'm_ncdm': 0.03, 'deg_ncdm': 1.5}],
    'normal')
The key 'Mnu' is not being used in the code, so its purpose is just to
describe the entry to the reader.
the value is a 2-tuple where the first entry [{},{},...,{}] is an array of
dictionaries containg the actual input to CLASS. The second entry is a keyword
which can be either 'normal' or 'power'. It tells the script how this input
will be combined with other inputs.

What does 'normal' and 'power' mean?
If an entry has the 'power' keyword, it will be combined with any other entry.
If an entry has the 'normal' keyword, it will not be combined with any other
entry having the 'normal' keyword, but it will be combined with all entries
carrying the 'power keyword.
Beware that the number of tests grow a lot when using the 'power' keyword.

In case 2, you should find the relevant entry and just add a new dictionary
to the array. E.g. if you want to test some negative mass model you should add
{'N_ncdm': 1, 'm_ncdm': -0.1, 'deg_ncdm': 1.0}

How are default parameters handled?
Any input array implicitly contains the empty dictionary. That means that if
Omega_k:0.0 is the default value, writing
CLASS_INPUT['Curvature'] = (
    [{'Omega_k': 0.01},
     {'Omega_k': -0.01}],
    'normal')
will test the default value Omega_k=0.0 along with the two specified models.

How to deal with inconsistent input?
Sometimes a specific feature requires the presence of another input parameter.
For instance, if we ask for tensor modes we must have temperature and/or
polarisation in the output. If not, CLASS is supposed to fail during the
evaluation of the input module and return an error message. This fail is the
correct behaviour of CLASS. To implement such a case, modify the function
has_incompatible_input(self)

Comparing output: When the flag 'COMPARE_OUTPUT_GAUGE' is set to true, the code will
rerun CLASS for each case under Newtonian gauge and then compare Cl's and
matter power spectrum. If the two are not close enough, it will generate a
PDF plot of this and save it in the 'fail' folder.
"""
import matplotlib as mpl
mpl.use('Agg')

import itertools
import matplotlib.pyplot as plt
import numpy as np
import os
import pytest
import re
import shutil
import tempfile
import unittest

from classy import Class
from classy import CosmoSevereError
from classy import CosmoComputationError
from math import log10
from matplotlib.offsetbox import AnchoredText
from parameterized import parameterized

# Customise test by reading environment variables
CLASS_VERBOSE = bool(int(os.getenv('CLASS_VERBOSE', '0'))) # Print output from CLASS?
COMPARE_OUTPUT_GAUGE = bool(int(os.getenv('COMPARE_OUTPUT_GAUGE', '0'))) # Compare synchronous and Newtonian gauge outputs?
COMPARE_OUTPUT_REF = bool(int(os.getenv('COMPARE_OUTPUT_REF', '0'))) # Compare classy with classyref?
POWER_ALL = bool(int(os.getenv('POWER_ALL', '0'))) # Combine every extension with each other? (Very slow!)
TEST_LEVEL = int(os.getenv('TEST_LEVEL', '0')) # 0 <= TEST_LEVEL <= 3

if COMPARE_OUTPUT_REF:
    try:
        import classyref
    except ImportError as exc:
        raise ImportError(
            "COMPARE_OUTPUT_REF=1 was requested, but classyref could not be imported. "
            "Install the reference wrapper as classyref before running reference comparisons."
        ) from exc
# Define bounds on the relative and absolute errors of C(l) and P(k)
# between reference, Newtonian and Synchronous gauge
COMPARE_CL_RELATIVE_ERROR = 3e-3
COMPARE_CL_RELATIVE_ERROR_GAUGE = 5*3e-3
COMPARE_CL_ABSOLUTE_ERROR = 1e-20
COMPARE_PK_RELATIVE_ERROR = 1e-2
COMPARE_PK_RELATIVE_ERROR_GAUGE = 5*1e-2
COMPARE_PK_ABSOLUTE_ERROR = 1e-20

# Dictionary of models to test the wrapper against. Each of these scenario will
# be run against all the possible output choices (nothing, tCl, mPk, etc...),
# with or without non-linearities.
# Never input the default value, as this will be **automatically** tested
# against. Indeed, when not specifying a field, CLASS takes the default input.
CLASS_INPUT = {}

CLASS_INPUT['Output_spectra'] = (
    [{'output': 'mPk', 'P_k_max_1/Mpc': 2},
     {'output': 'tCl'},
     {'output': 'tCl pCl lCl'},
     {'output': 'mPk tCl lCl', 'P_k_max_1/Mpc': 2},
     {'output': 'nCl sCl'},
     {'output': 'tCl pCl lCl nCl sCl'}],
    'power')

CLASS_INPUT['Nonlinear'] = (
    [{'non linear': 'halofit'}],
    'power')

CLASS_INPUT['Lensing'] = (
    [{'lensing': 'yes'}],
    'power')

if TEST_LEVEL > 0:
    CLASS_INPUT['Mnu'] = (
        [{'N_ur': 0.0, 'N_ncdm': 1, 'm_ncdm': 0.06, 'deg_ncdm': 3.0},
        {'N_ur': 1.5, 'N_ncdm': 1, 'm_ncdm': 0.03, 'deg_ncdm': 1.5}],
        'normal')

if TEST_LEVEL > 1:
    CLASS_INPUT['Curvature'] = (
        [{'Omega_k': 0.01},
        {'Omega_k': -0.01}],
        'normal')

    CLASS_INPUT['modes'] = (
        [{'modes': 't'},
        {'modes': 's, t'}],
        'power')

    CLASS_INPUT['Tensor_method'] = (
        [{'tensor method': 'exact'},
        {'tensor method': 'photons'}],
        'power')

    CLASS_INPUT['Scalar_field'] = (
        [{'Omega_scf': 0.1, 'attractor_ic_scf': 'yes',
        'scf_parameters': '10, 0, 0, 0'}],
        'normal')

if TEST_LEVEL > 2:
    CLASS_INPUT['Isocurvature_modes'] = (
        [{'ic': 'ad,nid,cdi', 'c_ad_cdi': -0.5}],
        'normal')

    CLASS_INPUT['Inflation'] = (
        [{'P_k_ini type': 'inflation_V'},
        {'P_k_ini type': 'inflation_H'}],
        'normal')
    # CLASS_INPUT['Inflation'] = (
    #     [{'P_k_ini type': 'inflation_V'},
    #     {'P_k_ini type': 'inflation_H'},
    #     {'P_k_ini type': 'inflation_V_end'}],
    #     'normal')

if POWER_ALL:
    for key, (models, state) in CLASS_INPUT.items():
        CLASS_INPUT[key] = (models, 'power')

INPUTPOWER = []
INPUTNORMAL = [{}]
for key, (models, state) in CLASS_INPUT.items():
    if state == 'power':
        INPUTPOWER.append([{}]+models)
    else:
        INPUTNORMAL.extend(models)

    PRODPOWER = list(itertools.product(*INPUTPOWER))

    DICTARRAY = []
    for normelem in INPUTNORMAL:
        for powelem in PRODPOWER:  # itertools.product(*modpower):
            temp_dict = normelem.copy()
            for elem in powelem:
                temp_dict.update(elem)
            DICTARRAY.append(temp_dict)


TUPLE_ARRAY = []
for e in DICTARRAY:
    TUPLE_ARRAY.append((e, ))


def powerset(iterable):
    xs = list(iterable)
    # note we return an iterator rather than a list
    return itertools.chain.from_iterable(
        itertools.combinations(xs, n) for n in range(1, len(xs)+1))

def custom_name_func(testcase_func, param_num, param):
    special_keys = ['N_ncdm']
    somekeys = []
    for key in param.args[0].keys():
        if key in special_keys:
            somekeys.append(key)
        elif 'mega' in key:
            somekeys.append(key)

    res = '{}_{:04d}_{}'.format(
        testcase_func.__name__,
        param_num,
        parameterized.to_safe_name('_'.join(somekeys))
    )
    return res.strip('_')

class TestClass(unittest.TestCase):
    """
    Testing Class and its wrapper classy on different cosmologies

    To run it, do
    ~] nosetest test_class.py

    It will run many times Class, on different cosmological scenarios, and
    everytime testing for different output possibilities (none asked, only mPk,
    etc..)

    """
    @classmethod
    def setUpClass(cls):
        cls.faulty_figs_path = os.path.join(
            os.path.sep.join(os.path.realpath(__file__).split(
                os.path.sep)[:-1]),
            'faulty_figs')

        if os.path.isdir(cls.faulty_figs_path):
            shutil.rmtree(cls.faulty_figs_path)

        os.mkdir(cls.faulty_figs_path)

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        """
        set up data used in the tests.
        setUp is called before each test function execution.
        """
        self.cosmo = Class()
        self.cosmo_newt = Class()

        if CLASS_VERBOSE:
            self.verbose = {
                'input_verbose': 1,
                'background_verbose': 1,
                'thermodynamics_verbose': 1,
                'perturbations_verbose': 1,
                'transfer_verbose': 1,
                'primordial_verbose': 1,
                'spectra_verbose': 1,
                'nonlinear_verbose': 1,
                'lensing_verbose': 1,
                'output_verbose': 1,
            }
        else:
            self.verbose = {}
        self.scenario = {}

    def tearDown(self):
        self.cosmo.struct_cleanup()
        self.cosmo.empty()
        self.cosmo = 0
        self.cosmo_newt.struct_cleanup()
        self.cosmo_newt.empty()
        self.cosmo_newt = 0
        del self.scenario

    def poormansname(self, somedict):
        string = "_".join(
            [k+'='+str(v)
             for k, v in list(somedict.items())])
        string = string.replace('/', '%')
        string = string.replace(',', '')
        string = string.replace(' ', '')
        return string
    
    def has_incompatible_input(self):

        should_fail = False

        # If we have tensor modes, we must have one tensor observable,
        # either tCl or pCl.
        if has_tensor(self.scenario):
            if 'output' not in list(self.scenario.keys()):
                should_fail = True
            else:
                output = self.scenario['output'].split()
                if 'tCl' not in output and 'pCl' not in output:
                    should_fail = True

        # If we have specified lensing, we must have lCl in output,
        # otherwise lensing will not be read (which is an error).
        if 'lensing' in list(self.scenario.keys()):
            if 'output' not in list(self.scenario.keys()):
                should_fail = True
            else:
                output = self.scenario['output'].split()
                if 'lCl' not in output:
                    should_fail = True
                elif 'tCl' not in output and 'pCl' not in output:
                    should_fail = True

        # If we have specified a tensor method, we must have tensors.
        if 'tensor method' in list(self.scenario.keys()):
            if not has_tensor(self.scenario):
                should_fail = True

        # If we have specified non linear, we must have some form of
        # perturbations output.
        if 'non linear' in list(self.scenario.keys()):
            if 'output' not in list(self.scenario.keys()):
                should_fail = True

        # If we ask for Cl's of lensing potential, we must have scalar modes.
        if 'output' in list(self.scenario.keys()) and 'lCl' in self.scenario['output'].split():
            if 'modes' in list(self.scenario.keys()) and self.scenario['modes'].find('s') == -1:
                should_fail = True

        # If we specify initial conditions (for scalar modes), we must have
        # perturbations and scalar modes.
        if 'ic' in list(self.scenario.keys()):
            if 'modes' in list(self.scenario.keys()) and self.scenario['modes'].find('s') == -1:
                should_fail = True
            if 'output' not in list(self.scenario.keys()):
                should_fail = True

        # If we use inflation module, we must have scalar modes,
        # tensor modes, no vector modes and we should only have adiabatic IC:
        if 'P_k_ini type' in list(self.scenario.keys()) and self.scenario['P_k_ini type'].find('inflation') != -1:
            if 'modes' not in list(self.scenario.keys()):
                should_fail = True
            else:
                if self.scenario['modes'].find('s') == -1:
                    should_fail = True
                if self.scenario['modes'].find('v') != -1:
                    should_fail = True
                if self.scenario['modes'].find('t') == -1:
                    should_fail = True
            if 'ic' in list(self.scenario.keys()) and self.scenario['ic'].find('i') != -1:
                should_fail = True


        return should_fail

    def compare_output(self, reference, reference_name, candidate, candidate_name, rtol_cl, rtol_pk):
        status_pass = True
        for elem in ['raw_cl', 'lensed_cl']:
            # Try to get the elem, but if they were not computed, a
            # CosmoComputeError should be raised. In this case, ignore the
            # whole block.
            try:
                to_test = getattr(candidate, elem)()
            except CosmoSevereError:
                continue
            ref = getattr(reference, elem)()
            for key, value in list(ref.items()):
                if key != 'ell':
                    # For all self spectra, try to compare allclose
                    if key[0] == key[1]:
                        # If it is a 'dd' or 'll', it is a dictionary.
                        if isinstance(value, dict):
                            for subkey in list(value.keys()):
                                try:
                                    np.testing.assert_allclose(
                                        value[subkey],
                                        to_test[key][subkey],
                                        rtol=rtol_cl,
                                        atol=COMPARE_CL_ABSOLUTE_ERROR)
                                except AssertionError:
                                    self.cl_faulty_plot(elem + "_" + key, value[subkey][2:], reference_name, to_test[key][subkey][2:], candidate_name, rtol_cl)
                                except TypeError:
                                    self.cl_faulty_plot(elem + "_" + key, value[subkey][2:], reference_name, to_test[key][subkey][2:], candidate_name, rtol_cl)
                        else:
                            try:
                                np.testing.assert_allclose(
                                    value,
                                    to_test[key],
                                    rtol=rtol_cl,
                                    atol=COMPARE_CL_ABSOLUTE_ERROR)
                            except (AssertionError, TypeError) as e:
                                self.cl_faulty_plot(elem + "_" + key, value[2:], reference_name, to_test[key][2:], candidate_name, rtol_cl)
                                status_pass = False
                    # For cross-spectra, as there can be zero-crossing, we
                    # compare an absolute difference against a physical scale.
                    # Pointwise relative errors are not meaningful near a zero
                    # crossing. For TE, sqrt(TT*EE) is the natural scale.
                    else:
                        try:
                            if key == 'te':
                                norm = np.maximum(
                                    np.sqrt(np.abs(ref['tt'] * ref['ee'])),
                                    np.sqrt(np.abs(
                                        to_test['tt'] * to_test['ee'])))
                                np.testing.assert_array_less(
                                    np.abs(value - to_test[key]),
                                    rtol_cl * norm + COMPARE_CL_ABSOLUTE_ERROR)
                            else:
                                norm = max(
                                    np.abs(value).max(),
                                    np.abs(to_test[key]).max())
                                if norm == 0.0:
                                    np.testing.assert_array_equal(
                                        value, to_test[key])
                                else:
                                    np.testing.assert_allclose(
                                        value / norm,
                                        to_test[key] / norm,
                                        rtol=0.0,
                                        atol=rtol_cl)
                        except AssertionError:
                            self.cl_faulty_plot(elem + "_" + key, value[2:], reference_name, to_test[key][2:], candidate_name, rtol_cl)
                            status_pass = False

        if 'output' in list(self.scenario.keys()):
            if self.scenario['output'].find('mPk') != -1:
                # testing equality of Pk
                k = np.logspace(-2, log10(self.scenario['P_k_max_1/Mpc']), 50)
                reference_pk = np.array([reference.pk(elem, 0) for elem in k])
                candidate_pk = np.array([candidate.pk(elem, 0) for elem in k])
                try:
                    np.testing.assert_allclose(
                        reference_pk,
                        candidate_pk,
                        rtol=rtol_pk,
                        atol=COMPARE_PK_ABSOLUTE_ERROR)
                except AssertionError:
                    self.pk_faulty_plot(k, reference_pk, reference_name, candidate_pk, candidate_name, rtol_pk)
                    status_pass = False

        return status_pass

    def store_ini_file(self, path):
        parameters = dict(self.verbose, **self.scenario)
        with open(path + '.ini', 'w') as param_file:
            param_file.write('# ' + str(parameters) + '\n')
            if len(parameters) == 0:
                # CLASS complains if the .ini file does not do anything.
                param_file.write('write warnings = yes\n')
            for key, value in list(parameters.items()):
                param_file.write(key + " = " + str(value)+ '\n')

    def cl_faulty_plot(self, cl_type, reference, reference_name, candidate, candidate_name, rtol):
        path = os.path.join(self.faulty_figs_path, self.name)
        fig, axes = plt.subplots(2, 1, sharex=True)
        ell = np.arange(max(np.shape(candidate))) + 2
        factor = ell*(ell + 1)/(2*np.pi) if cl_type[-2:] != 'pp' else ell**5
        axes[0].plot(ell, factor*reference, label=reference_name)
        axes[0].plot(ell, factor*candidate, label=candidate_name)
        axes[1].semilogy(ell, 100*abs(candidate/reference - 1), label=cl_type)
        axes[1].axhline(y=100*rtol, color='k', ls='--')

        axes[-1].set_xlabel(r'$\ell$')
        if cl_type[-2:] == 'pp':
            axes[0].set_ylabel(r'$\ell^5 C_\ell^\mathrm{{{_cl_type}}}$'.format(_cl_type=cl_type[-2:].upper()))
        else:
            axes[0].set_ylabel(r'$\ell(\ell + 1)/(2\pi)C_\ell^\mathrm{{{_cl_type}}}$'.format(_cl_type=cl_type[-2:].upper()))
        axes[1].set_ylabel('Relative error [%]')

        for ax in axes:
            ax.legend(loc='upper right')

        fig.tight_layout()
        fname = '{}_{}_{}_vs_{}.pdf'.format(path, cl_type, reference_name, candidate_name)
        fig.savefig(fname, bbox_inches='tight')
        plt.close(fig)

        # Store parameters (contained in self.scenario) to text file
        self.store_ini_file(path)

    def pk_faulty_plot(self, k, reference, reference_name, candidate, candidate_name, rtol):
        path = os.path.join(self.faulty_figs_path, self.name)

        fig, axes = plt.subplots(2, 1, sharex=True)
        axes[0].loglog(k, k**1.5*reference, label=reference_name)
        axes[0].loglog(k, k**1.5*candidate, label=candidate_name)
        axes[0].legend(loc='upper right')

        axes[1].loglog(k, 100*np.abs(candidate/reference - 1))
        axes[1].axhline(y=100*rtol, color='k', ls='--')

        axes[-1].set_xlabel(r'$k\quad [\mathrm{Mpc}^{-1}]$')
        axes[0].set_ylabel(r'$k^\frac{3}{2}P(k)$')
        axes[1].set_ylabel(r'Relative error [%]')

        fig.tight_layout()
        fname = path + '_pk_{}_vs_{}.pdf'.format(reference_name, candidate_name)
        fig.savefig(fname, bbox_inches='tight')
        plt.close(fig)

        # Store parameters (contained in self.scenario) to text file
        self.store_ini_file(path)

def has_tensor(input_dict):
    if 'modes' in list(input_dict.keys()):
        if input_dict['modes'].find('t') != -1:
            return True
    else:
        return False
    return False

@pytest.mark.dump_ini_files
class DumpIniFiles(TestClass):
    @parameterized.expand(TUPLE_ARRAY, doc_func=custom_name_func, custom_name_func=custom_name_func)
    def test_Valgrind(self, inputdict):
        """Dump files"""
        self.scenario.update(inputdict)
        self.name = self._testMethodName
        if self.has_incompatible_input():
            return
        path = os.path.join(self.faulty_figs_path, self.name)
        self.store_ini_file(path)
        self.scenario.update({'gauge':'Newtonian'})
        self.store_ini_file(path + 'N')


@pytest.mark.test_scenario
class TestScenario(TestClass):
    @parameterized.expand(TUPLE_ARRAY, doc_func=custom_name_func, custom_name_func=custom_name_func)
    def test_scenario(self, inputdict):
        """Test scenario"""
        self.scenario.update(inputdict)
        self.name = self._testMethodName
        self.cosmo.set(dict(itertools.chain(self.verbose.items(), self.scenario.items())))

        cl_dict = {
            'tCl': ['tt'],
            'lCl': ['pp'],
            'pCl': ['ee', 'bb'],
            'nCl': ['dens[1]-dens[1]'],
            'sCl': ['lens[1]-lens[1]'],
        }

        # 'lensing' is always set to yes. Therefore, trying to compute 'tCl' or
        # 'pCl' will fail except if we also ask for 'lCl'.
        if self.has_incompatible_input():
            self.assertRaises(CosmoSevereError, self.cosmo.compute)
            return
        else:
            self.cosmo.compute()

        self.assertTrue(
            self.cosmo.state,
            "Class failed to go through all __init__ methods")
        # Depending
        if 'output' in self.scenario.keys():
            # Positive tests of raw cls
            output = self.scenario['output']
            for elem in output.split():
                if elem in cl_dict.keys():
                    for cl_type in cl_dict[elem]:
                        cl = self.cosmo.raw_cl(100)
                        self.assertIsNotNone(cl, "raw_cl returned nothing")
                        self.assertEqual(
                            np.shape(cl[cl_type])[0], 101,
                            "raw_cl returned wrong size")
                if elem == 'mPk':
                    pk = self.cosmo.pk(0.1, 0)
                    self.assertIsNotNone(pk, "pk returned nothing")
            # Negative tests of output functions
            if not any([elem in list(cl_dict.keys()) for elem in output.split()]):
                # testing absence of any Cl
                self.assertRaises(CosmoSevereError, self.cosmo.raw_cl, 100)
            if 'mPk' not in output.split():
                # testing absence of mPk
                self.assertRaises(CosmoSevereError, self.cosmo.pk, 0.1, 0)

        if COMPARE_OUTPUT_REF or COMPARE_OUTPUT_GAUGE:
            # Now compute same scenario in Newtonian gauge
            self.cosmo_newt.set(dict(self.verbose, **self.scenario))
            self.cosmo_newt.set({'gauge': 'newtonian'})
            self.cosmo_newt.compute()

        if COMPARE_OUTPUT_GAUGE:
            # Compare synchronous and Newtonian gauge
            self.assertTrue(
                self.cosmo_newt.state,
                "Class failed to go through all __init__ methods in Newtonian gauge")

            status = self.compare_output(self.cosmo, "Synchronous", self.cosmo_newt, 'Newtonian', COMPARE_CL_RELATIVE_ERROR_GAUGE, COMPARE_PK_RELATIVE_ERROR_GAUGE)
            assert status, 'Gauge comparison failed!'

        if COMPARE_OUTPUT_REF:
            # Compute reference models in both gauges and compare
            cosmo_ref = classyref.Class()
            cosmo_ref.set(dict(self.verbose, **self.scenario))
            cosmo_ref.compute()
            status = self.compare_output(cosmo_ref, "Reference", self.cosmo, 'Synchronous', COMPARE_CL_RELATIVE_ERROR, COMPARE_PK_RELATIVE_ERROR)
            assert status, 'Reference comparison failed in Synchronous gauge!'

            # Mainline CLASS has the historical incomplete Newtonian scalar
            # KG equation, so it is not a valid Newtonian reference for scf.
            if 'Omega_scf' not in self.scenario:
                cosmo_ref = classyref.Class()
                cosmo_ref.set(dict(self.verbose, **self.scenario))
                cosmo_ref.set({'gauge': 'newtonian'})
                cosmo_ref.compute()
                status = self.compare_output(cosmo_ref, "Reference", self.cosmo_newt, 'Newtonian', COMPARE_CL_RELATIVE_ERROR, COMPARE_PK_RELATIVE_ERROR)
                assert status, 'Reference comparison failed in Newtonian gauge!'


class TestTensorMassiveNcdmRegression(TestClass):
    def _assert_reference_match(self, case_name, scenario):
        try:
            import classyref
        except ImportError:
            self.skipTest("classyref not available")

        self.scenario = dict(scenario)
        self.name = f"{self._testMethodName}_{case_name}"

        candidate = Class()
        reference = classyref.Class()
        try:
            candidate.set(dict(self.verbose, **scenario))
            candidate.compute()
            reference.set(dict(self.verbose, **scenario))
            reference.compute()
            status = self.compare_output(
                reference,
                "Reference",
                candidate,
                "Candidate",
                COMPARE_CL_RELATIVE_ERROR,
                COMPARE_PK_RELATIVE_ERROR)
            self.assertTrue(status, f"Reference comparison failed for {case_name}")
        finally:
            reference.struct_cleanup()
            reference.empty()
            candidate.struct_cleanup()
            candidate.empty()

    def test_tensor_massive_ncdm_matches_reference(self):
        base = {
            'N_ur': 0.0,
            'N_ncdm': 1,
            'm_ncdm': 0.06,
            'deg_ncdm': 3.0,
            'output': 'tCl',
            'modes': 't',
        }

        cases = [
            ('default', {}),
            ('exact', {'tensor method': 'exact'}),
        ]

        for case_name, extra in cases:
            with self.subTest(case_name=case_name):
                self._assert_reference_match(case_name, dict(base, **extra))


class TestReviewRegressions(TestClass):
    def _dot_syntax_base(self):
        return {
            'h': 0.67556,
            'omega_b': 0.022032,
            'omega_cdm': 0.12038,
            'A_s': 2.215e-9,
            'n_s': 0.9619,
            'tau_reio': 0.054,
            'YHe': 0.25,
            'output': 'tCl',
            'l_max_scalars': 100,
            'N_ur': 0.0,
        }

    def _assert_compute_succeeds(self, scenario):
        self.scenario = dict(scenario)
        self.name = self._testMethodName
        self.cosmo.set(dict(self.verbose, **scenario))
        self.cosmo.compute()
        self.assertTrue(self.cosmo.state)

    def _assert_compute_fails(self, scenario, message):
        self.scenario = dict(scenario)
        self.name = self._testMethodName
        self.cosmo.set(dict(self.verbose, **scenario))
        with self.assertRaises(CosmoSevereError) as ctx:
            self.cosmo.compute()
        self.assertIn(message, str(ctx.exception))

    def _assert_scenarios_match(self, scenario, reference_scenario, reference_name):
        self.scenario = dict(reference_scenario)
        self.name = self._testMethodName

        candidate = Class()
        reference = Class()
        try:
            candidate.set(dict(self.verbose, **scenario))
            candidate.compute()
            reference.set(dict(self.verbose, **reference_scenario))
            reference.compute()
            status = self.compare_output(
                reference,
                reference_name,
                candidate,
                "Candidate",
                COMPARE_CL_RELATIVE_ERROR,
                COMPARE_PK_RELATIVE_ERROR)
            self.assertTrue(status, f"Scenario mismatch against {reference_name}")
        finally:
            reference.struct_cleanup()
            reference.empty()
            candidate.struct_cleanup()
            candidate.empty()

    # The classyref comparisons below guard refactors, whose mistakes are O(1). A numerical
    # change on master moves them by more than rounding: the spline integrals of
    # tools/arrays.cpp moved rs_drag by 2e-6 and delta_idr by 3e-6 (2026-10-08).
    REFERENCE_RTOL = 1e-5

    def _compute_candidate_and_reference(self, scenario):
        try:
            import classyref
        except ImportError:
            self.skipTest("classyref not available")

        self.scenario = dict(scenario)
        self.name = self._testMethodName

        candidate = Class()
        reference = classyref.Class()
        try:
            candidate.set(dict(self.verbose, **scenario))
            candidate.compute()
            reference.set(dict(self.verbose, **scenario))
            reference.compute()
            return candidate, reference
        except Exception:
            reference.struct_cleanup()
            reference.empty()
            candidate.struct_cleanup()
            candidate.empty()
            raise

    def test_rs_drag_matches_reference(self):
        candidate, reference = self._compute_candidate_and_reference({
            'compute damping scale': 'yes',
        })
        try:
            self.assertAlmostEqual(candidate.rs_drag() / reference.rs_drag(), 1.,
                                   delta=self.REFERENCE_RTOL)
        finally:
            reference.struct_cleanup()
            reference.empty()
            candidate.struct_cleanup()
            candidate.empty()

    def test_dcdm_dr_matches_reference(self):
        scenario = {
            'Omega_dcdmdr': 0.12,
            'Gamma_dcdm': 10.0,
            'output': 'tCl',
            'l_max_scalars': 200,
        }
        candidate, reference = self._compute_candidate_and_reference(scenario)
        try:
            status = self.compare_output(
                reference,
                "Reference",
                candidate,
                "Candidate",
                COMPARE_CL_RELATIVE_ERROR,
                COMPARE_PK_RELATIVE_ERROR)
            self.assertTrue(status, "Reference comparison failed for DCDM_DR scenario")
        finally:
            reference.struct_cleanup()
            reference.empty()
            candidate.struct_cleanup()
            candidate.empty()

    def test_dncdm_dr_computes(self):
        """#308 guard: a decaying-NCDM (DNCDM_DR) run still computes after the
        relativistic-IC check is forwarded to the wrapped DNCDM child."""
        cosmo = Class()
        try:
            cosmo.set({
                'output': 'tCl',
                'l_max_scalars': 200,
                'N_ur': 3.046,
                'omega_b': 0.022032,
                'omega_cdm': 0.12038,
                # Pin YHe to skip BBN; this is a compute guard, not a precision test.
                'YHe': 0.25,
                'dncdm1.type': 'ncdm_decay_dr',
                'dncdm1.m': 1.0,
                'dncdm1.T': 0.71611,
                'dncdm1.Gamma': 1e3,
                'dncdm1.Omega_ini': 0.001,
                # A decaying parent needs one momentum grid for background and perturbations.
                'dncdm1.quadrature_strategy': 3,
                'dncdm1.momenta_bins': 8,
            })
            cosmo.compute()
            self.assertTrue(cosmo.raw_cl(100)['tt'].size > 0)
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()

    def test_dncdm_dr_collision_lmax_default_tracks_the_hierarchies(self):
        """l_max_dr_col has no default of its own: the DR collision at multipole
        l reads the parent NCDM's l-th multipole, so perturb_init refuses
        l_max_dr_col > l_max_ncdm. Its literal default of 17 was consistent
        until the #397 precision defaults lowered l_max_ncdm to 10, after which
        the DEFAULT configuration hit that check and every decaying-NCDM run
        aborted (test_dncdm_dr_computes above is the same regression). The
        default is now min(l_max_dr, l_max_ncdm); an explicit value is still
        validated, not clamped. l_max_ncdm is set to 10 here so that the minimum
        binds whatever the compiled defaults are (17 again since 2026-10-09)."""
        base = {
            'output': 'tCl',
            'l_max_scalars': 200,
            'N_ur': 3.046,
            'omega_b': 0.022032,
            'omega_cdm': 0.12038,
            'YHe': 0.25,
            'dncdm1.type': 'ncdm_decay_dr',
            'dncdm1.m': 1.0,
            'dncdm1.T': 0.71611,
            'dncdm1.Gamma': 1e3,
            'dncdm1.Omega_ini': 0.001,
            # A decaying parent needs one momentum grid for background and perturbations.
            'dncdm1.quadrature_strategy': 3,
            'dncdm1.momenta_bins': 8,
            'l_max_ncdm': 10,
        }
        self.scenario = dict(base)
        spectra = {}
        for label, extra in (('derived', {}), ('explicit', {'l_max_dr_col': 10})):
            cosmo = Class()
            cosmo.set(dict(self.verbose, **base, **extra))
            try:
                cosmo.compute()
                spectra[label] = cosmo.raw_cl(100)['tt'].copy()
            finally:
                cosmo.struct_cleanup()
                cosmo.empty()
        # The derived default must BE min(l_max_dr, l_max_ncdm) = 10, not merely
        # something that runs.
        np.testing.assert_array_equal(spectra['derived'], spectra['explicit'])

        # An explicit value beyond the hierarchies still fails loudly: the
        # collision would read multipoles the parent does not have.
        cosmo = Class()
        cosmo.set(dict(self.verbose, **base, **{'l_max_dr_col': 30}))
        try:
            with self.assertRaises(CosmoComputationError):
                cosmo.compute()
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()

    def test_tensor_with_dncdm_dr_is_rejected(self):
        """#345 guard: evolving tensor modes with the exact NCDM method
        together with a decaying-NCDM (DNCDM_DR) species is not supported and
        must be rejected up front. The DNCDM_DR composite is a CompositeSpecies
        (not an NCDMSpecies), so the tensor-NCDM loops would silently skip it;
        the guard in perturb_init prevents that silent, unsupported run. It
        fires during perturbation solving, so classy raises a
        CosmoComputationError (not the input-time CosmoSevereError)."""
        scenario = {
            'output': 'tCl',
            'modes': 't',
            'tensor method': 'exact',
            'N_ur': 3.046,
            'omega_b': 0.022032,
            'omega_cdm': 0.12038,
            # Pin YHe to skip BBN; this is a guard test, not a precision test.
            'YHe': 0.25,
            'dncdm1.type': 'ncdm_decay_dr',
            'dncdm1.m': 1.0,
            'dncdm1.T': 0.71611,
            'dncdm1.Gamma': 1e3,
            'dncdm1.Omega_ini': 0.001,
            # A decaying parent needs one momentum grid for background and perturbations.
            'dncdm1.quadrature_strategy': 3,
            'dncdm1.momenta_bins': 8,
        }
        self.scenario = dict(scenario)
        self.name = self._testMethodName
        self.cosmo.set(dict(self.verbose, **scenario))
        with self.assertRaises(CosmoComputationError) as ctx:
            self.cosmo.compute()
        self.assertIn(
            "Cannot evolve tensor modes with decaying NCDM species",
            str(ctx.exception))

    def test_theta_s_shooting_matches_reference(self):
        # The most common shoot: 100*theta_s varies h (module-level target in DoShooting).
        # Candidate (lazy DoShooting) must match the reference (ctor-shooting) within tol.
        scenario = {
            '100*theta_s': 1.041783,
            'output': 'tCl',
            'l_max_scalars': 200,
        }
        candidate, reference = self._compute_candidate_and_reference(scenario)
        try:
            status = self.compare_output(
                reference,
                "Reference",
                candidate,
                "Candidate",
                COMPARE_CL_RELATIVE_ERROR,
                COMPARE_PK_RELATIVE_ERROR)
            self.assertTrue(status, "Reference comparison failed for 100*theta_s shooting")
        finally:
            reference.struct_cleanup()
            reference.empty()
            candidate.struct_cleanup()
            candidate.empty()

    def test_idm_dr_idr_perturbations_match_reference(self):
        scenario = {
            'Omega_idm_dr': 0.12,
            'xi_idr': 0.3,
            'a_idm_dr': 1e-4,
            'output': 'tCl',
            'l_max_scalars': 100,
            'k_output_values': '0.1',
        }
        candidate, reference = self._compute_candidate_and_reference(scenario)
        try:
            candidate_scalar = candidate.get_perturbations()['scalar'][0]
            reference_scalar = reference.get_perturbations()['scalar'][0]
            np.testing.assert_allclose(
                candidate_scalar['delta_idr'][:10],
                reference_scalar['delta_idr'][:10],
                rtol=self.REFERENCE_RTOL,
                atol=1e-12)
        finally:
            reference.struct_cleanup()
            reference.empty()
            candidate.struct_cleanup()
            candidate.empty()

    def test_idr_without_idm_dr_computes(self):
        scenario = {
            'N_idr': 0.34,
            'Omega_idm_dr': 0.0,
            'output': 'tCl',
            'l_max_scalars': 100,
        }
        self.scenario = dict(scenario)
        self.cosmo.set(dict(self.verbose, **scenario))
        self.cosmo.compute()
        self.assertTrue(self.cosmo.state)

    def test_tensor_spectrum_does_not_depend_on_l_max_ncdm(self):
        """With tensor_method = massless_approximation the massive neutrinos are
        evolved through the pv-owned relativistic hierarchy, which runs to
        l_max_ur; the ncdm hierarchy depth is not used at all, so the tensor
        spectrum must be independent of l_max_ncdm.

        It was not. perturb_workspace_init sized the free-streaming array s_l
        from the species present -- l_max_ncdm + 1 entries when N_ur = 0 -- while
        the hierarchy indexed it to l_max_ur, so the spectrum was built partly
        from whatever lay past the end of the buffer. l_max_ncdm therefore chose
        how much heap garbage got read: NaN about half the time, different
        numbers otherwise."""
        scenario = {
            'N_ur': 0.0,
            'N_ncdm': 1,
            'm_ncdm': 0.06,
            'deg_ncdm': 3.0,
            'output': 'tCl',
            'modes': 't',
            'gauge': 'newtonian',
        }
        self.scenario = dict(scenario)

        shallow, deep = Class(), Class()
        try:
            shallow.set(dict(self.verbose, **scenario, **{'l_max_ncdm': 10}))
            shallow.compute()
            deep.set(dict(self.verbose, **scenario, **{'l_max_ncdm': 16}))
            deep.compute()

            cl_shallow = shallow.raw_cl(100)['tt']
            cl_deep = deep.raw_cl(100)['tt']
            self.assertTrue(np.all(np.isfinite(cl_shallow)), "tensor TT is not finite")
            self.assertTrue(np.all(np.isfinite(cl_deep)), "tensor TT is not finite")
            np.testing.assert_allclose(
                cl_shallow, cl_deep, rtol=1e-10,
                err_msg="tensor TT depends on l_max_ncdm, which it does not use")
        finally:
            deep.struct_cleanup()
            deep.empty()
            shallow.struct_cleanup()
            shallow.empty()

    def test_drmd_without_idr_drmd_computes(self):
        scenario = {
            'z_stop': 1.0e4,
            'G_over_aH_drmd_ini': 1.0,
            'f_idm_drmd': 0.1,
            'delta_Neff_drmd': 0.0,
            'output': 'tCl',
            'l_max_scalars': 100,
        }
        self.scenario = dict(scenario)
        self.cosmo.set(dict(self.verbose, **scenario))
        self.cosmo.compute()
        self.assertTrue(self.cosmo.state)

    # ---- recombination = hyrec -------------------------------------------
    # Until HYREC-2 landed, nothing in the test suite set `recombination` at all
    # (issue #396), which is how a HyRec that was 1% off its own RECFAST went
    # unnoticed. These pin the combinations that were broken or rejected before.

    def test_hyrec_with_massive_neutrino_computes(self):
        """HyRec used to rebuild H(z) itself and count the massive neutrino as
        both matter and radiation, a 0.29% error at z = 1100 (issues #396, #369).
        HYREC-2 is handed CLASS's background instead."""
        scenario = {
            'recombination': 'HyRec',
            'N_ur': 2.0308,
            'N_ncdm': 1,
            'm_ncdm': 0.06,
            'output': 'tCl,pCl',
            'l_max_scalars': 500,
        }
        self.scenario = dict(scenario)
        self.cosmo.set(dict(self.verbose, **scenario))
        self.cosmo.compute()
        self.assertTrue(self.cosmo.state)

    def test_hyrec_with_scalar_field_computes(self):
        """Scalar-field dark energy is not a Fluid, so HyRec's CPL reconstruction
        silently ignored it. HYREC-2 reads the true background, so it cannot."""
        scenario = {
            'recombination': 'HyRec',
            'Omega_fld': 0,
            'Omega_scf': 0.1,
            'attractor_ic_scf': 'yes',
            'scf_parameters': '10.0, 0.0, 0.0, 0.0, 100.0, 0.0',
            'output': 'tCl',
            'l_max_scalars': 500,
        }
        self.scenario = dict(scenario)
        self.cosmo.set(dict(self.verbose, **scenario))
        self.cosmo.compute()
        self.assertTrue(self.cosmo.state)

    def test_baryon_clumping_recombines_earlier(self):
        """Three-zone baryon clumping (H0 Olympics 2107.10291 sec. 2.4.1, its
        model M1): a clumpy plasma recombines earlier, while the fully ionized
        one, and its Thomson rate, is that of the mean density. b = 0 is an
        ordinary point, shape keys and all."""
        def thermodynamics(**clumping):
            cosmo = Class(clumping)
            try:
                cosmo.compute(level=['thermodynamics'])
                return (cosmo.get_current_derived_parameters(['z_rec'])['z_rec'],
                        cosmo.get_thermodynamics())
            finally:
                cosmo.struct_cleanup()
                cosmo.empty()

        z_rec0, thermo0 = thermodynamics()
        self.assertEqual(thermodynamics(baryon_clumping_b=0.,
                                        baryon_clumping_Delta_1=0.2)[0], z_rec0)
        z_rec, thermo = thermodynamics(baryon_clumping_b=0.5)
        self.assertAlmostEqual(z_rec - z_rec0, 16.6, delta=0.5)
        for key in ('x_e', "kappa' [Mpc^-1]"):
            self.assertAlmostEqual(np.interp(1e4, thermo['z'], thermo[key])
                                   / np.interp(1e4, thermo0['z'], thermo0[key]), 1.,
                                   places=10, msg=key)

    def test_hyrec_and_recfast_agree_to_a_few_tenths_of_a_percent(self):
        """The point of the upgrade. Before it, HyRec sat ~1% from CLASS's own
        RECFAST in TT and cost dchi2 ~ +25 against Planck; the two codes should
        differ by a few tenths of a percent at most."""
        scenario = {
            'recombination': 'HyRec',
            'output': 'tCl',
            'l_max_scalars': 1000,
        }
        reference = dict(scenario, **{'recombination': 'RECFAST'})

        candidate, ref = Class(), Class()
        try:
            candidate.set(dict(self.verbose, **scenario))
            candidate.compute()
            ref.set(dict(self.verbose, **reference))
            ref.compute()

            cl_hyrec = candidate.raw_cl(1000)['tt'][2:]
            cl_recfast = ref.raw_cl(1000)['tt'][2:]
            rel = np.abs(cl_hyrec/cl_recfast - 1.0)
            self.assertLess(
                np.max(rel), 0.005,
                "HyRec and RECFAST TT differ by more than 0.5%%: max %.3e" % np.max(rel))
        finally:
            ref.struct_cleanup()
            ref.empty()
            candidate.struct_cleanup()
            candidate.empty()

    def test_hyrec_calibration_brings_helium_to_cosmorec(self):
        """hyrec_calibration = cosmorec (the default) scales HYREC-2's He I rate by a
        fit to CosmoRec (docs/superpowers/specs/2026-10-08-hyrec-helium-cosmorec-design.md).
        RECFAST's own CosmoRec calibration is an independent route to the same
        reference, so through He I recombination the calibrated HyRec must sit far
        closer to it than HYREC-2 as distributed does, and hydrogen must not move."""
        def x_e(**extra):
            cosmo = Class(extra)
            try:
                cosmo.compute(level=['thermodynamics'])
                thermo = cosmo.get_thermodynamics()
                return thermo['z'], thermo['x_e']
            finally:
                cosmo.struct_cleanup()
                cosmo.empty()

        plain = x_e(recombination='HyRec', hyrec_calibration='none')
        calibrated = x_e(recombination='HyRec')
        recfast = x_e(recombination='RECFAST')

        def gap(a, b, z):
            return np.max(np.abs(np.interp(z, *a) / np.interp(z, *b) - 1.))

        helium = np.arange(1601., 3000.)
        self.assertGreater(gap(plain, recfast, helium), 3e-3)
        self.assertLess(gap(calibrated, recfast, helium), 1e-3)
        self.assertLess(gap(calibrated, plain, np.arange(801., 1500.)), 1e-5)

    @parameterized.expand([
        ('omega_b_0.02425', {'omega_b': 0.02425, 'omega_cdm': 0.11122, 'h': 0.67763, 'YHe': 0.27884}),
        ('omega_b_0.02012', {'omega_b': 0.02012, 'omega_cdm': 0.1388, 'h': 0.68257, 'YHe': 0.22695}),
    ])
    def test_hyrec_at_tight_thermo_tolerance(self, name, cosmology):
        """With tol_thermo_integration = 1e-9, HyRec failed in ndf15 on 36 of 100
        random LCDM cosmologies -- a singular LU, or a bogus z_rec < 500. The cause was
        the zero derivative the wrapper returned below XHEII_MIN, a trap for stiff
        steps; helium now relaxes to zero there. These two failed that way.
        Completion alone would miss the bogus-z_rec mode, so z_rec is checked against
        RECFAST's at the same cosmology (they differ by ~0.2)."""
        def z_rec(**extra):
            cosmo = Class(dict(cosmology, N_ur=3.044, tol_thermo_integration=1e-9, **extra))
            try:
                cosmo.compute(level=['thermodynamics'])
                return cosmo.get_current_derived_parameters(['z_rec'])['z_rec']
            finally:
                cosmo.struct_cleanup()
                cosmo.empty()

        reference = z_rec(recombination='RECFAST')
        for calibration in ('cosmorec', 'none'):
            self.assertAlmostEqual(z_rec(recombination='HyRec', hyrec_calibration=calibration),
                                   reference, delta=1.0, msg=calibration)

    def test_varying_constants_follow_their_power_laws(self):
        """varying_fundamental_constants (Hart & Chluba, arXiv:1705.03925; the
        H0 Olympics 'varying m_e' model) moves every atomic energy level, and so
        recombination, by alpha^2 m_e, and the Thomson rate by alpha^2/m_e^2,
        above varying_transition_redshift only. RECFAST and HyRec apply this
        through different atomic physics, so each is checked."""
        def thermodynamics(alpha, me, **extra):
            cosmo = Class(dict({'varying_fundamental_constants': 'instantaneous',
                                'varying_alpha': alpha, 'varying_me': me}, **extra))
            try:
                cosmo.compute(level=['thermodynamics'])
                return (cosmo.get_current_derived_parameters(['z_rec', 'YHe']),
                        cosmo.get_thermodynamics())
            finally:
                cosmo.struct_cleanup()
                cosmo.empty()

        for recombination in ('RECFAST', 'HyRec'):
            # A fixed Y_He, since BBN would move n_e with alpha (checked below).
            fixed = {'recombination': recombination, 'YHe': 0.25}
            derived0, thermo0 = thermodynamics(1., 1., **fixed)
            for alpha, me in ((1., 1.05), (1.02, 1.)):
                msg = f"{recombination}, alpha = {alpha}, m_e = {me}"
                derived, thermo = thermodynamics(alpha, me, **fixed)
                self.assertAlmostEqual((1. + derived['z_rec'])/(1. + derived0['z_rec']),
                                       alpha**2*me, delta=3e-3, msg=msg)
                # Reionized below the transition, fully ionized at z = 1e4.
                kappa_dot = [np.interp([2., 1e4], th['z'], th["kappa' [Mpc^-1]"])
                             for th in (thermo, thermo0)]
                np.testing.assert_allclose(kappa_dot[0]/kappa_dot[1],
                                           [1., alpha**2/me**2], rtol=1e-5, err_msg=msg)

        # Only alpha moves the helium yield of BBN, linearly as in class_public.
        self.assertAlmostEqual(thermodynamics(1.02, 1.)[0]['YHe']
                               / thermodynamics(1., 1.)[0]['YHe'], 1.02, places=10)

    def test_halofit_tail_accepts_non_analytic_primordial_table(self):
        """Halofit's extended sigma tail must not query non-analytic primordial
        tables above their computed k range."""
        external_pk = os.path.join(
            os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
            'external_Pk',
            'Pk_example.dat')
        scenario = {
            'P_k_ini type': 'external_Pk',
            'command': 'cat {}'.format(external_pk),
            'output': 'mPk',
            'P_k_max_1/Mpc': 2,
            'non linear': 'halofit',
        }
        self.scenario = dict(scenario)
        self.cosmo.set(dict(self.verbose, **scenario))
        self.cosmo.compute()
        self.assertTrue(self.cosmo.state)

    def test_external_primordial_table_reaches_the_continued_lensing_source(self):
        """The CMB lensing source is continued past the last perturbation mode up to
        k_max_extrapolated = 70/Mpc, and the spectra integrate it there. An external
        table need only cover the computed modes (the bundled example stops at 10/Mpc),
        so the primordial module continues it, linearly in ln k - ln P. A power law
        tabulated to 10/Mpc must then give the analytic C_l: before that, the run stopped
        with k = 10.4 out of the table's range."""
        k = np.geomspace(1e-6, 10., 1401)
        with tempfile.TemporaryDirectory() as tmp:
            table = os.path.join(tmp, 'power_law.dat')
            np.savetxt(table, np.c_[k, 2.1e-9 * (k / 0.05)**(0.965 - 1.)], fmt='%.17e')
            base = {'output': 'tCl,pCl,lCl', 'lensing': 'yes', 'non linear': 'halofit'}
            cls = []
            for primordial in ({'P_k_ini type': 'external_Pk', 'command': 'cat {}'.format(table)},
                               {'A_s': 2.1e-9, 'n_s': 0.965, 'k_pivot': 0.05}):
                cosmo = Class(dict(self.verbose, **base, **primordial))
                try:
                    cosmo.compute()
                    cls.append((cosmo.raw_cl(2500)['pp'][2:], cosmo.lensed_cl(2500)['bb'][2:]))
                finally:
                    cosmo.struct_cleanup()
                    cosmo.empty()
        for external, analytic, name in zip(cls[0], cls[1], ('pp', 'lensed BB')):
            np.testing.assert_allclose(external, analytic, rtol=1e-9, err_msg=name)

    def test_inflation_v_default_step_size_computes(self):
        """The default inflation_V integration must not reject a valid final step."""
        scenario = {
            'P_k_ini type': 'inflation_V',
            'output': 'tCl',
            'modes': 's, t',
        }
        self.scenario = dict(scenario)
        self.cosmo.set(dict(self.verbose, **scenario))
        self.cosmo.compute()
        self.assertTrue(self.cosmo.state)

    def test_axion_monodromy_reproduces_its_analytic_delta_ns(self):
        """'potential = monodromy' is the axion-monodromy potential of
        arXiv:0907.2916, V = mu^3 (x + b f cos(x/f)) with x = V_4 - phi.

        Two things are pinned here. First the plumbing: before this shape
        existed the 'potential' key was read and thrown away for
        'P_k_ini type = inflation_V' ("only polynomial coded so far"), so
        anything but polynomial was silently ignored and V_1 was reinterpreted
        as the linear Taylor coefficient. Second the physics: the modulation
        must show up in P_s(k) with the fractional amplitude the paper derives
        in its eq. (2.9),

            delta n_s = 12 b / sqrt(1+(3 f phi_*)^2)
                        * sqrt(pi/8 * coth(pi/(2 f phi_*)) * f phi_*) ,

        with f and phi_* in REDUCED Planck units, i.e. sqrt(8 pi) times the
        V_1 and V_4 this module is given. See
        docs/superpowers/specs/2026-09-06-axion-monodromy-inflation-design.md.
        """
        import math

        sq8pi = math.sqrt(8. * math.pi)
        phi_star, f_red, b = 11., 0.02, 0.08
        x_star = phi_star/sq8pi
        # normalise to A_s ~ 2.1e-9 through the first-order slow-roll relation
        V_0 = 2.1e-9*3./(128.*math.pi*x_star**3)

        common = {
            'P_k_ini type': 'inflation_V',
            'potential': 'monodromy',
            'V_0': V_0, 'V_1': f_red/sq8pi, 'V_3': 1., 'V_4': x_star,
            'output': 'tCl',
            'modes': 's, t',
            'k_per_decade_primordial': 200,
            'primordial_inflation_tol_integration': 1e-5,
        }

        spectra = {}
        for label, b_value in (('smooth', 0.), ('modulated', b)):
            cosmo = Class()
            cosmo.set(dict(self.verbose, **common, **{'V_2': b_value}))
            cosmo.compute(level=['primordial'])
            primordial = cosmo.get_primordial()
            spectra[label] = (primordial['k [1/Mpc]'].copy(),
                              primordial['P_scalar(k)'].copy())
            cosmo.struct_cleanup()
            cosmo.empty()

        k, P_modulated = spectra['modulated']
        _, P_smooth = spectra['smooth']
        ratio = P_modulated/P_smooth - 1.

        # Project onto the paper's template cos(phi_k/f), phi_k =
        # sqrt(phi_*^2 - 2 ln(k/k_*)), with the phase free (the cos and sin
        # columns) and a quadratic in ln k to absorb any residual trend. phi_*
        # is scanned rather than fixed: k <-> phi_k is a slow-roll relation and
        # a 1% error in it is many radians of phase when f is small.
        ln_k = np.log(k/0.05)
        window = np.abs(ln_k) < 2.
        ln_k, ratio = ln_k[window], ratio[window]
        best_amplitude = 0.
        for phi_star_fit in np.linspace(10.5, 11.5, 501):
            theta = np.sqrt(phi_star_fit**2 - 2.*ln_k)/f_red
            design = np.column_stack([np.cos(theta), np.sin(theta),
                                      np.ones_like(ln_k), ln_k, ln_k**2])
            coefficients, *_ = np.linalg.lstsq(design, ratio, rcond=None)
            best_amplitude = max(best_amplitude,
                                 math.hypot(coefficients[0], coefficients[1]))

        y = f_red*phi_star
        delta_ns = (12.*b/math.sqrt(1. + (3.*y)**2)
                    * math.sqrt(math.pi/8./math.tanh(math.pi/(2.*y))*y))
        self.assertAlmostEqual(best_amplitude/delta_ns, 1., delta=0.06)

        # A wrong template frequency must NOT score: without this the test would
        # pass on any spectrum with enough structure in it.
        theta = np.sqrt(phi_star**2 - 2.*ln_k)/(f_red/3.)
        design = np.column_stack([np.cos(theta), np.sin(theta),
                                  np.ones_like(ln_k), ln_k, ln_k**2])
        coefficients, *_ = np.linalg.lstsq(design, ratio, rcond=None)
        self.assertLess(math.hypot(coefficients[0], coefficients[1]),
                        0.2*delta_ns)

    def test_pheno_axion_reports_its_peak_ede_fraction(self):
        """f_EDE at its peak is the headline derived parameter of a published
        EDE analysis (Poulin et al. 1811.04083), and the pheno-axion fluid's
        internally derived scales (m_fld, alpha_fld, omega_axion) were likewise
        unreachable from python. Both now come through the generic
        "<species key>.<field>" derived-parameter form. Resolves #367."""
        scenario = {
            'output': 'tCl',
            'fluid_equation_of_state': 'pheno_axion',
            'n_pheno_axion': 3,
            'log10_axion_ac': -3.5,
            'fraction_fld_ac': 0.1,
            'Theta_initial_fld': 2.8,
        }
        self.scenario = dict(scenario)
        self.cosmo.set(dict(self.verbose, **scenario))
        self.cosmo.compute(level=['background'])

        names = ['Fluid.f_peak', 'Fluid.a_peak', 'Fluid.z_peak', 'Fluid.a_c',
                 'Fluid.n_axion', 'Fluid.m_fld', 'Fluid.alpha_fld',
                 'Fluid.omega_axion']
        derived = self.cosmo.get_current_derived_parameters(names)

        # f_peak/a_peak/z_peak are read off the same table the background
        # columns are printed from, so they must agree with it row for row.
        background = self.cosmo.get_background()
        f_ede = background['(.)rho_fld']/background['(.)rho_crit']
        i_peak = np.argmax(f_ede)
        self.assertAlmostEqual(derived['Fluid.f_peak']/f_ede[i_peak], 1.,
                               places=12)
        self.assertAlmostEqual(derived['Fluid.z_peak']/background['z'][i_peak],
                               1., places=12)
        self.assertAlmostEqual(derived['Fluid.a_peak']
                               * (1. + background['z'][i_peak]), 1., places=12)

        # The injection peaks just after the transition, a little above the
        # fraction the input fixes AT a_c (where w has only turned halfway).
        self.assertGreater(derived['Fluid.f_peak'], 0.1)
        self.assertLess(derived['Fluid.f_peak'], 0.15)
        self.assertGreater(derived['Fluid.a_peak'], 10**-3.5)
        self.assertLess(derived['Fluid.a_peak'], 10**-3.)

        self.assertAlmostEqual(derived['Fluid.a_c']/10**-3.5, 1., places=12)
        self.assertEqual(derived['Fluid.n_axion'], 3.)
        for name in ('Fluid.m_fld', 'Fluid.alpha_fld', 'Fluid.omega_axion'):
            self.assertGreater(derived[name], 0., msg=name)

        # GetSpeciesParam answers an unknown name with 0.0, so the wrapper must
        # keep refusing names it does not know rather than reporting that 0:
        # an unlisted field, a species that is not there, and -- the case the
        # whitelist alone cannot catch -- a listed field this species does not
        # answer. Raised in Copilot's review of PR #428.
        for bad in ('Fluid.f_pea', 'NoSuchSpecies.f_peak', 'CDM.m_fld',
                    'Fluid.M_2'):
            with self.assertRaises(CosmoSevereError, msg=bad):
                self.cosmo.get_current_derived_parameters([bad])

    def test_pede_follows_its_closed_form(self):
        """(G)PEDE (H0 Olympics 2107.10291 sec. 2.5.2-3) is defined by
        rho_de(z) = rho_de,0 [1 - tanh(Delta log10(1+z))]. The fluid integrates
        rho from w(a), starting from the closed-form integral at a_ini, so the
        table matching the definition checks both. At Delta = 2.5,
        1 - tanh(x) is 1e-30 at a_ini, where the naive form cancels to zero
        and the fluid would never switch on. Delta = 0 is LambdaCDM."""
        for delta, eos, use_ppf in ((0., 'PEDE', 'yes'), (1., 'PEDE', 'no'),
                                    (2.5, 'GPEDE', 'yes')):
            cosmo = Class({'Omega_Lambda': 0., 'fluid_equation_of_state': eos,
                           'Delta_pede': delta, 'use_ppf': use_ppf})
            try:
                cosmo.compute(level=['background'])
                background = cosmo.get_background()
            finally:
                cosmo.struct_cleanup()
                cosmo.empty()
            x = delta*np.log10(1. + background['z'])
            rho = background['(.)rho_fld']
            np.testing.assert_allclose(
                np.log(rho/rho[-1]), np.log(2.) - 2.*x - np.log1p(np.exp(-2.*x)),
                atol=1e-3, err_msg=f"Delta_pede = {delta}")
            self.assertAlmostEqual(background['(.)w_fld'][-1],
                                   -1. - delta/(3.*np.log(10.)), places=10)

    def test_axion_scalar_field_peak_fraction_is_its_shooting_target(self):
        """The peak energy fraction is derived from the background table alone,
        so it is not a fluid feature: the exact Klein-Gordon axion answers it
        too, which is what makes the two axion treatments comparable. A field
        still frozen today peaks today, at the Omega_scf the shooting hit."""
        scenario = {
            'output': 'tCl',
            'Omega_scf': 0.05,
            'scf_potential': 'axion',
            'f_axion': 0.5,
            'n_axion': 1,
            'Theta_initial_scf': 2.0,
        }
        self.scenario = dict(scenario)
        self.cosmo.set(dict(self.verbose, **scenario))
        self.cosmo.compute(level=['background'])

        derived = self.cosmo.get_current_derived_parameters(
            ['ScalarField.f_peak', 'ScalarField.z_peak'])
        self.assertAlmostEqual(derived['ScalarField.f_peak'], 0.05, places=4)
        self.assertLess(derived['ScalarField.z_peak'], 1e-3)

    def test_composite_species_reports_its_peak_fraction(self):
        """A composite registers background slots for its children and none for
        itself, so its density exists only through Rho(). Reading the table by
        the species' own rho index would silently report 0 for exactly the
        species whose total has to be asked for. Raised in Copilot's review of
        PR #428."""
        scenario = {
            'output': 'tCl',
            'Omega_dcdmdr': 0.12,
            'Gamma_dcdm': 10.0,
        }
        self.scenario = dict(scenario)
        self.cosmo.set(dict(self.verbose, **scenario))
        self.cosmo.compute(level=['background'])

        derived = self.cosmo.get_current_derived_parameters(
            ['DCDM_DR.f_peak', 'DCDM_DR.a_peak', 'DCDM_DR.z_peak'])
        background = self.cosmo.get_background()
        # The composite total is its children summed, which is what Rho() does.
        f_total = ((background['(.)rho_dcdm'] + background['(.)rho_dr'])
                   / background['(.)rho_crit'])
        i_peak = np.argmax(f_total)
        self.assertGreater(derived['DCDM_DR.f_peak'], 0.)
        self.assertAlmostEqual(derived['DCDM_DR.f_peak']/f_total[i_peak], 1.,
                               places=12)
        self.assertAlmostEqual(derived['DCDM_DR.z_peak']
                               / background['z'][i_peak], 1., places=12)
        self.assertAlmostEqual(derived['DCDM_DR.a_peak']
                               * (1. + background['z'][i_peak]), 1., places=12)

    def _nede_trigger(self, params):
        """Cold NEDE with f_NEDE = 0.1, H_* = 0.2 m_NEDE on top of params:
        rho_NEDE/(f_NEDE H_*^2) just before the transition, and the rest of the
        universe at z_star over H_*^2 (1 - f_NEDE). Both are 1 where the
        transition is placed right."""
        f_NEDE, m_NEDE = 0.1, 10**2.5
        H_star = 0.2*m_NEDE
        cosmo = Class(dict({'f_NEDE': f_NEDE, 'm_NEDE': m_NEDE}, **params))
        try:
            cosmo.compute(level=['background'])
            z_star = cosmo.get_current_derived_parameters(['NEDE.z_star'])['NEDE.z_star']
            background = cosmo.get_background()
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()
        z = background['z']
        nede = background['(.)rho_NEDE']
        others = background['(.)rho_tot'] - nede
        # The two table rows just before the transition, where `others` is smooth.
        before = np.where(z > z_star)[0]
        i, j = before[np.argsort(z[before])[:2]]
        slope = np.log(others[j]/others[i])/np.log((1. + z[j])/(1. + z[i]))
        return (nede[i]/(f_NEDE*H_star**2),
                others[i]*((1. + z_star)/(1. + z[i]))**slope/(H_star**2*(1. - f_NEDE)))

    def test_nede_transition_is_where_H_reaches_its_trigger(self):
        """Cold NEDE (H0 Olympics 2107.10291 sec. 2.4.5) decays where H falls to
        H_* = H_over_m_NEDE m_NEDE, and f_NEDE is its share of the density
        there: rho_NEDE = f_NEDE H_*^2 before the transition, and the rest of
        the universe reaches H_*^2 (1 - f_NEDE) at z_star. The massive neutrino
        is still relativistic at z_star; counting it as matter would move the
        transition."""
        nede, others = self._nede_trigger({'N_ur': 2.0328, 'N_ncdm': 1, 'm_ncdm': 0.06})
        self.assertAlmostEqual(nede, 1., places=12)
        self.assertAlmostEqual(others, 1., delta=1e-4)

    def test_nede_transition_counts_the_majoron_before_it_heats_the_neutrinos(self):
        """NEDE places its transition before the background is solved, so the
        majoron gives its density from the fit. Today's density scaled back as
        a^-4 counts the Delta N_eff ~ 0.11 of the majoron decays from the start,
        and put the rest of the universe 0.8% below H_*^2 (1 - f_NEDE) at z_star
        (T_gamma ~ m_phi = 1 eV, before most of the heating)."""
        _, others = self._nede_trigger({'m_majoron': 1., 'Gamma_eff_majoron': 10., 'N_ur': 0})
        self.assertAlmostEqual(others, 1., delta=1e-4)

    def test_nede_spectra_match_triggerclass(self):
        """NEDE/LambdaCDM in lensed TT and EE against TriggerCLASS (commit
        174227ad, z_decay_NEDE = 4589.54 -- this z_star -- with three_eos_NEDE =
        three_ceff2_NEDE = 2, no viscosity and a negligible trigger), for the
        same cosmology. The two agree to 2.1e-3 over 2 <= l <= 2500. This
        exercises the trigger, the junction conditions, and carrying NEDE
        across the later tca/rsa switches: dropping it there costs 1e-2."""
        cosmology = {'h': 0.6781, 'omega_b': 0.0223828, 'omega_cdm': 0.1201075,
                     'N_ur': 3.046, 'YHe': 0.245, 'A_s': 2.1e-9, 'n_s': 0.9660499,
                     'z_reio': 11.357, 'output': 'tCl,pCl,lCl', 'lensing': 'yes',
                     'l_max_scalars': 2500}
        spectra = []
        for nede in ({}, {'f_NEDE': 0.1, 'm_NEDE': 10**2.5, 'w_NEDE': 2./3.}):
            cosmo = Class(dict(cosmology, **nede))
            try:
                cosmo.compute()
                spectra.append(cosmo.lensed_cl(1500))
            finally:
                cosmo.struct_cleanup()
                cosmo.empty()
        triggerclass = {'tt': [0.99554, 1.03684, 1.06018],
                        'ee': [1.16450, 1.00847, 0.88217]}
        ell = [500, 1000, 1500]
        for key, expected in triggerclass.items():
            ratio = spectra[1][key][ell]/spectra[0][key][ell]
            np.testing.assert_allclose(ratio, expected, atol=3e-3, err_msg=key)

    def test_scalar_tensor_closes_on_H0(self):
        """The scalar-tensor species shoots its potential's constant so that its
        effective density vanishes today, i.e. H(z=0) = H0, beside Lambda's usual
        closure. The massless strongly coupled case (xi = 0.3, lambda = 0) grows
        f to ~70 by today, far from the guess: the unknown once lived in 1/Mpc^2,
        where Newton's step tolerance of 1e-3 declared convergence after one step
        and left H(0) 34% off. 100*theta_s shoots h at the same time."""
        emg = {'xi_st': 0.3, 'phi_ini_st': 0.49, 'lambda_st': 3.42e-108}
        massless = {'xi_st': 0.3, 'phi_ini_st': 0.49}
        for model, extra in ((emg, {'h': 0.6781}), (massless, {'h': 0.6781}),
                             (emg, {'100*theta_s': 1.04092})):
            cosmo = Class(dict(model, **extra))
            try:
                cosmo.compute(level=['thermodynamics'])
                H_today = cosmo.get_background()['H [1/Mpc]'][-1]
                h = cosmo.h()
            finally:
                cosmo.struct_cleanup()
                cosmo.empty()
            self.assertAlmostEqual(H_today * 299792.458 / (100. * h), 1., delta=1e-3,
                                   msg=str(model))

    def test_scalar_tensor_refuses_what_it_does_not_solve(self):
        """The species replaces GR's synchronous-gauge scalar equations only. A
        Newtonian gauge, tensor modes, curvature, a PPF fluid (another closure of
        the Einstein equations), NEDE (which locates its transition from the
        densities before the background exists), or N-body-gauge sources would
        silently run GR's or unchecked equations; each is refused instead."""
        from classy import CosmoSevereError
        emg = {'xi_st': 0.3, 'phi_ini_st': 0.49, 'lambda_st': 3.42e-108}
        refused = [{'gauge': 'newtonian'},
                   {'modes': 's,t', 'r': 0.1},
                   {'Omega_k': 0.01},
                   {'Omega_Lambda': 0., 'w0_fld': -0.9, 'use_ppf': 'yes'},
                   {'f_NEDE': 0.1, 'm_NEDE': 10**2.5},
                   # N-body-gauge sources use GR's Poisson equation and p_tot'; both
                   # routes that switch them on are refused.
                   {'output': 'mTk', 'Nbody gauge transfer functions': 'yes'},
                   {'output': 'dTk', 'extra metric transfer functions': 'yes'}]
        for extra in refused:
            with self.assertRaises(CosmoSevereError, msg=str(extra)):
                cosmo = Class(dict(dict(emg, output='tCl'), **extra))
                try:
                    cosmo.compute()
                finally:
                    cosmo.struct_cleanup()
                    cosmo.empty()

    def test_scalar_tensor_solar_system_parameters(self):
        """G_eff_today (Cavendish, over the action's G) and gamma_PPN - 1 today, for
        a massless field (Boisseau et al. gr-qc/0001066). EMG's field rolls to 0,
        so both are GR's; induced gravity keeps a coupling today and must match
        the formulas at the last background row. A massive neutrino rides along:
        it is minimally coupled and nothing about the species should notice it."""
        emg = {'xi_st': 0.3, 'phi_ini_st': 0.49, 'lambda_st': 3.42e-108,
               'N_ur': 2.0328, 'N_ncdm': 1, 'm_ncdm': 0.06}
        xi = 1./200.
        induced = {'N_st': 0., 'xi_st': xi, 'phi_ini_st': (1./xi)**0.5}
        names = ['ScalarTensor.G_eff_today', 'ScalarTensor.gamma_PPN_minus_1']
        for model in (emg, induced):
            cosmo = Class(dict(model, output='tCl'))
            try:
                cosmo.compute()
                derived = cosmo.get_current_derived_parameters(names)
                phi = cosmo.get_background()['phi_st'][-1]
            finally:
                cosmo.struct_cleanup()
                cosmo.empty()
            f = model.get('N_st', 1.) + model['xi_st']*phi**2
            f_phi = 2.*model['xi_st']*phi
            self.assertAlmostEqual(derived[names[0]],
                                   (2.*f + 4.*f_phi**2)/(f*(2.*f + 3.*f_phi**2)), places=12)
            self.assertAlmostEqual(derived[names[1]], -f_phi**2/(f + 2.*f_phi**2), places=12)
            if model is emg:
                self.assertAlmostEqual(derived[names[0]], 1., places=6)
                self.assertAlmostEqual(derived[names[1]], 0., places=6)

    def test_scalar_tensor_matches_hi_class_brans_dicke(self):
        """Brans-Dicke against hi_class (public hi_class, commit ece81a57,
        gravity_model = brans_dicke, omega_BD = 50, phi_ini = 1, closure through
        Omega_smg = -1, run at tol_perturbations_integration = 1e-6 and
        perturbations_sampling_stepsize = 0.02). The scalar-tensor species is
        induced gravity there: N = 0, xi = 1/(4 omega_BD), phi_BD = xi phi^2.
        hi_class solves it with alpha functions and a V_X variable, written
        independently of this code; the two agree to 2e-7 in f = M*^2, 3.5e-5 in
        H(z), 1e-4 in P(k) and 7e-4 in the lensed C_l pinned here, on effects of
        20-110%. The C_l are compared at hi_class's source sampling: at the
        default perturb_sampling_stepsize = 0.1 the ratio is 6e-3 off at l < 200
        and 3e-3 at l = 220 (7e-3 with no transfer cut), which says nothing
        about the scalar-tensor species.

        The C_l pins are hi_class's converged in l (l_linstep 5; at its default
        40 its cubic spline moves them by up to 7e-4). Converged in l, this
        code is within 2e-4 of them except at EE l = 810, where the ratio is
        steepest: 1.1e-3, a difference between the codes (the scalar-tensor
        spec measured EE to 1.1e-3 rms). Hence the C_l tolerance 1.5e-3.
        hi_class also ran at M2_tuning_smg = no, pert_initial_conditions_smg =
        zero, with the cosmology below."""
        cosmology = {'h': 0.6781, 'omega_b': 0.02238, 'omega_cdm': 0.1201,
                     'N_ur': 3.044, 'YHe': 0.245, 'tau_reio': 0.0544,
                     'A_s': 2.1e-9, 'n_s': 0.9660, 'output': 'tCl,pCl,lCl,mPk',
                     'lensing': 'yes', 'l_max_scalars': 2000, 'P_k_max_1/Mpc': 1.,
                     'perturb_sampling_stepsize': 0.02}
        brans_dicke = {'N_st': 0., 'xi_st': 1./(4.*50.), 'phi_ini_st': (1./(1./200.))**0.5}
        h = cosmology['h']
        results = []
        for model in ({}, brans_dicke):
            cosmo = Class(dict(cosmology, **model))
            try:
                cosmo.compute()
                results.append({'cl': cosmo.lensed_cl(2000),
                                'pk': [cosmo.pk_lin(k*h, 0.) for k in (0.01, 0.03, 0.2, 0.5)],
                                'background': cosmo.get_background()})
            finally:
                cosmo.struct_cleanup()
                cosmo.empty()
        hi_class = {'tt': ([220, 810, 1000, 1500, 2000], [0.99949, 0.95516, 1.16681, 1.15800, 1.02913]),
                    'ee': ([100, 220, 810, 1000], [0.93748, 0.84849, 1.42870, 0.94745])}
        for key, (ell, expected) in hi_class.items():
            ratio = results[1]['cl'][key][ell]/results[0]['cl'][key][ell]
            np.testing.assert_allclose(ratio, expected, atol=1.5e-3, err_msg=key)
        np.testing.assert_allclose(np.array(results[1]['pk'])/np.array(results[0]['pk']),
                                   [1.21923, 1.21351, 1.22084, 1.22645], atol=1e-3, err_msg='P(k)')
        background = results[1]['background']
        z = background['z'][::-1]
        f = np.interp(np.log1p([1100., 100., 0.]), np.log1p(z), background['f_st'][::-1])
        np.testing.assert_allclose(f, [1.022987, 1.064279, 1.171768], rtol=1e-5, err_msg='f')

    def test_monodromy_parameter_domain_is_rejected(self):
        """The monodromy potential divides by the decay constant V_1, so V_1 = 0
        reaches primordial_inflation_check_potential as a NaN -- and its V<=0 and
        dV>=0 tests are ordered comparisons, which do not catch one. The other
        three bound the branch the field is supposed to roll down. Raised in
        review of PR #415."""
        import math

        sq8pi = math.sqrt(8. * math.pi)
        good = {'P_k_ini type': 'inflation_V', 'potential': 'monodromy',
                'V_0': 1.4830727019323157e-12, 'V_1': 0.02/sq8pi, 'V_2': 0.08,
                'V_3': 1., 'V_4': 11./sq8pi,
                'output': 'tCl', 'modes': 's, t'}
        for key, bad_value in (('V_1', 0.),      # decay constant -> division by zero
                               ('V_1', -1e-3),   # same model as (b, f) -> (-b, -f)
                               ('V_3', 0.),      # monodromy power
                               ('V_0', 0.),      # monodromy scale
                               ('V_4', 0.)):     # no branch left to roll down
            scenario = dict(self.verbose, **good)
            scenario[key] = bad_value
            cosmo = Class()
            cosmo.set(scenario)
            with self.assertRaises(CosmoComputationError,
                                   msg=f'{key} = {bad_value} was accepted'):
                cosmo.compute(level=['primordial'])
            cosmo.struct_cleanup()
            cosmo.empty()

    def test_unknown_inflation_potential_is_rejected(self):
        """The 'potential' key used to be parsed and discarded, so a typo in it
        ran the polynomial shape with whatever V_i happened to be set."""
        self.cosmo.set(dict(self.verbose, **{
            'P_k_ini type': 'inflation_V',
            'potential': 'no_such_potential',
            'output': 'tCl',
            'modes': 's, t',
        }))
        self.assertRaises(CosmoSevereError, self.cosmo.compute)

    def test_z_max_pk_above_the_thermodynamics_table_computes(self):
        """thermodynamics_at_z extrapolates analytically above the tabulated
        range, and that branch used to be the one path through the function that
        never wrote its *last_index out-parameter (#380). The no-CMB source
        sampler seeds its interpolation cursor through exactly that call, so the
        cursor stayed indeterminate and the first inter_closeby handed garbage to
        array_hunt_growing_closeby.

        recfast_z_initial = z_max_pk + 1 puts the seed lookup at the table edge,
        which is what walks it onto the unguarded path -- it is what
        notebooks/many_times.ipynb does."""
        self.cosmo.set(dict(self.verbose, **{
            'output': 'mTk',
            'z_max_pk': 46000,
            'recfast_Nz0': 92000,
            'recfast_z_initial': 46001.,
            'k_per_decade_for_pk': 400,
            'k_per_decade_for_bao': 400,
            'k_min_tau0': 40.,
            'perturb_sampling_stepsize': 0.05,
            'P_k_max_1/Mpc': 1.0,
            'compute damping scale': 'yes',
            'gauge': 'newtonian',
        }))
        self.cosmo.compute()
        self.assertIn('mTk', self.cosmo.pars['output'])

    def test_lensed_polarisation_without_temperature(self):
        """output = pCl,lCl with lensing is a documented combination. lensing_init()
        used to read the unlensed TT spectrum through an index that is only set when
        TT is computed (AddressSanitizer: BUS at lensing_module.cpp). Lensed EE, BB
        and phiphi do not depend on TT, so they must equal those of tCl,pCl,lCl."""
        spectra = {}
        for output in ('tCl,pCl,lCl', 'pCl,lCl'):
            self.cosmo.set(dict(self.verbose, output=output, lensing='yes', l_max_scalars=1500))
            self.cosmo.compute()
            spectra[output] = self.cosmo.lensed_cl(1500)
            self.cosmo.struct_cleanup()
            self.cosmo.empty()
        reference, candidate = spectra['tCl,pCl,lCl'], spectra['pCl,lCl']
        self.assertNotIn('tt', candidate)
        for key in ('ee', 'bb', 'pp'):
            np.testing.assert_allclose(candidate[key], reference[key], rtol=1e-12, atol=0., err_msg=key)

    def test_cmb_tt_at_small_l_max(self):
        """The k range of the scalar CMB C_l's was k_max tau0 = 1.8 l_max: 360 at
        l_max_scalars = 200, and unlensed TT came out 11% low at l = 200. It is now at
        least l_max + k_max_tau0_above_l_max (2026-10-27): the same C_l's asked for up to
        l = 200 or up to l = 1000 must agree there (5e-5)."""
        spectra = {}
        for l_max in (200, 1000):
            self.cosmo.set(dict(self.verbose, output='tCl', l_max_scalars=l_max))
            self.cosmo.compute()
            spectra[l_max] = np.asarray(self.cosmo.raw_cl(200)['tt'])[2:]
            self.cosmo.struct_cleanup()
            self.cosmo.empty()
        np.testing.assert_allclose(spectra[200], spectra[1000], rtol=1e-3, atol=0.)
    def test_pk_equal_with_lensing(self):
        """halofit's Pk_equal method tabulates w0 and Omega_m from z = 0 to pk_eq_z_max, but the
        lensing potential asks for halofit at higher z, to where sigma(R) falls below 1 (beyond
        z = 7). That lookup failed in every Pk_equal run with lensing (array_search_bracket) until
        2026-10-27; beyond the table the values at pk_eq_z_max stand in (C_L^phiphi moves by 8.5e-7
        when the table reaches z = 10 instead, too slow to repeat here: 44 s)."""
        self.cosmo.set(dict(self.verbose, output='tCl,lCl', lensing='yes', l_max_scalars=1500,
                            Omega_Lambda=0., w0_fld=-0.9, wa_fld=-0.3,
                            **{'non linear': 'halofit', 'pk_eq': 'yes'}))
        self.cosmo.compute()
        pp = np.asarray(self.cosmo.lensed_cl(1500)['pp'])[2:]
        self.cosmo.struct_cleanup()
        self.cosmo.empty()
        self.assertTrue(np.all(np.isfinite(pp)) and np.all(pp > 0.))

    def test_kept_bessel_table_does_not_change_results(self):
        """A classy object keeps the transfer module's Bessel table from one cosmology to the
        next (include/bessel_table_cache.h). The table is reused only when it reaches far enough
        in x, and tables are built on a fixed spacing, so a cosmology computed after another
        must give exactly what it gives when computed first. The earlier cosmology here needs a
        larger table (reused), a smaller one (rebuilt), or is a slightly different cosmology."""
        base = dict(self.verbose, output='tCl,pCl,lCl', lensing='yes', l_max_scalars=1500)

        def lensed(cosmo, extra):
            cosmo.set(dict(base, **extra))
            cosmo.compute()
            cl = cosmo.lensed_cl(1500)
            cosmo.struct_cleanup()
            cosmo.empty()
            return cl

        fresh = lensed(Class(), {})
        for before in ({'k_max_tau0_over_l_max': 2.5}, {'k_max_tau0_over_l_max': 1.2}, {'h': 0.70}):
            cosmo = Class()
            lensed(cosmo, before)
            again = lensed(cosmo, {})
            for key in ('tt', 'ee', 'te', 'bb', 'pp'):
                np.testing.assert_array_equal(again[key], fresh[key], err_msg=f'{key} after {before}')

    def test_retired_hyrec_file_parameters_are_rejected(self):
        """HYREC-2 takes one directory, so the three per-file paths are gone. A
        user who set them should hear about it rather than be ignored."""
        self.cosmo.set(dict(self.verbose, **{
            'recombination': 'HyRec',
            'Alpha_inf hyrec file': '/nowhere/Alpha_inf.dat',
            'output': 'tCl',
            'l_max_scalars': 100,
        }))
        with self.assertRaises(CosmoSevereError) as ctx:
            self.cosmo.compute()
        self.assertIn('hyrec_path', str(ctx.exception))

    def test_dot_syntax_standard_partial_field_uses_legacy_defaults(self):
        scenario = {
            **self._dot_syntax_base(),
            'nu1.type': 'ncdm_standard',
            'nu1.m': 0.06,
            'nu1.deg': 2.0,
            'nu2.type': 'ncdm_standard',
            'nu2.m': 0.08,
        }
        reference = dict(scenario, **{'nu2.deg': 1.0})
        self._assert_scenarios_match(scenario, reference, "Explicit legacy default")

    def test_dot_syntax_standard_fluid_approximation_must_match(self):
        scenario = {
            **self._dot_syntax_base(),
            'nu1.type': 'ncdm_standard',
            'nu1.m': 0.06,
            'nu1.fluid_approximation': 2,
            'nu2.type': 'ncdm_standard',
            'nu2.m': 0.08,
            'nu2.fluid_approximation': 3,
        }
        self._assert_compute_fails(scenario, "must be identical for all dot-syntax NCDM-family species")

    def test_dot_syntax_unrecognised_species_type_is_rejected(self):
        """#430: a misspelled <instance>.type built no species. classy already refused
        the run for its unread keys; the error now names the cause and the fix."""
        scenario = {
            **self._dot_syntax_base(),
            'nuM.type': 'ncdm',
            'nuM.m': 0.06,
        }
        self._assert_compute_fails(scenario, "'ncdm' is not a species type")

    def test_dot_syntax_fluid_with_pk_eq_matches_legacy(self):
        """pk_eq rebuilds the input with an effective w0_fld/wa_fld. With the fluid in
        dot syntax, that override used to abort as "input sets both legacy key 'w0_fld'
        and dot-syntax 'f.w0'", while the legacy spelling ran."""
        base = {
            'output': 'mPk',
            'non_linear': 'halofit',
            'pk_eq': 'yes',
            'Omega_Lambda': 0,
            'P_k_max_1/Mpc': 1.0,
        }
        scenario = dict(base, **{'f.type': 'fluid', 'f.w0': -0.9, 'f.wa': 0.1})
        reference = dict(base, **{'w0_fld': -0.9, 'wa_fld': 0.1})
        self._assert_scenarios_match(scenario, reference, "Legacy fluid keys")

    def test_dot_syntax_standard_psd_filenames_follow_true_flags(self):
        scenario = {
            **self._dot_syntax_base(),
            'nu1.type': 'ncdm_standard',
            'nu1.m': 0.06,
            'nu1.use_psd_file': 1,
            'nu1.psd_filename': '../psd_FD_single.dat',
            'nu2.type': 'ncdm_standard',
            'nu2.m': 0.08,
            'nu2.use_psd_file': 0,
        }
        self._assert_compute_succeeds(scenario)

    def test_dot_syntax_interacting_partial_field_uses_legacy_defaults(self):
        scenario = {
            **self._dot_syntax_base(),
            'nu1.type': 'ncdm_self_interacting',
            'nu1.m': 0.06,
            'nu1.deg': 2.0,
            'nu2.type': 'ncdm_self_interacting',
            'nu2.m': 0.08,
        }
        reference = dict(scenario, **{'nu2.deg': 1.0})
        self._assert_scenarios_match(scenario, reference, "Explicit legacy default")

    def test_dot_syntax_interacting_single_species_both_geff_forms_rejected(self):
        # A *single* species may not set both G_eff and log10G_eff (ambiguous).
        # Distinct species are free to use different representations: each stores
        # its own G_eff locally, so there is no conflict to reject.
        scenario = {
            **self._dot_syntax_base(),
            'nu1.type': 'ncdm_self_interacting',
            'nu1.m': 0.06,
            'nu1.G_eff': 1e-4,
            'nu1.log10G_eff': -4.0,
        }
        self._assert_compute_fails(scenario, "specify exactly one of G_eff or log10G_eff")

    def test_dot_syntax_interacting_mixed_geff_representations_across_species_allowed(self):
        # nu1 uses the linear G_eff, nu2 uses log10G_eff. Species are configured
        # independently, so this must succeed -- AND nu2's log10G_eff=-4.0 must
        # resolve to exactly the same physics as G_eff=1e-4. Compare against a
        # reference where nu2 uses the linear form directly, so the test fails if
        # log10G_eff is ever parsed wrong or silently ignored.
        scenario = {
            **self._dot_syntax_base(),
            'nu1.type': 'ncdm_self_interacting',
            'nu1.m': 0.06,
            'nu1.G_eff': 1e-4,
            'nu2.type': 'ncdm_self_interacting',
            'nu2.m': 0.08,
            'nu2.log10G_eff': -4.0,
        }
        reference = dict(scenario)
        del reference['nu2.log10G_eff']
        reference['nu2.G_eff'] = 1e-4
        self._assert_scenarios_match(scenario, reference, "Linear G_eff for nu2")

    def test_scalar_field_synchronous_gauge_computes(self):
        scenario = {
            'output': 'tCl',
            'gauge': 'synchronous',
            'Omega_fld': 0,
            'Omega_scf': 0.1,
            'attractor_ic_scf': 'yes',
            'scf_parameters': '10, 0, 0, 0',
        }
        self._assert_compute_succeeds(scenario)

    def test_majoron_heats_the_neutrinos_and_damps_their_free_streaming(self):
        """Majoron (Escudero & Witte 1909.04044, 2103.03249; H0 Olympics
        2107.10291 sec. 2.3.6): on Escudero & Witte's fitted background, majoron
        decays raise N_eff by ~0.11 near T ~ m_phi, and the l >= 2 multipoles
        relax at the inverse-decay rate. For Gamma_eff -> 0 the fluid is three
        free-streaming neutrinos."""
        def run(**params):
            cosmo = Class(dict({'output': 'tCl', 'N_ur': 0}, **params))
            try:
                cosmo.compute()
                return cosmo.get_background(), cosmo.raw_cl(1000)['tt']
            finally:
                cosmo.struct_cleanup()
                cosmo.empty()

        per_neutrino = 7./8.*(4./11.)**(4./3.)
        background, coupled = run(m_majoron=1., Gamma_eff_majoron=10.)
        N_eff = background['(.)rho_nuphi']/background['(.)rho_g']/per_neutrino
        self.assertAlmostEqual(N_eff[-1] - N_eff[0], 0.1127, delta=2e-3)

        background, weak = run(m_majoron=1., Gamma_eff_majoron=1e-4)
        N_eff = background['(.)rho_nuphi'][0]/background['(.)rho_g'][0]/per_neutrino
        _, free = run(N_ur=N_eff)
        np.testing.assert_allclose(weak[2:], free[2:], rtol=1e-3)

        # The damping of free streaming, pinned at m_phi = 1 eV, Gamma_eff = 10.
        self.assertAlmostEqual(coupled[600]/weak[600], 1.145, delta=0.01)

    def test_type3_synchronous_computes(self):
        # Type-3 composite (coupled CDM + beta scalar field) end-to-end, in the
        # paper regime (arXiv:1604.04222): the composite builds an injectable
        # 1EXP V0*exp(-lambda*phi) scalar field with non-attractor FROZEN ICs and
        # shoots V0 to hit Omega_scf, plus an uncoupled CDM whose theta stays 0
        # (no coupling terms yet -> Task 5). Frozen ICs keep the field negligible
        # during radiation domination, so beta=-0.5 (the (1-2*beta)=2 kinetic
        # boost) no longer trips the radiation-start Omega_r check.
        # scf_parameters = 'V0_placeholder, lambda'; tuning index 0 shoots V0.
        scenario = {
            'output': 'tCl mPk',
            'gauge': 'synchronous',
            'Omega_fld': 0,
            'Omega_scf': 0.7,
            'attractor_ic_scf': 'no',
            'scf_parameters': '1, 1.22',
            'scf_veta': -0.5,
        }
        self._assert_compute_succeeds(scenario)

    def test_type3_newtonian_rejected(self):
        # Type3Species::SupportsNewtonianGauge() is false, so the perturbation module
        # refuses Newtonian gauge before integrating anything, for any nonzero beta.
        scenario = {
            'output': 'tCl',
            'gauge': 'newtonian',
            'Omega_fld': 0,
            'Omega_scf': 0.7,
            'attractor_ic_scf': 'no',
            'scf_parameters': '1, 1.22',
            'scf_veta': -0.5,
        }
        self.scenario = dict(scenario)
        self.cosmo.set(dict(self.verbose, **scenario))
        with self.assertRaises(Exception):
            self.cosmo.compute()

    def test_type3_beta_zero_unset_is_plain_scf(self):
        # scf_veta absent => no composite => plain CDM + scalar field via the
        # default exponential-quintessence potential + attractor ICs (unchanged).
        scenario = {
            'output': 'tCl',
            'gauge': 'synchronous',
            'Omega_fld': 0,
            'Omega_scf': 0.1,
            'attractor_ic_scf': 'yes',
            'scf_parameters': '10, 0, 0, 0',
        }
        self._assert_compute_succeeds(scenario)

    def test_type3_drives_theta_cdm(self):
        # Coupled synchronous CDM carries theta_cdm as a dynamical variable; with
        # the momentum transfer OFF (beta=0) it stays identically 0. The coupling
        # (beta=-0.5) must drive it away from zero. Uses the proven frozen-IC
        # injectable potential ('V0, lambda'); the brief's attractor + 4-param form
        # does not compute at beta=-0.5 (singular matrix), which is exactly why
        # test_type3_synchronous_computes also uses frozen ICs.
        scenario = {
            'output': 'mPk',
            'gauge': 'synchronous',
            'Omega_fld': 0,
            'Omega_scf': 0.7,
            'attractor_ic_scf': 'no',
            'scf_parameters': '1, 1.22',
            'scf_veta': -0.5,
            'k_output_values': '0.1',
        }
        self.scenario = dict(scenario)
        self.cosmo.set(dict(self.verbose, **scenario))
        self.cosmo.compute()
        scalar = self.cosmo.get_perturbations()['scalar'][0]
        assert np.max(np.abs(scalar['theta_cdm'])) > 1e-3

    def test_type3_suppresses_growth(self):
        # Paper (arXiv:1604.04222, Fig. 1/2): negative beta suppresses growth.
        # Controlled comparison: two Type-3 composites with the IDENTICAL frozen-IC
        # injectable potential, differing ONLY in beta. beta=-1e-4 is the
        # effectively-uncoupled baseline (CreateAll drops the composite at beta==0,
        # so an exact-zero-beta composite is not constructible). The brief's
        # plain-scf (no scf_veta) baseline is NOT a controlled comparison: it takes
        # a different potential/IC path. Larger |beta| => smaller P(k).
        base = {
            'output': 'mPk', 'P_k_max_1/Mpc': 1.0, 'z_pk': 0,
            'gauge': 'synchronous', 'Omega_fld': 0, 'Omega_scf': 0.7,
            'attractor_ic_scf': 'no', 'scf_parameters': '1, 1.22',
        }
        uncoupled = Class()
        uncoupled.set(dict(self.verbose, scf_veta=-1e-4, **base))
        uncoupled.compute()
        coupled = Class()
        coupled.set(dict(self.verbose, scf_veta=-0.5, **base))
        coupled.compute()
        k = 0.1
        try:
            ratio = coupled.pk(k, 0) / uncoupled.pk(k, 0)
            assert ratio < 1.0  # negative beta suppresses growth (paper Fig. 1/2)
        finally:
            for c in (uncoupled, coupled):
                c.struct_cleanup()
                c.empty()

    def test_type3_pk_ratio_turnup_feature(self):
        # Paper (arXiv:1604.04222, Fig. 2 right): as |beta| increases from 0 the
        # suppression of P(k) first deepens (monotonic at fixed k), then the knee
        # shifts to higher k so that a fixed intermediate-k point recovers — the
        # "turn-up" feature.  We assert two coarse proxies:
        #   (i)  moderate coupling (beta=-0.5) suppresses P(k) relative to the
        #        effectively-uncoupled reference (beta=-1e-4);
        #   (ii) very large |beta| (beta=-1e4) gives LESS suppression at k=0.2
        #        than the intermediate coupling (beta=-1e2), confirming the
        #        turn-up direction.
        # Reference baseline: scf_veta=-1e-4 keeps the run on the identical
        # Type-3 composite code path while making the momentum-transfer negligible
        # (CreateAll drops the composite for beta==0).
        base = {
            'output': 'mPk', 'P_k_max_1/Mpc': 1.0, 'z_pk': 0,
            'gauge': 'synchronous', 'Omega_fld': 0, 'Omega_scf': 0.7,
            'attractor_ic_scf': 'no', 'scf_parameters': '1, 1.22',
        }
        ref = Class()
        ref.set(dict(self.verbose, scf_veta=-1e-4, **base))
        ref.compute()
        k = 0.2
        pk0 = ref.pk(k, 0)
        ratios = {}
        try:
            for beta in (-0.5, -1e2, -1e4):
                c = Class()
                c.set({**self.verbose, **base, 'scf_veta': beta})
                try:
                    c.compute()
                    ratios[beta] = c.pk(k, 0) / pk0
                finally:
                    c.struct_cleanup()
                    c.empty()
        finally:
            ref.struct_cleanup()
            ref.empty()
        # (i) Suppression at moderate |beta|
        assert ratios[-0.5] < 1.0, f"Expected suppression, got ratio={ratios[-0.5]:.4f}"
        # (ii) Turn-up: very large |beta| shows less suppression than intermediate |beta|
        assert ratios[-1e4] > ratios[-1e2], (
            f"Expected turn-up: ratio[-1e4]={ratios[-1e4]:.4f} should exceed "
            f"ratio[-1e2]={ratios[-1e2]:.4f}"
        )

    def test_dot_syntax_interacting_psd_filenames_follow_true_flags(self):
        scenario = {
            **self._dot_syntax_base(),
            'nu1.type': 'ncdm_self_interacting',
            'nu1.m': 0.06,
            'nu1.use_psd_file': 1,
            'nu1.psd_filename': '../psd_FD_single.dat',
            'nu2.type': 'ncdm_self_interacting',
            'nu2.m': 0.08,
            'nu2.use_psd_file': 0,
        }
        self._assert_compute_succeeds(scenario)

    def test_dot_syntax_decay_dr_partial_field_uses_legacy_defaults(self):
        scenario = {
            **self._dot_syntax_base(),
            'nu1.type': 'ncdm_decay_dr',
            'nu1.m': 0.06,
            'nu1.Gamma': 1e-3,
            'nu1.deg': 2.0,
            'nu2.type': 'ncdm_decay_dr',
            'nu2.m': 0.08,
            'nu2.Gamma': 2e-3,
            # A decaying parent needs one momentum grid for background and perturbations.
            'nu1.quadrature_strategy': 3,
            'nu1.momenta_bins': 8,
            'nu2.quadrature_strategy': 3,
            'nu2.momenta_bins': 8,
        }
        reference = dict(scenario, **{'nu2.deg': 1.0})
        self._assert_scenarios_match(scenario, reference, "Explicit legacy default")


class TestLazyLifecycle(unittest.TestCase):
    """Wrapper lifecycle: methods must always reflect the current parameter
    dict, with construction = validity and set()/compute() as legacy shims.
    See docs/superpowers/specs/2026-07-03-classy-lazy-wrapper-design.md.

    """

    def _fresh_angular_distance(self, h, z=1.0):
        cosmo = Class({'h': h})
        try:
            return cosmo.angular_distance(z)
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()

    def test_set_without_compute_returns_fresh_results(self):
        # The historic footgun: set() then a method call WITHOUT compute()
        # must return results for the new parameters, not stale ones.
        cosmo = Class({'h': 0.67})
        try:
            d_before = cosmo.angular_distance(1.0)
            cosmo.set({'h': 0.70})
            d_after = cosmo.angular_distance(1.0)  # no compute() in between
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()
        d_expected = self._fresh_angular_distance(0.70)
        self.assertNotEqual(d_before, d_after)
        self.assertAlmostEqual(d_after / d_expected, 1.0, places=10)

    def test_unread_parameter_raises_at_first_access(self):
        # After set() with a bogus parameter, the first access (not just
        # compute()) must surface the input error instead of silently
        # serving results for the previous parameters.
        cosmo = Class()
        try:
            cosmo.set({'this_parameter_does_not_exist': 1})
            with self.assertRaises(CosmoSevereError):
                cosmo.angular_distance(1.0)
        finally:
            cosmo.empty()
            cosmo.struct_cleanup()

    def test_methods_work_without_compute(self):
        # Construction = validity: no compute() call is ever needed.
        cosmo = Class({'output': 'tCl', 'l_max_scalars': 100})
        try:
            cl = cosmo.raw_cl(100)
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()
        self.assertTrue(np.all(np.isfinite(cl['tt'])))
        self.assertGreater(np.max(np.abs(cl['tt'])), 0.0)

    def test_bad_parameter_raises_at_construction(self):
        # Constructor-style use validates input at construction (parse-time).
        with self.assertRaises(CosmoSevereError):
            Class({'this_parameter_does_not_exist': 1})

    def test_unread_parameter_error_names_it(self):
        # The error comes from C++ (InputModule::RejectUnreadParameters), the same check
        # the binary runs, and it names the parameter and the way out.
        with self.assertRaisesRegex(CosmoSevereError,
                                    "this_parameter_does_not_exist.*allow_unread_parameters"):
            Class({'this_parameter_does_not_exist': 1})

    def test_allow_unread_parameters_runs_anyway(self):
        cosmo = Class({'this_parameter_does_not_exist': 1, 'allow_unread_parameters': 'yes',
                       'output': 'tCl', 'l_max_scalars': 100})
        try:
            self.assertGreater(np.max(np.abs(cosmo.raw_cl(100)['tt'])), 0.0)
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()

    def test_unread_parameter_raises_at_compute_not_set(self):
        # MontePython compatibility: set() never raises; the error surfaces
        # at compute(), inside the caller's try block.
        cosmo = Class()
        try:
            cosmo.set({'this_parameter_does_not_exist': 1})  # must not raise
            with self.assertRaises(CosmoSevereError):
                cosmo.compute()
        finally:
            cosmo.empty()
            cosmo.struct_cleanup()

    def test_montepython_flow(self):
        # The canonical MontePython sequence: repeated set(), compute(),
        # likelihood-style method calls, then a second point.
        cosmo = Class()
        try:
            cosmo.set({'output': 'tCl', 'l_max_scalars': 100})
            cosmo.set({'h': 0.70})
            cosmo.compute()
            cl1 = cosmo.raw_cl(100)
            self.assertTrue(cosmo.state)
            self.assertGreater(np.max(np.abs(cl1['tt'])), 0.0)

            cosmo.struct_cleanup()
            cosmo.set({'h': 0.68})
            cosmo.compute()
            cl2 = cosmo.raw_cl(100)
            # atol=0: raw dimensionless Cl values are ~1e-10, far below
            # np.allclose's default atol=1e-8, which would compare any two
            # Cl arrays as "close".
            self.assertFalse(np.allclose(cl1['tt'][2:], cl2['tt'][2:],
                                         rtol=1e-5, atol=0.0))
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()

    def test_error_is_sticky_after_failed_compute(self):
        # If reset() clears parameters_changed before the unread-parameter
        # check, a caught error leaves the gate believing it is clean while
        # _thisptr holds a Cosmology that silently ignored the bogus parameter.
        # The second access must re-raise rather than serving stale results.
        cosmo = Class()
        try:
            cosmo.set({'this_parameter_does_not_exist': 1})
            with self.assertRaises(CosmoSevereError):
                cosmo.compute()
            # Gate must remain dirty: subsequent access must still raise.
            with self.assertRaises(CosmoSevereError):
                cosmo.angular_distance(1.0)
            # Recovery: clear the bad parameter, set a good one, verify results.
            cosmo.empty()
            cosmo.set({'h': 0.70})
            d = cosmo.angular_distance(1.0)
            self.assertTrue(np.isfinite(d) and d > 0)
        finally:
            cosmo.empty()
            cosmo.struct_cleanup()

    def test_no_stale_serve_after_failed_rebuild(self):
        # If reset() clears parameters_changed before the Cosmology constructor
        # runs, a thrown constructor error leaves the gate believing it is clean
        # while _thisptr still holds the previous Cosmology. The next access
        # must re-raise rather than silently returning d1.
        cosmo = Class({'h': 0.67})
        try:
            d1 = cosmo.angular_distance(1.0)
            # Adding 100*theta_s while h is already in the dict makes the
            # Cosmology constructor throw (cannot specify both).
            cosmo.set({'100*theta_s': 1.042})
            with self.assertRaises(CosmoSevereError):
                cosmo.compute()
            # Gate must remain dirty: subsequent access must still raise.
            with self.assertRaises(CosmoSevereError):
                cosmo.angular_distance(1.0)
        finally:
            cosmo.empty()
            cosmo.struct_cleanup()

    def test_pars_property_mutation_is_inert(self):
        # pars exposes a snapshot, not the live dict: mutating the returned
        # dict must neither change the cosmology nor leak into the internal
        # parameter dict (set() is the only mutation API — anything else
        # would bypass the parameters_changed gate). The constructor must
        # likewise not alias or pollute the caller's dict.
        caller_dict = {'h': 0.67}
        cosmo = Class(caller_dict)
        try:
            d1 = cosmo.angular_distance(1.0)
            cosmo.pars['h'] = 0.99
            self.assertEqual(cosmo.pars['h'], 0.67)
            self.assertAlmostEqual(cosmo.angular_distance(1.0) / d1, 1.0,
                                   places=12)
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()
        self.assertEqual(caller_dict, {'h': 0.67})


class TestCobayaCompatibility(unittest.TestCase):
    """What Cobaya's `classy` theory (cobaya/theories/classy/classy.py) needs
    from the wrapper, pinned without importing Cobaya. Each of these stopped a
    Cobaya run cold; test_cobaya.py drives Cobaya itself end to end.

    """

    def test_version_is_exposed_and_matches_pyproject(self):
        # Cobaya refuses a classy without __version__, or older than v3.3.3.
        import classy
        pyproject = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 os.pardir, 'pyproject.toml')
        with open(pyproject) as f:
            expected = re.search(r'^version\s*=\s*"([^"]+)"', f.read(),
                                 re.M).group(1)
        self.assertEqual(classy.__version__, expected)
        parts = tuple(int(p) for p in classy.__version__.split('.')[:3])
        self.assertGreaterEqual(parts, (3, 3, 3))

    def test_set_accepts_keyword_arguments(self):
        # Cobaya calls classy.set(**args).
        by_keyword, by_dict = Class(), Class()
        try:
            by_keyword.set(h=0.70)
            by_dict.set({'h': 0.70})
            self.assertEqual(by_keyword.angular_distance(1.0),
                             by_dict.angular_distance(1.0))
        finally:
            for cosmo in (by_keyword, by_dict):
                cosmo.struct_cleanup()
                cosmo.empty()

    def test_set_merges_dict_and_keywords(self):
        # A CLASS key that is not a Python identifier can only travel in the
        # dict, so both forms must be usable in one call.
        merged, reference = Class(), Class({'h': 0.70, 'Omega_k': 0.01})
        try:
            merged.set({'Omega_k': 0.01}, h=0.70)
            self.assertEqual(merged.angular_distance(1.0),
                             reference.angular_distance(1.0))
        finally:
            for cosmo in (merged, reference):
                cosmo.struct_cleanup()
                cosmo.empty()

    def test_non_linear_underscore_spelling_is_read(self):
        # class_public v3 spells the key non_linear, and Cobaya always sends
        # that spelling. The legacy 'non linear' must keep working.
        base = {'output': 'mPk', 'P_k_max_1/Mpc': 3.0}
        underscore = Class(dict(base, non_linear='halofit'))
        legacy = Class(dict(base, **{'non linear': 'halofit'}))
        try:
            pk_underscore = underscore.pk(1.0, 0.0)
            self.assertEqual(pk_underscore, legacy.pk(1.0, 0.0))
            # and the correction is really applied, not silently dropped
            self.assertGreater(pk_underscore / underscore.pk_lin(1.0, 0.0), 1.5)
        finally:
            for cosmo in (underscore, legacy):
                cosmo.struct_cleanup()
                cosmo.empty()

    def test_both_non_linear_spellings_is_an_error(self):
        with self.assertRaises(CosmoSevereError):
            Class({'output': 'mPk', 'non_linear': 'halofit',
                   'non linear': 'hmcode'})

    # Cobaya rejects the point on CosmoComputationError and aborts the run on
    # CosmoSevereError, so a module method that throws must arrive as one of
    # the two, not as the ValueError/RuntimeError of a plain `except +`.

    def test_numeric_failure_in_a_module_method_rejects_the_point(self):
        # class_test (std::runtime_error) in NonlinearModule::nonlinear_pk_at_k_and_z
        cosmo = Class({'output': 'mPk', 'P_k_max_1/Mpc': 1.0})
        try:
            with self.assertRaises(CosmoComputationError):
                cosmo.pk(100.0, 0.0)
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()

    def test_structural_failure_in_a_module_method_is_severe(self):
        # class_test_severe (std::invalid_argument) in LensingModule::cl_output
        cosmo = Class({'output': 'tCl,lCl', 'lensing': 'yes',
                       'l_max_scalars': 300})
        try:
            with self.assertRaises(CosmoSevereError):
                cosmo.lensed_cl(-5)
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()


class TestQueriesTheRunCannotAnswer(unittest.TestCase):
    """Asking for something the run did not compute is a CosmoSevereError.

    The C++ method refuses, and its exception crosses into Python with its
    severity intact, so the wrapper needs no guard of its own. Without the C++
    check these reads index empty tables.
    """

    COSMOLOGIES = {
        'no_pk': {'output': 'tCl', 'l_max_scalars': 300},
        'pk': {'output': 'mPk', 'P_k_max_1/Mpc': 1.0},
        'lensed': {'output': 'tCl,lCl', 'lensing': 'yes', 'l_max_scalars': 300},
        'halofit_for_lensing': {'output': 'tCl,lCl', 'lensing': 'yes', 'non_linear': 'halofit',
                                'l_max_scalars': 300},
        'sigma8_input_no_pk': {'output': 'tCl', 'sigma8': 0.8, 'l_max_scalars': 300},
    }

    @classmethod
    def setUpClass(cls):
        cls.cosmo = {name: Class(params) for name, params in cls.COSMOLOGIES.items()}

    @classmethod
    def tearDownClass(cls):
        for cosmo in cls.cosmo.values():
            cosmo.struct_cleanup()
            cosmo.empty()

    @parameterized.expand([
        ('pk_without_mPk', 'no_pk', lambda c: c.pk(0.1, 0.)),
        ('pk_grid_without_mPk', 'no_pk',
         lambda c: c.get_pk(np.full((1, 1, 1), 0.1), np.zeros(1), 1, 1, 1)),
        ('sigma_without_mPk', 'no_pk', lambda c: c.sigma(8., 0.)),
        ('sigma8_without_mPk', 'no_pk', lambda c: c.sigma8()),
        ('sigma8_cb_without_mPk', 'no_pk', lambda c: c.sigma8_cb()),
        ('sigma8_input_without_mPk', 'sigma8_input_no_pk', lambda c: c.A_s()),
        ('non_linear_pk_array_without_a_method', 'pk',
         lambda c: c.get_pk_array(np.array([0.1]), np.zeros(1), 1, 1, 1)),
        ('pk_tilt_without_mPk', 'no_pk', lambda c: c.pk_tilt(0.1, 0.)),
        ('pk_tilt_beyond_the_k_table', 'pk', lambda c: c.pk_tilt(100., 0.)),
        ('raw_cl_without_Cls', 'pk', lambda c: c.raw_cl()),
        ('raw_cl_no_copy_without_Cls', 'pk', lambda c: c.raw_cl_no_copy()),
        ('lensed_cl_without_Cls', 'pk', lambda c: c.lensed_cl()),
        ('lensed_cl_without_lensing', 'no_pk', lambda c: c.lensed_cl()),
        ('raw_cl_above_l_max', 'lensed', lambda c: c.raw_cl(801)),
        ('raw_cl_no_copy_above_l_max', 'lensed', lambda c: c.raw_cl_no_copy(801)),
        ('lensed_cl_above_l_max', 'lensed', lambda c: c.lensed_cl(301)),
    ])
    def test_is_severe(self, name, cosmology, query):
        with self.assertRaises(CosmoSevereError):
            query(self.cosmo[cosmology])

    def test_sigma8_needs_a_matter_spectrum_not_the_mPk_output(self):
        # Halofit for lensing computes P_m, and sigma8 with it, without 'mPk'.
        self.assertGreater(self.cosmo['halofit_for_lensing'].sigma8(), 0.)


LSS_BASE = {
    'omega_b': 0.02237,
    'omega_cdm': 0.1200,
    'h': 0.6736,
    'A_s': 2.1e-9,
    'n_s': 0.9649,
    'tau_reio': 0.0544,
    'P_k_max_1/Mpc': 3.0,
    'z_max_pk': 3.0,
    # halofit needs k_NL well above P_k_max to reach the top of the z grid; the
    # extension costs integrand points, not perturbation modes.
    'nonlinear_min_k_max': 30.0,
}


class TestLargeScaleStructureProducts(unittest.TestCase):
    """The grid and growth products Cobaya's LSS requirements call for (#435).

    Everything is checked against the wrapper's own scalar calls on the same
    instance, so these pin the plumbing rather than agreement with another code.
    """

    @classmethod
    def setUpClass(cls):
        cls.cosmo = Class(dict(LSS_BASE, output='mPk,dTk,vTk',
                               non_linear='halofit'))
        cls.h = cls.cosmo.h()
        cls.pk, cls.k, cls.z = cls.cosmo.get_pk_and_k_and_z(nonlinear=True)
        cls.pk_lin, _, _ = cls.cosmo.get_pk_and_k_and_z(nonlinear=False)
        cls.tk, cls.k_tk, cls.z_tk = cls.cosmo.get_transfer_and_k_and_z()

    @classmethod
    def tearDownClass(cls):
        cls.cosmo.struct_cleanup()
        cls.cosmo.empty()

    def sample_indices(self):
        return ([0, len(self.k) // 3, len(self.k) // 2, len(self.k) - 1],
                [0, len(self.z) // 2, len(self.z) - 1])

    def test_every_node_of_the_k_axis_is_callable(self):
        # The scalar API used to bound k by exp(ln_k_[-1]) while the grid hands
        # out k_[-1]; exp(log(k)) lands an ulp low, so the topmost node the grid
        # advertised was rejected by pk().
        for k in (self.k[0], self.k[-1]):
            self.assertGreater(self.cosmo.pk(k, 0.0), 0.0)
            self.assertGreater(self.cosmo.pk_lin(k, 0.0), 0.0)

    def test_pk_grid_shape_and_axes(self):
        self.assertEqual(self.pk.shape, (len(self.k), len(self.z)))
        self.assertTrue(np.all(np.diff(self.k) > 0))
        # redshifts run downwards and end exactly at today
        self.assertTrue(np.all(np.diff(self.z) < 0))
        self.assertEqual(self.z[-1], 0.0)
        self.assertTrue(np.all(np.isfinite(self.pk)) and np.all(self.pk > 0))

    def test_pk_grid_matches_the_scalar_calls_at_its_own_nodes(self):
        for ik in self.sample_indices()[0]:
            for iz in self.sample_indices()[1]:
                self.assertAlmostEqual(
                    self.pk[ik, iz] / self.cosmo.pk(self.k[ik], self.z[iz]), 1.0,
                    places=10)
                self.assertAlmostEqual(
                    self.pk_lin[ik, iz] / self.cosmo.pk_lin(self.k[ik], self.z[iz]), 1.0,
                    places=10)

    def test_non_linear_grid_is_boosted_over_the_linear_one(self):
        # otherwise a silently-linear table would pass every other test here
        self.assertGreater(self.pk[-2, -1] / self.pk_lin[-2, -1], 1.5)

    def test_pk_grid_h_units_rescales_k_and_leaves_pk_alone(self):
        pk_h, k_h, z_h = self.cosmo.get_pk_and_k_and_z(nonlinear=True, h_units=True)
        np.testing.assert_allclose(k_h, self.k / self.h, rtol=1e-14)
        np.testing.assert_array_equal(pk_h, self.pk)
        np.testing.assert_array_equal(z_h, self.z)

    def test_non_linear_grid_without_a_non_linear_method_is_an_error(self):
        cosmo = Class(dict(LSS_BASE, output='mPk'))
        try:
            cosmo.get_pk_and_k_and_z(nonlinear=False)  # linear is fine
            with self.assertRaises(CosmoSevereError):
                cosmo.get_pk_and_k_and_z(nonlinear=True)
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()

    def test_grid_without_a_redshift_range_is_an_error(self):
        # no z_max_pk: the table is stored at z=0 only, so there is no grid
        cosmo = Class({'output': 'mPk', 'P_k_max_1/Mpc': 1.0})
        try:
            with self.assertRaises(CosmoSevereError):
                cosmo.get_pk_and_k_and_z(nonlinear=False)
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()

    def test_only_clustering_species_excludes_massive_neutrinos(self):
        cosmo = Class(dict(LSS_BASE, output='mPk', N_ur=2.0308,
                           **{'nu1.type': 'ncdm_standard', 'nu1.m': 0.3}))
        try:
            total, k, z = cosmo.get_pk_and_k_and_z(nonlinear=False)
            clustering, _, _ = cosmo.get_pk_and_k_and_z(nonlinear=False,
                                                        only_clustering_species=True)
            for ik in (0, len(k) // 2, len(k) - 2):
                self.assertAlmostEqual(
                    total[ik, -1] / cosmo.pk_lin(k[ik], z[-1]), 1.0, places=10)
                self.assertAlmostEqual(
                    clustering[ik, -1] / cosmo.pk_cb_lin(k[ik], z[-1]), 1.0, places=10)
            # free-streaming neutrinos suppress the total spectrum below P_cb
            self.assertGreater(clustering[-2, -1] / total[-2, -1], 1.01)
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()

    def test_get_pk_array_follows_its_nonlinear_argument(self):
        cosmo = Class(dict(LSS_BASE, output='mPk', non_linear='halofit'))
        k, z = np.array([0.01, 0.1, 1.0]), np.array([0., 0.5])
        try:
            for nonlinear, pk in ((0, cosmo.pk_lin), (1, cosmo.pk)):
                with self.subTest(nonlinear=nonlinear):
                    np.testing.assert_allclose(
                        cosmo.get_pk_array(k, z, len(k), len(z), nonlinear),
                        [pk(k_i, z_j) for z_j in z for k_i in k], rtol=1e-3)
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()

    def test_cb_spectra_without_warm_matter_are_the_total_matter_spectra(self):
        # delta_cb is the cold part of delta_m, so without warm matter the two
        # are one spectrum, which the module stores once, as P_m.
        cosmo = Class(dict(LSS_BASE, output='mPk', non_linear='halofit'))
        k, z = 0.1, 0.5
        k_grid, z_grid = np.full((1, 1, 1), k), np.full(1, z)
        k_vec, z_vec = np.array([0.01, 0.1, 1.0]), np.array([0., 0.5])
        pairs = {
            'pk': (lambda c: c.pk_cb(k, z), lambda c: c.pk(k, z)),
            'pk_lin': (lambda c: c.pk_cb_lin(k, z), lambda c: c.pk_lin(k, z)),
            'get_pk': (lambda c: c.get_pk_cb(k_grid, z_grid, 1, 1, 1),
                       lambda c: c.get_pk(k_grid, z_grid, 1, 1, 1)),
            'get_pk_lin': (lambda c: c.get_pk_cb_lin(k_grid, z_grid, 1, 1, 1),
                           lambda c: c.get_pk_lin(k_grid, z_grid, 1, 1, 1)),
            'get_pk_array': (lambda c: c.get_pk_cb_array(k_vec, z_vec, 3, 2, 0),
                             lambda c: c.get_pk_array(k_vec, z_vec, 3, 2, 0)),
            'sigma': (lambda c: c.sigma_cb(8., z), lambda c: c.sigma(8., z)),
            'sigma8': (lambda c: c.sigma8_cb(), lambda c: c.sigma8()),
        }
        try:
            for name, (cb, m) in pairs.items():
                with self.subTest(name):
                    np.testing.assert_array_equal(cb(cosmo), m(cosmo))
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()

    def test_sigma_h_units_is_the_same_radius_in_other_units(self):
        for z in (0.0, 1.0):
            self.assertAlmostEqual(
                self.cosmo.sigma(8.0, z, h_units=True) / self.cosmo.sigma(8.0 / self.h, z),
                1.0, places=12)
        self.assertAlmostEqual(
            self.cosmo.sigma(8.0, 0.0, h_units=True) / self.cosmo.sigma8(), 1.0, places=12)

    def test_effective_f_sigma8_is_its_own_finite_difference(self):
        for z, step in ((1.5, 0.1), (0.5, 0.1), (0.05, 0.1), (0.0, 0.1)):
            got = self.cosmo.effective_f_sigma8(z, z_step=step)
            if z < step / 10.0:
                step = step / 10.0
                want = (self.cosmo.sigma(8, z, h_units=True)
                        - self.cosmo.sigma(8, z + step, h_units=True)) / step * (1 + z)
            else:
                step = min(step, z) if z < step else step
                want = (self.cosmo.sigma(8, z - step, h_units=True)
                        - self.cosmo.sigma(8, z + step, h_units=True)) / (2 * step) * (1 + z)
            self.assertAlmostEqual(got / want, 1.0, places=12)

    def test_effective_f_sigma8_tracks_f_times_sigma8(self):
        # a loose check that it is fsigma8 and not, say, its negative
        for z in (0.0, 0.5, 1.5):
            approx = (self.cosmo.scale_independent_growth_factor_f(z)
                      * self.cosmo.sigma(8, z, h_units=True))
            self.assertAlmostEqual(self.cosmo.effective_f_sigma8(z) / approx,
                                   1.0, places=2)

    def test_angular_distance_from_to(self):
        # flat space: d_A(z1,z2) = (chi2 - chi1)/(1 + z2)
        chi = self.cosmo.z_of_r(np.array([0.5, 1.5]))[0]
        self.assertAlmostEqual(
            self.cosmo.angular_distance_from_to(0.5, 1.5) / ((chi[1] - chi[0]) / 2.5),
            1.0, places=12)
        # from here it is the ordinary angular diameter distance
        self.assertAlmostEqual(
            self.cosmo.angular_distance_from_to(0.0, 1.5) / self.cosmo.angular_distance(1.5),
            1.0, places=10)
        # Cobaya hands it unordered pairs and expects zero for the wrong order
        self.assertEqual(self.cosmo.angular_distance_from_to(1.5, 0.5), 0.0)
        self.assertEqual(self.cosmo.angular_distance_from_to(1.5, 1.5), 0.0)

    def test_transfer_grid_shares_the_axes_of_the_pk_grid(self):
        np.testing.assert_array_equal(self.k_tk, self.k)
        np.testing.assert_array_equal(self.z_tk, self.z)
        for name in ('d_m', 'd_tot', 'phi', 'psi', 'd_b', 'd_cdm'):
            self.assertIn(name, self.tk)
            self.assertEqual(self.tk[name].shape, (len(self.k), len(self.z)))
        self.assertNotIn('k (h/Mpc)', self.tk)

    def test_transfer_grid_matches_get_transfer_at_a_node(self):
        index_z = 3
        one = self.cosmo.get_transfer(self.z[index_z])
        for name in ('d_m', 'phi', 'psi'):
            np.testing.assert_allclose(self.tk[name][:, index_z], one[name], rtol=1e-9)

    def test_d_m_is_the_contrast_the_matter_spectrum_is_built_from(self):
        # P_lin(k,z) = P_primordial(k) * d_m(k,z)^2, so the ratio must not depend
        # on z. That pins d_m to the right quantity without fixing a convention.
        ratio = self.pk_lin / self.tk['d_m'] ** 2
        for index_z in range(len(self.z)):
            np.testing.assert_allclose(ratio[:, index_z], ratio[:, -1], rtol=1e-6)

    def test_weyl_grid_is_the_matter_grid_rescaled(self):
        weyl, k, z = self.cosmo.get_Weyl_pk_and_k_and_z(nonlinear=False)
        np.testing.assert_array_equal(k, self.k)
        np.testing.assert_array_equal(z, self.z)
        expected = (self.pk_lin
                    * ((self.tk['phi'] + self.tk['psi']) / 2.0 / self.tk['d_m']) ** 2
                    * self.k[:, None] ** 4)
        np.testing.assert_allclose(weyl, expected, rtol=1e-12)
        self.assertTrue(np.all(np.isfinite(weyl)))

    def test_weyl_grid_h_units_carries_through_the_k4(self):
        weyl, k, _ = self.cosmo.get_Weyl_pk_and_k_and_z(nonlinear=False)
        weyl_h, k_h, _ = self.cosmo.get_Weyl_pk_and_k_and_z(nonlinear=False, h_units=True)
        np.testing.assert_allclose(k_h, k / self.h, rtol=1e-14)
        np.testing.assert_allclose(weyl_h, weyl / self.h ** 4, rtol=1e-12)

    def test_transfer_grid_without_transfers_is_an_error(self):
        cosmo = Class(dict(LSS_BASE, output='mPk'))
        try:
            with self.assertRaises(CosmoSevereError):
                cosmo.get_transfer_and_k_and_z()
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()

    def test_transfer_grid_without_a_redshift_range_is_an_error(self):
        # No z_max_pk means no late-time sampling, and late_sources_ is null in
        # that case -- this has to raise rather than reach the node reader.
        cosmo = Class({'output': 'mPk,dTk', 'P_k_max_1/Mpc': 1.0})
        try:
            with self.assertRaises(CosmoSevereError):
                cosmo.get_transfer_and_k_and_z()
            # the z=0 path stays available
            self.assertIn('d_m', cosmo.get_transfer(0.0))
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()



class TestTabulatedSelection(unittest.TestCase):
    """selection = tabulated: the dN/dz of each bin read in a file.

    The file has the redshift in the first column and one column per bin. The
    reference is the same cosmology with analytic gaussian bins, tabulated here
    on a grid fine enough for the spline to reproduce them.
    """

    MEANS = (0.5, 1.0)
    WIDTHS = (0.07, 0.1)
    BASE = {'output': 'nCl', 'number count contributions': 'density',
            'l_max_lss': 300, 'non_diagonal': 1}

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.cosmo = Class()

    def tearDown(self):
        self.cosmo.struct_cleanup()
        self.cosmo.empty()
        shutil.rmtree(self.tmp)

    def write(self, name, text):
        path = os.path.join(self.tmp, name)
        with open(path, 'w') as f:
            f.write(text)
        return path

    def gaussian_file(self):
        """The gaussian bins of CLASS, cut at 5 sigma like selection_cut_at_sigma."""
        z = np.arange(0., 2.0001, 0.002)
        columns = [z]
        for mean, width in zip(self.MEANS, self.WIDTHS):
            x = (z - mean) / width
            columns.append(np.where(np.abs(x) <= 5., np.exp(-0.5 * x**2), 0.))
        lines = ['# z  dN/dz_1  dN/dz_2']
        lines += [' '.join('%.10e' % c[i] for c in columns) for i in range(len(z))]
        return self.write('bins.dat', '\n'.join(lines) + '\n')

    def density_cl(self, params):
        self.cosmo.set(dict(self.BASE, **params))
        self.cosmo.compute()
        cl = self.cosmo.raw_cl(300)
        self.cosmo.struct_cleanup()
        self.cosmo.empty()
        return {key: cl[key][2:] for key in cl if key.startswith('dens')}

    def test_tabulated_bins_reproduce_gaussian_bins(self):
        gaussian = self.density_cl({'selection': 'gaussian',
                                    'selection_mean': '0.5,1.0',
                                    'selection_width': '0.07,0.1'})
        tabulated = self.density_cl({'selection': 'tabulated',
                                     'selection_file': self.gaussian_file()})
        for key in ('dens[1]-dens[1]', 'dens[2]-dens[2]', 'dens[1]-dens[2]'):
            scale = np.max(np.abs(gaussian[key]))
            error = np.max(np.abs(tabulated[key] - gaussian[key])) / scale
            self.assertLess(error, 1e-3, key)

    def test_bias_of_each_tabulated_bin_scales_its_spectra(self):
        path = self.gaussian_file()
        unbiased = self.density_cl({'selection': 'tabulated', 'selection_file': path})
        biased = self.density_cl({'selection': 'tabulated', 'selection_file': path,
                                  'selection_bias': '1.,2.'})
        np.testing.assert_allclose(biased['dens[1]-dens[1]'], unbiased['dens[1]-dens[1]'], rtol=1e-10)
        np.testing.assert_allclose(biased['dens[1]-dens[2]'], 2. * unbiased['dens[1]-dens[2]'], rtol=1e-10)
        np.testing.assert_allclose(biased['dens[2]-dens[2]'], 4. * unbiased['dens[2]-dens[2]'], rtol=1e-10)

    def test_malformed_selection_files_are_rejected(self):
        cases = {
            'missing': (None, 'could not open'),
            'single_column': ('0.0\n0.1\n0.2\n', 'at least two columns'),
            'one_line': ('0.0 1.0 1.0\n', 'at least two lines'),
            'ragged': ('0.0 0.0 0.0\n0.1 1.0\n0.2 0.0 0.0\n', 'contains 2 numbers'),
            'unreadable': ('0.0 0.0 0.0\n0.1 1.0 one\n0.2 0.0 0.0\n', 'could not read'),
            'decreasing_z': ('0.0 0.0 0.0\n0.2 1.0 1.0\n0.1 0.0 0.0\n', 'strictly growing'),
            'negative_z': ('-0.1 0.0 0.0\n0.1 1.0 1.0\n0.2 0.0 0.0\n', 'is negative'),
            'negative_dNdz': ('0.0 0.0 0.0\n0.1 -1.0 1.0\n0.2 0.0 0.0\n', 'negative at z'),
            'empty_bin': ('0.0 0.0 0.0\n0.1 1.0 0.0\n0.2 0.0 0.0\n', 'vanishes everywhere'),
            'unordered_bins': ('0.0 0.0 0.0\n0.1 0.0 1.0\n0.2 1.0 0.0\n0.3 0.0 0.0\n',
                               'growing mean redshift'),
            'infinite': ('0.0 0.0 0.0\n0.1 inf 1.0\n0.2 0.0 0.0\n', 'not a finite number'),
            'nan': ('0.0 0.0 0.0\n0.1 nan 1.0\n0.2 0.0 0.0\n', 'not a finite number'),
            'overflow': ('0.0 0.0 0.0\n0.1 1e999 1.0\n0.2 0.0 0.0\n', 'not a finite number'),
        }
        for name, (text, message) in cases.items():
            path = (os.path.join(self.tmp, 'nowhere.dat') if text is None
                    else self.write(name + '.dat', text))
            self.cosmo.set(dict(self.BASE, selection='tabulated', selection_file=path))
            with self.assertRaises(CosmoSevereError, msg=name) as ctx:
                self.cosmo.compute()
            self.assertIn(message, str(ctx.exception), name)
            self.cosmo.struct_cleanup()
            self.cosmo.empty()

    def test_bin_with_weight_only_at_z0_has_a_positive_mean(self):
        """The first moment of the linear interpolant is z1/3 > 0. The trapezoidal rule
        applied to z dN/dz gave a mean of exactly zero, and k_max = l_max/(tau0 - tau(mean))
        diverged: the run never finished."""
        path = self.write('z0.dat', '0 1\n0.1 0\n')
        self.cosmo.set({'output': 'nCl', 'number count contributions': 'density',
                        'l_max_lss': 100, 'selection': 'tabulated', 'selection_file': path})
        self.cosmo.compute()
        cl = self.cosmo.raw_cl(100)['dens[1]-dens[1]'][2:]
        self.assertTrue(np.all(np.isfinite(cl)))
        self.assertTrue(np.all(cl > 0.))

    def test_inconsistent_selection_inputs_are_rejected(self):
        path = self.gaussian_file()
        cases = {
            'no_file': ({'selection': 'tabulated'}, "selection_file"),
            'with_mean': ({'selection': 'tabulated', 'selection_file': path,
                           'selection_mean': '0.5,1.0'}, 'do not pass'),
            'with_dNdz': ({'selection': 'tabulated', 'selection_file': path,
                           'dNdz_selection': 'analytic'}, 'dNdz_selection'),
        }
        for name, (params, message) in cases.items():
            self.cosmo.set(dict(self.BASE, **params))
            with self.assertRaises(CosmoSevereError, msg=name) as ctx:
                self.cosmo.compute()
            self.assertIn(message, str(ctx.exception), name)
            self.cosmo.struct_cleanup()
            self.cosmo.empty()

    def test_out_of_range_sampling_parameters_are_rejected(self):
        """Range checks reject the point (CosmoComputationError). A negative lensing
        factor and fraction used to pass the check on their product and loop forever,
        and so did a selection_mean of zero (k_max diverges)."""
        cases = {
            'negative_lensing': ({'output': 'tCl,lCl', 'lensing': 'yes',
                                  'cmb_lensing_sampling_factor': -1.,
                                  'cmb_lensing_sampling_fraction': -0.001},
                                 'cmb_lensing_sampling_factor'),
            'lensing_fraction': ({'output': 'tCl,lCl', 'lensing': 'yes',
                                  'cmb_lensing_sampling_fraction': 1.5},
                                 'cmb_lensing_sampling_fraction'),
            'mean_at_z0': (dict(self.BASE, selection='gaussian', selection_mean='0.',
                                selection_width='0.1', non_diagonal=0),
                           'mean redshift'),
        }
        for name, (params, message) in cases.items():
            self.cosmo.set(params)
            with self.assertRaises(CosmoComputationError, msg=name) as ctx:
                self.cosmo.compute()
            self.assertIn(message, str(ctx.exception), name)
            self.cosmo.struct_cleanup()
            self.cosmo.empty()


class TestHMcode(unittest.TestCase):
    """HMcode 2016 and 2020 (Mead et al. 2021), with and without baryonic feedback.

    The references are the boosts P_nl/P_lin of CAMB 2.0.4 (halofit_version mead2016,
    mead2020 and mead2020_feedback), which runs Mead's own HMcode, for the same cosmology.
    The boost cancels most of the differences between the linear spectra of the two codes;
    what remains is below 1e-3 at these k, for both versions. The tolerance of 3e-3 is ten
    times smaller than the differences between the versions at k = 0.1, 1 and 3/Mpc.
    """

    BASE = {'output': 'mPk', 'non_linear': 'hmcode', 'h': 0.6736, 'omega_b': 0.02237,
            'omega_cdm': 0.12, 'A_s': 2.1e-9, 'n_s': 0.9649, 'tau_reio': 0.0544,
            'N_ur': 3.044, 'P_k_max_1/Mpc': 20., 'z_max_pk': 1.}
    K = (0.1, 0.3, 1., 3.)  # 1/Mpc
    CAMB = {
        'mead2016': {0.: (1.11534, 2.10407, 9.63662, 27.18714),
                     1.: (1.04123, 1.50186, 5.03742, 16.37017)},
        'mead2020': {0.: (1.05982, 2.10976, 9.35522, 26.20411),
                     1.: (1.04047, 1.51836, 5.01971, 16.49750)},
        'mead2020_feedback, log10T_AGN = 8.3': {0.: (1.04300, 1.94759, 7.63480, 19.30947),
                                                1.: (1.03575, 1.45914, 4.24738, 12.87991)},
        # dark energy w0 = -0.9, wa = 0.1 (PPF in CAMB) and curvature Omega_k = 0.02, which
        # differ through the growth and the Dolag correction of the concentrations
        'mead2016, w0wa': {0.: (1.10864, 2.04693, 9.35441, 27.38723),
                           1.: (1.03883, 1.48323, 4.92397, 16.52766)},
        'mead2016, Omega_k': {0.: (1.11204, 2.07562, 9.49323, 27.22785),
                              1.: (1.03969, 1.48970, 4.95468, 16.35080)},
        'mead2020, w0wa': {0.: (1.06014, 2.06238, 9.21963, 26.97426),
                           1.: (1.03954, 1.50009, 4.90261, 16.63141)},
        'mead2020, Omega_k': {0.: (1.05990, 2.08530, 9.26173, 26.39346),
                              1.: (1.03992, 1.50672, 4.93988, 16.50623)},
    }
    W0WA = {'Omega_Lambda': 0., 'w0_fld': -0.9, 'wa_fld': 0.1}
    OMEGA_K = {'Omega_k': 0.02}

    def setUp(self):
        self.cosmo = Class()

    def tearDown(self):
        self.cosmo.struct_cleanup()
        self.cosmo.empty()

    def boost(self, params, k=K, z=(0., 1.)):
        self.cosmo.struct_cleanup()
        self.cosmo.empty()
        self.cosmo.set(dict(self.BASE, **params))
        self.cosmo.compute()
        return {zz: np.array([self.cosmo.pk(kk, zz) / self.cosmo.pk_lin(kk, zz) for kk in k])
                for zz in z}

    def assert_matches_camb(self, params, reference):
        boost = self.boost(params)
        for z, expected in self.CAMB[reference].items():
            np.testing.assert_allclose(boost[z], expected, rtol=3e-3,
                                       err_msg=f'{reference} at z = {z}')

    def test_2016_matches_camb(self):
        self.assert_matches_camb({'hmcode_version': '2016'}, 'mead2016')

    def test_2020_matches_camb(self):
        self.assert_matches_camb({'hmcode_version': '2020'}, 'mead2020')

    def test_2020_baryonic_feedback_matches_camb(self):
        self.assert_matches_camb({'hmcode_version': '2020_baryonic_feedback',
                                  'log10T_heat_hmcode': 8.3},
                                 'mead2020_feedback, log10T_AGN = 8.3')

    def test_dark_energy_and_curvature_match_camb(self):
        """The Dolag correction of the concentrations, in both versions: it applies to every
        halo and in every cosmology, relative to a flat LCDM with the same Omega_m."""
        for version in ('2016', '2020'):
            for name, params in (('w0wa', self.W0WA), ('Omega_k', self.OMEGA_K)):
                self.assert_matches_camb(dict(params, hmcode_version=version),
                                         f'mead{version}, {name}')

    def test_default_version_is_2020_and_2015_is_2016(self):
        for pair in (({}, {'hmcode_version': '2020'}),
                     ({'hmcode_version': '2015'}, {'hmcode_version': '2016'})):
            first, second = (self.boost(params) for params in pair)
            for z in first:
                np.testing.assert_array_equal(first[z], second[z], err_msg=str(pair))

    def test_feedback_suppression_grows_with_the_heating_temperature(self):
        """Feedback leaves large scales alone and suppresses small ones, the more so the
        hotter the AGN heating."""
        k = (0.01, 3.)
        fitted = self.boost({}, k=k, z=(0.,))[0.]
        ratios = [self.boost({'hmcode_version': '2020_baryonic_feedback',
                              'log10T_heat_hmcode': t}, k=k, z=(0.,))[0.] / fitted
                  for t in (7.3, 7.8, 8.3)]
        for ratio in ratios:
            self.assertAlmostEqual(ratio[0], 1., delta=1e-3)
        self.assertTrue(1. > ratios[0][1] > ratios[1][1] > ratios[2][1] > 0.5,
                        [r[1] for r in ratios])

    def test_2020_runs_with_massive_neutrinos_and_a_w0wa_fluid(self):
        for params in ({'N_ur': 2.0308, 'N_ncdm': 1, 'm_ncdm': 0.06},
                       {'Omega_Lambda': 0., 'w0_fld': -0.9, 'wa_fld': 0.1}):
            boost = self.boost(params)
            for z in boost:
                self.assertTrue(np.all(np.isfinite(boost[z])), params)
                self.assertGreater(boost[z][2], 3., params)

    def test_inconsistent_hmcode_inputs_are_rejected(self):
        """Structural errors: a version that does not exist, and the inputs of one version
        given to another."""
        cases = {
            'unknown_version': ({'hmcode_version': '2018'}, 'is not one of'),
            'eta_0_with_2020': ({'eta_0': 0.6}, 'parameters of HMcode 2016'),
            'c_min_with_2020': ({'c_min': 3.}, 'parameters of HMcode 2016'),
            'feedback_model_with_2020': ({'feedback model': 'owls_agn'},
                                         'parameters of HMcode 2016'),
            'T_heat_without_feedback': ({'log10T_heat_hmcode': 8.},
                                        'only used with hmcode_version = 2020_baryonic'),
            'T_heat_with_2016': ({'hmcode_version': '2016', 'log10T_heat_hmcode': 8.},
                                 'only used with hmcode_version = 2020_baryonic'),
            'nk_wiggle_with_2016': ({'hmcode_version': '2016', 'nk_wiggle': 256},
                                    'only used by HMcode 2020'),
        }
        for name, (params, message) in cases.items():
            self.cosmo.set(dict(self.BASE, **params))
            with self.assertRaises(CosmoSevereError, msg=name) as ctx:
                self.cosmo.compute()
            self.assertIn(message, str(ctx.exception), name)
            self.cosmo.struct_cleanup()
            self.cosmo.empty()

    def test_unsupported_hmcode_2020_points_are_rejected(self):
        """Range errors: too few wavenumbers for the wiggle split, and dark energy that the
        growth without radiation of HMcode 2020 cannot represent."""
        cases = {
            'nk_wiggle': ({'nk_wiggle': 8}, 'nk_wiggle = 8 is too small'),
            'scalar_field': ({'Omega_scf': 0.1, 'attractor_ic_scf': 'yes',
                              'scf_parameters': '10, 0, 0, 0'},
                             'HMcode 2020 knows dark energy as Lambda and the fluid only'),
        }
        for name, (params, message) in cases.items():
            self.cosmo.set(dict(self.BASE, **params))
            with self.assertRaises(CosmoComputationError, msg=name) as ctx:
                self.cosmo.compute()
            self.assertIn(message, str(ctx.exception), name)
            self.cosmo.struct_cleanup()
            self.cosmo.empty()

    def test_non_finite_linear_spectrum_is_rejected(self):
        """A linear P(k) that is not finite stops the run before the sigma integrals.
        extrapolation_method = 0 (extrap_zero) continues the sources by zero beyond k_max,
        so ln P(k) is -inf there. Without the check HMcode stopped with x_min = nan in its
        search for nu = 1, and halofit, which met the same NaN from the default
        extrapolation at dark acoustic oscillations, could leave out the non-linear
        correction without a word."""
        self.cosmo.set(dict(self.BASE, extrapolation_method=0))
        with self.assertRaises(CosmoComputationError) as ctx:
            self.cosmo.compute()
        self.assertIn('extrapolated beyond k_max', str(ctx.exception))

if __name__ == '__main__':
    toto = TestClass()
    unittest.main()
