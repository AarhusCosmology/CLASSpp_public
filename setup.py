import os
import runpy
from concurrent.futures import ThreadPoolExecutor
import distutils.ccompiler

from pathlib import Path

def _parallel_compile(self, sources, output_dir=None, macros=None,
                      include_dirs=None, debug=0, extra_preargs=None,
                      extra_postargs=None, depends=None):
    """Drop-in replacement for CCompiler.compile that compiles in parallel."""
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(
        output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)

    def _compile_obj(obj):
        try:
            src, ext = build[obj]
        except KeyError as exc:
            raise KeyError(
                "Missing build mapping for object {!r} during compilation".format(obj)
            ) from exc
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)

    try:
        max_jobs = int(os.environ.get('MAX_JOBS', 0))
    except ValueError:
        max_jobs = 0
    max_jobs = max(max_jobs, 0) or os.cpu_count() or 1
    if max_jobs == 1:
        for obj in objects:
            _compile_obj(obj)
    else:
        with ThreadPoolExecutor(max_workers=max_jobs) as executor:
            futures = [executor.submit(_compile_obj, obj) for obj in objects]
            for future in futures:
                future.result()  # re-raises any compilation error

    return objects


distutils.ccompiler.CCompiler.compile = _parallel_compile

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy


def run_generate_wrapper(project_root):
    runpy.run_path(Path(project_root) / 'generate_wrapper.py', run_name='__main__')


def my_cythonize(*args, **kwargs):
    run_generate_wrapper(Path(__file__).resolve().parent)
    return cythonize(*args, **kwargs)


prepend = lambda dir, *fnames: list(map(lambda s: os.path.join(dir, s), fnames))
source_files = (
    prepend(
        'hyrec',
        'helium.c',
        'history.c',
        'hydrogen.c',
        'hyrec.c',
        'hyrectools.c',
    )
    + prepend(
        'source',
        'background_module.cpp',
        'cosmology.cpp',
        'input_module.cpp',
        'lensing_module.cpp',
        'nonlinear_module.cpp',
        'perturbations_module.cpp',
        'primordial_module.cpp',
        'spectra_module.cpp',
        'thermodynamics_module.cpp',
        'transfer_module.cpp',
    )
    + prepend(
        'species',
        'baryons.cpp',
        'base_species.cpp',
        'cdm.cpp',
        'composite_species.cpp',
        'dark_radiation_species.cpp',
        'dcdm.cpp',
        'dcdm_dr_species.cpp',
        'dncdm_dr_species.cpp',
        'dncdm_species.cpp',
        'fluid.cpp',
        'greybody_ncdm_species.cpp',
        'idm_dr_idr_species.cpp',
        'idm_drmd_idr_drmd_species.cpp',
         'lambda.cpp',
         'ncdm_base_species.cpp',
         'ncdm_species.cpp',
         'ncdm_interacting_species.cpp',
         'perturb_column_writer.cpp',
         'photons.cpp',
         'scalar_field.cpp',
         'species_collection.cpp',
         'species_input.cpp',
         'ultra_relativistic.cpp',
        'interacting_species.cpp',
    )
    + prepend(
        'tools',
        'arrays.cpp',
        'common.cpp',
        'dei_rkck.cpp',
        'evolver_ndf15.cpp',
        'evolver_rkck.cpp',
        'evolver_rkdp45.cpp',
        'exceptions.cpp',
        'growTable.cpp',
        'hyperspherical.cpp',
        'parser.cpp',
        'quadrature.cpp',
        'sparse.cpp',
        'trigonometric_integrals.cpp',
    )
)
c_source_files   = [s for s in source_files if s.endswith('.c')]
cpp_source_files = [s for s in source_files if s.endswith('.cpp')]


include_dirs = [numpy.get_include()]
root_folder = '.'
for sub_folder in ['include', 'main', 'source', 'tools', 'species', '.']:
    include_dirs.append(os.path.join(root_folder, sub_folder))

# Define cython extension and fix Python version.
#
# Symbol visibility: on Linux/macOS (gcc/clang) the default for shared
# libraries is to export every symbol, so when both classy and a parallel
# build (e.g. classyref in the lvl2 nose tests) are loaded into the same
# Python process the dynamic linker can collapse duplicate C++ symbols
# (NCDMSpecies::CreateAll, SpeciesBuildContext, etc.) onto whichever .so
# was loaded first — silently mismatching ABI between the two builds and
# segfaulting. Building with -fvisibility=hidden hides everything except
# the Cython PyInit_* entry point (which Cython tags PyMODINIT_FUNC =
# default visibility), so the two .so files stay independent. Windows
# (MSVC) hides DLL symbols by default, so no flag is needed there.
if os.name == 'nt':
    extra_compile_args = ['/std:c++17', '/O2']
else:
    extra_compile_args = ['-std=c++17', '-O3', '-fvisibility=hidden', '-fvisibility-inlines-hidden']

classy_ext = Extension('classy', ['classy.pyx'] + cpp_source_files,
                           include_dirs=include_dirs,
                           libraries=['m'] if not os.name == 'nt' else [],
                           library_dirs=[root_folder],
                           language="c++",
                           extra_compile_args=extra_compile_args,
                           define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
                                          ("__CLASSDIR__", '"{}"'.format(os.path.abspath(root_folder)))]);
myclib = ('myclib', {'sources': c_source_files,
                     'include_dirs':include_dirs})

setup(
    libraries=[myclib],
    ext_modules=my_cythonize(
        classy_ext,
        language_level=3,
        annotate=False,
    ),
    packages=[
        'classy.bbn',
        'classy.hyrec',
    ],
    package_dir={
        'classy': '',
    },
    package_data={
        'classy.bbn': ['*.dat'],
        'classy.hyrec': ['*.dat'],
    },
)
