from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import os

def my_cythonize(*args, **kwargs):
    with open('generate_wrapper.py', 'r', encoding='utf-8') as f:
        exec(f.read())
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
        'cdm.cpp',
        'dark_radiation_species.cpp',
        'dcdm.cpp',
        'fluid.cpp',
        'lambda.cpp',
        'ncdm_species.cpp',
        'photons.cpp',
        'scalar_field.cpp',
        'ultra_relativistic.cpp',
    )
    + prepend(
        'tools',
        'arrays.cpp',
        'common.cpp',
        'dark_radiation.cpp',
        'dei_rkck.cpp',
        'evolver_ndf15.cpp',
        'evolver_rkck.cpp',
        'exceptions.cpp',
        'growTable.cpp',
        'hyperspherical.cpp',
        'non_cold_dark_matter.cpp',
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

# Define cython extension and fix Python version
classy_ext = Extension('classy', ['classy.pyx'] + cpp_source_files,
                           include_dirs=include_dirs,
                           libraries=['m'] if not os.name == 'nt' else [],
                           library_dirs=[root_folder],
                           language="c++",
                           extra_compile_args=(['-std=c++17'] if os.name != 'nt' else ['/std:c++17']),
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
