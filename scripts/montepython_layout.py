"""Recreate the python/build/lib.* layout that MontePython >= 3.6 expects.

setuptools used to leave the built extension in build/lib.<platform>-<tag>/,
which `make classy` renamed to lib.<sys.version>.<platform>-<tag>/ (see
https://github.com/brinckmann/montepython_public/issues/371). scikit-build-core
has no such directory, so this script reproduces the same end result from the
module installed in site-packages.
"""
import os
import shutil
import sys
import sysconfig

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Never import a stray in-repo classy/ directory (dev symlink hack).
sys.path = [p for p in sys.path
            if os.path.abspath(p if p else os.getcwd()) != repo_root]
import classy  # noqa: E402

libdir = os.path.join(
    repo_root, 'python', 'build',
    'lib.{}.{}-{}'.format(sys.version, sysconfig.get_platform(),
                          sys.implementation.cache_tag))
os.makedirs(libdir, exist_ok=True)
shutil.copy2(classy.__file__, libdir)
data_dir = os.path.join(os.path.dirname(classy.__file__), 'classy')
shutil.copytree(data_dir, os.path.join(libdir, 'classy'), dirs_exist_ok=True)
print('MontePython layout created at', libdir)
