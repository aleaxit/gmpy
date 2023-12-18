"""
Run some tests for Cython integration

To run the test just do

    python runtests.py

The test suite consists of three files:

- runtests.py: a master file
- setup_cython.py: the setup file to control Cython compilation
- test_cython.pyx: the cython file which uses the C-API of gmpy2
"""

import gmpy2
import shutil
import subprocess
import sys
import tempfile
import os

try:
    import Cython
except ImportError:
    sys.stderr.write('Cython is not installed... skipping cython tests')
    sys.exit(0)

print()
print("Unit tests for gmpy2 {0} with Cython {1}".format(gmpy2.version(), Cython.__version__))
print()


try:
    old_path = os.getcwd()
    tempdir_path = tempfile.mkdtemp()

    dirname = os.path.dirname(__file__)
    if dirname != '':
        os.chdir(dirname)
    if sys.version_info >= (3, 8):
        shutil.copytree('./', tempdir_path, dirs_exist_ok=True)
    else:
        from distutils.dir_util import copy_tree

        copy_tree('./', tempdir_path)
    os.chdir(tempdir_path)

    if subprocess.call([sys.executable, 'setup_cython.py', 'build_ext', '--inplace']):
        raise SystemExit('compilation failed')


    if subprocess.call([sys.executable, '-c', 'import gmpy2; import test_cython; test_cython.run()']):
        raise SystemExit('cython test failed')

finally:
    os.chdir(old_path)
    shutil.rmtree(tempdir_path)
