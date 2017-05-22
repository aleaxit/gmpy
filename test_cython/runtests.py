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
from distutils.dir_util import copy_tree

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

    os.chdir(os.path.dirname(__file__))
    copy_tree('./', tempdir_path)
    os.chdir(tempdir_path)

    if subprocess.call([sys.executable, 'setup_cython.py', 'build_ext', '--inplace']):
        raise SystemExit('compilation failed')

    if subprocess.call([sys.executable, '-c', 'import test_cython']):
        raise SystemExit('cython test failed')

finally:
    os.chdir(old_path)
    shutil.rmtree(tempdir_path)
