from __future__ import print_function

import sys

print("\nImportant test information!\n")
print("Please use 'test/runtests.py' to run the new test suite.\n")

if sys.version_info[0] == 3:
    print("Please use 'test3/gmpy_test.py' to run legacy tests with Python 3.x.\n")
else:
    print("Please use 'test2/gmpy_test.py' to run legacy tests with Python 2.x.\n")
