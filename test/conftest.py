import pytest

import gmpy2


collect_ignore_glob = ['*.txt']


@pytest.fixture(autouse=True, scope='function')
def _set_default_context():
    gmpy2.set_context(gmpy2.context())


def pytest_report_header(config):
    print("""
  Mutliple-precision library:     {0}
  Floating-point library:         {1}
  Complex library:                {2}
""".format(gmpy2.mp_version(),
           gmpy2.mpfr_version(),
           gmpy2.mpc_version()))
