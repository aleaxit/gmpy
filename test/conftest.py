import pytest

import gmpy2


collect_ignore_glob = ['*.txt']


@pytest.fixture(autouse=True, scope='module')
def _set_default_context():
    gmpy2.set_context(gmpy2.context())
