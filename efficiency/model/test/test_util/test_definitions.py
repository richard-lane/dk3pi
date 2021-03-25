from util import definitions
from os.path import exists


def test_libs_built():
    """
    We have imported the definitions module: the CF and DCS C wrapper libraries should be built.

    """
    assert exists(definitions.CF_LIB)
    assert exists(definitions.DCS_LIB)

