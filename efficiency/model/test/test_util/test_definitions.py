from util import definitions


def test_assign_cf():
    """
    Test assigning to the CF fcn handle

    """
    assert definitions.CF_FCN is None

    definitions.assign_cf("hello")
    assert definitions.CF_FCN == "hello"


def test_assign_cf_persistence():
    """
    Test that our CF fcn handle has maintained its state

    """
    assert definitions.CF_FCN == "hello"


def test_assign_dcs():
    """
    Test assigning to the DCS fcn handle

    """
    assert definitions.DCS_FCN is None

    definitions.assign_dcs("bye")
    assert definitions.DCS_FCN == "bye"


def test_assign_dcs_persistence():
    """
    Test that our DCS fcn handle has maintained its state

    """
    assert definitions.DCS_FCN == "bye"
