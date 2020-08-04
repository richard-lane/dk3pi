def cut_d_mass(D0_M):
    """
    Returns bool for our D mass cut
    True if d_masses[i] is in the allowed range

    d_masses in MeV

    """
    return 1840.83 < D0_M < 1888.83


def cut_delta_m(Dstar_M, D0_M):
    """
    Returns bool for our DeltaM cut
    True if deltam[i] is in the allowed range

    masses in MeV

    """
    return 139.3 < Dstar_M - D0_M < 152


def cut_d_impact_param(D0_IPCHI2_OWNPV):
    """
    Returns bool for impact parameter cut
    True if our chisq is below the allowed value

    Takes in D0_IPCHI2_OWNPV

    """
    return D0_IPCHI2_OWNPV < 9


def not_a_cut_function(f):
    assert False


def uses_helper():
    """
    Check if this works

    """
    not_a_cut_function(2)
