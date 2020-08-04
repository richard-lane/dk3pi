BRANCHES = [
    "Dstar_M",
    "D0_IPCHI2_OWNPV",
    "D0_PE",
    "D0_PX",
    "D0_PY",
    "D0_PZ",
    "D0_M",
    "D0_TAU",
    "Kminus_PE",
    "Kminus_PX",
    "Kminus_PY",
    "Kminus_PZ",
    "pi1plus_PE",
    "pi1plus_PX",
    "pi1plus_PY",
    "pi1plus_PZ",
    "pi2plus_PE",
    "pi2plus_PX",
    "pi2plus_PY",
    "pi2plus_PZ",
    "pi3minus_PE",
    "pi3minus_PX",
    "pi3minus_PY",
    "pi3minus_PZ",
]


# Cuts should start with "cut"
def cut_d_impact_param(D0_IPCHI2_OWNPV):
    return D0_IPCHI2_OWNPV < 9


# They can accept multiple parameters
def cut_delta_m(Dstar_M, D0_M):
    return 139.3 < Dstar_M - D0_M < 152


# Functions not matching the pattern cut* will be ignored
def not_a_cut_function(f):
    assert False


def d_mass_helper(D0_M):
    return 1840.83 < D0_M < 1888.83


# Functions can call other functions, in case your cut is complicated
def cut_d_mass(D0_M):
    return d_mass_helper(2)

