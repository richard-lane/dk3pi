"""
Locations of analysis productions

"""


def get(magnetisation: str, year: int, data_type: str):
    """
    Find an analysis production location as a (glob on lxplus) str.

    :param magnetisation: Magnet direction, either "MagUp" or "MagDown"
    :param year: Data-taking year, one of (2015, 2016, 2017, 2018) as an int
    :param data_type: Data type: one of ("RS", "WS", "phsp")

    :raises AssertionError: if one of the arguments doesn't match the expected values
    :raises KeyError: if no production is found
    :return: absolute location of analysis production on lxplus, as a glob (i.e. path ending in *.root)

    """
    allowed_years_range = range(2015, 2019)
    assert year in allowed_years_range
    assert magnetisation in {"MagUp", "MagDown"}
    assert data_type in {"RS", "WS", "phsp"}

    # By magnetisation
    mag_up = {
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123278/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123280/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123276/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123308/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123306/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123310/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123286/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123288/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123290/0000/*.root",
    }
    mag_down = {
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123316/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123314/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123296/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123318/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123320/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123322/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123274/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123272/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123270/0000/*.root",
    }

    # By year
    year_2015 = set()
    year_2016 = {
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123278/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123280/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123276/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123316/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123314/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123296/0000/*.root",
    }
    year_2017 = {
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123308/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123306/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123310/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123318/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123320/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123322/0000/*.root",
    }
    year_2018 = {
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123286/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123288/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123290/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123274/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123272/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123270/0000/*.root",
    }

    # By decay type
    phsp = {
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123278/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123316/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123308/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123318/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123286/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123274/0000/*.root",
    }
    rs = {
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123280/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123314/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123306/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123320/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123288/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123272/0000/*.root",
    }
    ws = {
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123276/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2016/CHARM_D02HHHH_DVNTUPLE.ROOT/00123296/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123310/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2017/CHARM_D02HHHH_DVNTUPLE.ROOT/00123322/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123290/0000/*.root",
        "/eos/lhcb/grid/prod/lhcb/MC/2018/CHARM_D02HHHH_DVNTUPLE.ROOT/00123270/0000/*.root",
    }

    # sort of a hack with locals()
    mag_set = mag_down if magnetisation == "MagDown" else mag_up
    year_set = locals()[f"year_{year}"]
    type_set = locals()[data_type.lower()]

    # Find the intersection of our year/type/magnetisation sets
    # Hopefully, this is either a set containing 1 path (mainline case), or the empty set (nothing found)
    intersection = mag_set & year_set & type_set

    assert len(intersection) < 2

    return intersection.pop()
