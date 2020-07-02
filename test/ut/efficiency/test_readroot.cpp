#ifndef TEST_READROOT_CPP
#define TEST_READROOT_CPP

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <TFile.h>

#include "ReadRoot.h"

/*
 * Test that the correct error gets thrown when the wrong number of branches is specified when trying to read a ROOT
 * file
 */
BOOST_AUTO_TEST_CASE(test_wrong_number_branches)
{
    std::vector<std::string> threeStrs = {"a", "b", "c"};
    std::vector<std::string> fiveStrs  = {"a", "b", "c", "d", "e"};

    BOOST_CHECK_THROW(ReadRoot(nullptr, "foo", threeStrs), InvalidBranchesException);
    BOOST_CHECK_THROW(ReadRoot(nullptr, "foo", fiveStrs), InvalidBranchesException);
    // Don't check for valid branches as this is really slow; the constructor reads all the data from the ROOT file.
    // Maybe this means I should devolve responsibility more...
}

/*
 * Test that the correct error gets thrown when the specified branch doesn't exist on the ROOT file
 * file
 */
BOOST_AUTO_TEST_CASE(test_branch_not_found)
{
    boost::filesystem::path currentDir = boost::filesystem::path(__FILE__).parent_path();
    boost::filesystem::path rootFile("../../it/efficiency/dBarCf.root");

    // Unit test probably shouldn't allocate memory but meh
    TFile* cfFile = new TFile((currentDir / rootFile).string().c_str());

    std::vector<std::string> realBranchNames = {"_1_K~", "_2_pi#", "_3_pi#", "_4_pi~"};
    BOOST_CHECK_NO_THROW(ReadRoot(cfFile, "DalitzEventList", realBranchNames));

    std::vector<std::string> fakeBranchNames = {"greg", "bob", "fred", "laura"};
    BOOST_CHECK_THROW(ReadRoot(cfFile, "DalitzEventList", fakeBranchNames), BranchNotFoundException);

    delete cfFile;
}

#endif // TEST_READROOT_CPP
