#include <boost/test/unit_test.hpp>

#include <TFile.h>
#include <boost/filesystem.hpp>

#include "ReadRoot.h"

/*
 * Test that the correct error gets thrown when the wrong number of branches is specified when trying to read a ROOT
 * file
 */
BOOST_AUTO_TEST_CASE(test_wrong_number_branches)
{
    boost::filesystem::path currentDir = boost::filesystem::path(__FILE__).parent_path();
    boost::filesystem::path rootFile("../../it/efficiency/dBarCf.root");
    std::unique_ptr<TFile>  cfFile(new TFile((currentDir / rootFile).string().c_str()));

    std::vector<std::string> threeStrs = {"a", "b", "c"};
    std::vector<std::string> fiveStrs  = {"a", "b", "c", "d", "e"};

    BOOST_CHECK_THROW(ReadRoot(cfFile.get(), "DalitzEventList", threeStrs), InvalidBranchesException);
    BOOST_CHECK_THROW(ReadRoot(cfFile.get(), "DalitzEventList", fiveStrs), InvalidBranchesException);
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

/*
 * Test we get a filenotfound error when the root file doesnt exist
 */
BOOST_AUTO_TEST_CASE(test_rootfile_not_found)
{
    boost::filesystem::path currentDir = boost::filesystem::path(__FILE__).parent_path();
    boost::filesystem::path rootFile("../../it/efficiency/dBarCf.root");
    boost::filesystem::path nonExistentRootFile("bleh.root");

    // Unit test probably shouldn't allocate memory but meh
    std::unique_ptr<TFile> cfFile(new TFile((currentDir / rootFile).string().c_str()));
    std::unique_ptr<TFile> notFile(new TFile(nonExistentRootFile.string().c_str()));

    std::vector<std::string> realBranchNames = {"_1_K~", "_2_pi#", "_3_pi#", "_4_pi~"};
    BOOST_CHECK_NO_THROW(ReadRoot(cfFile.get(), "DalitzEventList", realBranchNames));
    BOOST_CHECK_THROW(ReadRoot(notFile.get(), "DalitzEventList", realBranchNames), InvalidRootFile);
}

/*
 * Test we get a treenotfound error when the tree doesnt exist
 */
BOOST_AUTO_TEST_CASE(test_tree_not_found)
{
    boost::filesystem::path currentDir = boost::filesystem::path(__FILE__).parent_path();
    boost::filesystem::path rootFile("../../it/efficiency/dBarCf.root");

    // Unit test probably shouldn't allocate memory but meh
    std::unique_ptr<TFile> cfFile(new TFile((currentDir / rootFile).string().c_str()));
    std::string            realTree{"DalitzEventList"};
    std::string            fakeTree{"boiboiboib"};

    std::vector<std::string> realBranchNames = {"_1_K~", "_2_pi#", "_3_pi#", "_4_pi~"};
    BOOST_CHECK_NO_THROW(ReadRoot(cfFile.get(), realTree, realBranchNames));
    BOOST_CHECK_THROW(ReadRoot(cfFile.get(), fakeTree, realBranchNames), InvalidTree);
}
