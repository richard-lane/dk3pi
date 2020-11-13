/*
 * Add a DELTA_M branch to a prompt ROOT file
 *
 * Creates a copy of the original ROOT file ./copy_<original file name>
 */

#include <experimental/filesystem>
#include <iostream>
#include <regex>
#include <string>

#include <TFile.h>
#include <TTree.h>

#include "D2K3PiError.h"

/*
 * Extract name of a ROOT file from its path
 *
 * Regex matchs a filename ending in .root, which may look like e.g, "file.root" "./file.root",
 * "directory/file.root", etc.
 * Should definitely be unit tested, but it isn't because i don't want to
 */
static std::string rootFileName(const std::string& path)
{
    std::smatch matches;
    std::regex_search(path, matches, std::regex(R"(.*?\/?([\w\-. ]+\.root))"));
    if (matches.empty()) {
        std::cerr << "No root file found at " << path << std::endl;
        throw D2K3PiException();
    }

    return matches[1].str();
}

/*
 * If no branch containing the D mass difference exists, we may have to create one
 *
 * This does that
 *
 * Copies the old root file to copy_<filename>, in case of things going bad
 */
[[maybe_unused]] static void addDeltaMBranch(const std::string& path,
                                             const std::string& treeName,
                                             const std::string& dStarMassBranchName,
                                             const std::string& dMassBranchName,
                                             const std::string& deltaMBranchName)
{
    // Find the filename from the path so we know what to copy to
    std::string newPath{"copy_" + rootFileName(path)};

    // Copy old file before we do anything else, in case we break it
    std::cout << "Copying " << path << " to " << newPath << std::endl;
    std::experimental::filesystem::copy(path, newPath);

    std::cout << "Reading Tree" << std::endl;
    TFile  f(path.c_str(), "update");
    TTree* tree = f.Get<TTree>(treeName.c_str());
    assert(tree);

    std::cout << "Reading branches" << std::endl;
    double deltaM{0};
    double dMass{0};
    double dStarMass{0};
    // ROOT uses C-style return codes
    // a code of 0 or 1 are both ok
    int ok{0};
    ok = tree->SetBranchAddress(dMassBranchName.c_str(), &dMass);
    assert(ok == 0 || ok == 1);
    ok = tree->SetBranchAddress(dStarMassBranchName.c_str(), &dStarMass);
    assert(ok == 0 || ok == 1);

    // Populate the branch
    std::cout << "Creating branch " << deltaMBranchName << std::endl;
    TBranch* deltaMBranch{tree->Branch(deltaMBranchName.c_str(), &deltaM, (deltaMBranchName + "/D").c_str())};
    auto     numEvents = tree->GetEntries();
    for (decltype(numEvents) i = 0; i < numEvents; ++i) {
        tree->GetEntry(i);
        deltaM = dStarMass - dMass;
        deltaMBranch->Fill();
    }

    std::cout << "Writing tree" << std::endl;
    tree->Write("", TObject::kOverwrite);
}

int main(const int argc, const char* argv[])
{
    // Need to be passed file name, tree name, and the names of the branches
    if (argc != 6) {
        std::cerr
            << "Usage: add_delta_m_branch.exe <filepath> <treeName> <D* mass branch> <D mass branch> <delta M branch>"
            << std::endl;
        throw D2K3PiException();
    }

    addDeltaMBranch(argv[1], argv[2], argv[3], argv[4], argv[5]);
}
