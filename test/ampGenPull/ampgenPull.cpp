#include "D2K3PiError.h"

/*
 * Read in two ROOT files d.root and dBar.root, take their ratios of decay times
 */
void ampgenFit(const char* dFile, const char* dBarFile)
{
    // Read in ROOT files
    

    // Take ratios

    // Perform a fit

    // Save the fit parameters and errors to a text file
}

int main(int argc, char* argv[])
{
    // Should probably also check that D and Dbar have been passed in the right order
    // but what are the chances of that
    if (argc != 4) {
        std::cerr << "Usage: ./ampgenpull <d root file> <d bar root file>" << std::endl;
        throw D2K3PiException();
    }

    ampgenFit(argv[1], argv[2]);
    return 0;
}
