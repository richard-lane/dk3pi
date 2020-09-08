#include <dlfcn.h>
#include <vector>

#include "amplitudes.h"

amplitudeFcnPtr readFromSharedLib(const std::string& library, const std::string& name)
{
    void*           libHandle{nullptr};
    amplitudeFcnPtr fcn = nullptr;

    // Open the shared library
    libHandle = dlopen(library.c_str(),
                       RTLD_NOW); // Special magic flag that tells us to resolve all symbols in the library before
                                  // returning, else throw an error

    // Set our function pointer and check that it was found properly
    fcn = (amplitudeFcnPtr)dlsym(libHandle, name.c_str());
    if (!fcn) {
        throw FunctionNotFoundException(name, library);
    }

    return fcn;
}

std::complex<double> amplitude(const dDecay_t& dDecayParameters, const amplitudeFcnPtr amplitudeFcn)
{
    // Translate our dDecayParameters into the right format for the amplitude function
    std::vector<double> eventData{
        dDecayParameters.kParams.px,
        dDecayParameters.kParams.py,
        dDecayParameters.kParams.pz,
        dDecayParameters.kParams.energy,

        dDecayParameters.pi1Params.px,
        dDecayParameters.pi1Params.py,
        dDecayParameters.pi1Params.pz,
        dDecayParameters.pi1Params.energy,

        dDecayParameters.pi2Params.px,
        dDecayParameters.pi2Params.py,
        dDecayParameters.pi2Params.pz,
        dDecayParameters.pi2Params.energy,

        dDecayParameters.pi3Params.px,
        dDecayParameters.pi3Params.py,
        dDecayParameters.pi3Params.pz,
        dDecayParameters.pi3Params.energy,
    };

    // Find the charge of the kaon
    // Either returns (2 * 1 - 1) = +1 or (2 * 0 - 1) = -1
    short int kCharge = 2 * dDecayParameters.kPlus - 1;

    return (*amplitudeFcn)(eventData.data(), kCharge);
}
