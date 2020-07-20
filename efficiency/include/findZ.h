#ifndef FINDZ_H
#define FINDZ_H

#include <complex>
#include <vector>

#include "amplitudes.h"
#include "efficiencyUtil.h"

struct PhaseOutOfRangeException : public std::exception {
    PhaseOutOfRangeException(const double phase) : phase(phase) { ; }

    const char* what() const throw()
    {
        return ("Phase " + std::to_string(phase) + " must be in range -pi to +pi for calculating DCS fiddle factor")
            .c_str();
    }

    const double phase;
};

/*
 * NB our DCS/CF events are D0 -> K+3pi / Dbar0 -> K+3pi respectively
 *
 * In order to calculate Z, we need to apply a fiddle factor to the one of the amplitudes that AmpGen returns; this is
 * calculated here + should be applied to the DCS amplitude.
 *
 * Since AmpGen amplitudes are normalised to the same thing (presumably 1), but we want our DCS amplitude to be much
 * smaller than our DCS amplitude, we will need to scale them. We can do this by scaling the DCS amplitude by something
 * that depends on the number of DCS and CF events.
 *
 * We also want our DCS and CF amplitudes to have a relative strong phase of 0, so we'll need to know the average strong
 * phase between the two amplitudes before scaling. This is defined as the phase of Integral(AB*) over phase space,
 * where A is the CF amplitude and B* the complex conjuage of the DCS amplitude. In radians pls between -pi and pi
 */
std::complex<double>
amplitudeFiddleFactor(const size_t numDcsEvents, const size_t numCfEvents, const double strongPhaseInRadians);

/*
 * Find the charm interference parameter Z given a set of D -> K3pi decays
 *
 * Takes in an "offset" the DCS and CF amplitudes
 * This should be such that we get DCS/CF amplitude ratio ~ 0.055 and average relative strong phase ~ 0
 */
std::complex<double> charmInterferenceZ(const std::vector<dDecay_t>& dDecayKinematics,
                                        const amplitudeFcnPtr        dcsAmplitudeFcn,
                                        const amplitudeFcnPtr        cfAmplitudeFcn,
                                        const std::complex<double>&  dcsOffset = std::polar<double>(1, 0.0));

#endif // FINDZ_H
