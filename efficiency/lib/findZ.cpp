#include <cmath>

#include "amplitudes.h"
#include "efficiencyUtil.h"
#include "findZ.h"

std::complex<double>
amplitudeFiddleFactor(const size_t numDcsEvents, const size_t numCfEvents, const double strongPhaseInRadians)
{
    if (strongPhaseInRadians < -M_PI_2 || strongPhaseInRadians > M_PI_2) {
        throw PhaseOutOfRangeException(strongPhaseInRadians);
    }

    double magnitude = std::sqrt((double)numDcsEvents / (double)numCfEvents);
    return std::polar<double>(magnitude, strongPhaseInRadians);
}

std::complex<double> charmInterferenceZ(const std::vector<dDecay_t>& dDecayKinematics,
                                        const amplitudeFcnPtr        dcsAmplitudeFcn,
                                        const amplitudeFcnPtr        cfAmplitudeFcn,
                                        const std::complex<double>&  dcsOffset)
{
    // Initialises to 0 by default
    std::complex<double> z;
    double               cf{0.0};
    double               dcs{0.0};

    for (auto it = dDecayKinematics.begin(); it != dDecayKinematics.end(); ++it) {
        std::complex<double> cfAmplitude = amplitude(*it, cfAmplitudeFcn);
        std::complex<double> dcsAmplitude =
            amplitude(*it, dcsAmplitudeFcn) * dcsOffset; // Our DCS amplitude must be modulated by this offset

        cf += std::norm(cfAmplitude);
        dcs += std::norm(dcsAmplitude);
        z += cfAmplitude * std::conj(dcsAmplitude);
    }

    // Scale our interference parameter by the total CF and DCS parameters
    return z / (sqrt(cf * dcs));
}
