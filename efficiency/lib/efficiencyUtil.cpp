#include "efficiencyUtil.h"

double invariantMass(const std::vector<kinematicParams_t>& systemKinematics)
{
    double energy{0};
    double px{0};
    double py{0};
    double pz{0};
    for (auto it = systemKinematics.begin(); it != systemKinematics.end(); ++it) {
        energy += it->energy;
        px += it->px;
        py += it->py;
        pz += it->pz;
    }

    double pSquared    = px * px + py * py + pz * pz;
    double massSquared = energy * energy - pSquared;

    return std::sqrt(massSquared);
}

std::vector<double> event2invariantMasses(const dDecay_t& event)
{
    (void)event;
    return std::vector<double>{1, 2, 3};
}
