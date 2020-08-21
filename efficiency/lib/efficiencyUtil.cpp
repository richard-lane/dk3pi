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
    return std::vector<double>{
        invariantMass(std::vector<kinematicParams_t>{event.kParams, event.pi1Params}),
        invariantMass(std::vector<kinematicParams_t>{event.pi1Params, event.pi2Params}),
        invariantMass(std::vector<kinematicParams_t>{event.pi2Params, event.pi3Params}),
        invariantMass(std::vector<kinematicParams_t>{event.kParams, event.pi1Params, event.pi2Params}),
        invariantMass(std::vector<kinematicParams_t>{event.pi1Params, event.pi2Params, event.pi3Params})};
}

double pT(const std::vector<kinematicParams_t>& particles)
{
    // Add all the px and py's together for our particles and find the pT of this combination
    kinematicParams_t aggregate;
    for (const kinematicParams_t& particle : particles) {
        aggregate.px += particle.px;
        aggregate.py += particle.py;
    }
    return pT(aggregate);
}

double pT(const dDecay_t& dDecay)
{
    const std::vector<kinematicParams_t> particles{
        dDecay.kParams, dDecay.pi1Params, dDecay.pi2Params, dDecay.pi3Params};
    return pT(particles);
}
