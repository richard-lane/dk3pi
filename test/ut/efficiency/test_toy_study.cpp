#ifndef TEST_TOY_STUDY_CPP
#define TEST_TOY_STUDY_CPP

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <TFile.h>
#include <boost/filesystem.hpp>

#include <random>

#include "efficiencyUtil.h"
#include "toyStudy.h"

/*
 * Test that the correct error gets thrown when our detection probability is less than 0 or greater than 1
 */
BOOST_AUTO_TEST_CASE(test_bad_detection_prob)
{
    // These return probabilities outside the range 0,1 so an error should be thrown if they are called when working out
    // if an event is detected
    auto minus2 = [](const dDecay_t& event) {
        (void)event;
        return -2.0;
    };
    auto plus2 = [](const dDecay_t& event) {
        (void)event;
        return 2.0;
    };

    // This returns a probability within the range 0,1 so shouldn't cause an error when we use it to work out if an
    // event is detected
    auto half = [](const dDecay_t& event) {
        (void)event;
        return 0.5;
    };

    std::mt19937          gen(1);
    std::vector<dDecay_t> events(1);

    BOOST_CHECK_THROW(applyEfficiency(&gen, minus2, events), EventDetectionProbNotNormalised);
    BOOST_CHECK_THROW(applyEfficiency(&gen, plus2, events), EventDetectionProbNotNormalised);
    BOOST_CHECK_NO_THROW(applyEfficiency(&gen, half, events));
}

/*
 * Test that we get the correct after applying the detection filter thing
 */
BOOST_AUTO_TEST_CASE(test_toy_efficiency)
{
    // Returns 1 if k_px = 1, 0 otherwise
    // 1 is exactly representable in binary (i hope??), so float equality is ok here
    auto detected = [](const dDecay_t& event) { return (event.kParams.px == 1) / (double)true; };

    std::mt19937 gen(1);
    dDecay_t     detectedEvent{.dParams   = {},
                           .kParams   = kinematicParams{.px = 1, .py = 0, .pz = 0, .energy = 0},
                           .pi1Params = {},
                           .pi2Params = {},
                           .pi3Params = {}};

    dDecay_t notDetectedEvent{.dParams   = {},
                              .kParams   = kinematicParams{.px = 2, .py = 0, .pz = 0, .energy = 0},
                              .pi1Params = {},
                              .pi2Params = {},
                              .pi3Params = {}};

    std::vector<dDecay_t> events{notDetectedEvent,
                                 detectedEvent,
                                 detectedEvent,
                                 notDetectedEvent,
                                 notDetectedEvent,
                                 detectedEvent,
                                 notDetectedEvent};
    std::vector<dDecay_t> expectedEvents{detectedEvent, detectedEvent, detectedEvent};

    applyEfficiency(&gen, detected, events);
    BOOST_CHECK(events == expectedEvents);
}

#endif // TEST_TOY_STUDY_CPP
