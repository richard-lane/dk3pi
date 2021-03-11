#include "wrapper_utils.h"

#include "dcs.cpp"

/*
 * Return the DCS amplitude of an event as an AmpGenWrapper::Complex_t
 *
 * event:      array of 16 doubles (K_px, K_py, K_py, K_E, ...) etc. for K, pi1, pi2, pi3
 * kaonCharge: charge of the kaon, either +1 or -1
 *
 */
extern "C" AmpGenWrapper::Complex_t dcs_wrapper(double event[16], const int kCharge)
{
    // The definition of AMP is provided in dcs.cpp
    return AmpGenWrapper::wrapper(event, kCharge, AMP);
}
