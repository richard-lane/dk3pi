/*
 * Functions for calculating the DCS and CF amplitude components of a D->K3pi event
 *
 */

#ifndef AMPLITUDES_H
#define AMPLITUDES_H

#include <complex>
#include <dlfcn.h>
#include <exception>
#include <iostream>
#include <string>

#include "efficiencyUtil.h"

/*
 * Couldn't find the specified function in library
 */
struct FunctionNotFoundException : public std::exception {
    FunctionNotFoundException(const std::string& name, const std::string& library) : name(name), library(library) { ; }

    const char* what() const throw()
    {
        std::cerr << dlerror() << std::endl;
        return ("Failed to find function " + name + " in library " + library).c_str();
    }

    const std::string name;
    const std::string library;
};

/*
 * Define a type for a function that will return an amplitude at a given point in phase space
 */
typedef std::complex<double> (*amplitudeFcnPtr)(double const* decayKinematics, const int& kaonCharge);

/*
 * Open a shared library and return a function from it by name
 *
 * The default function to look for in the library is named AMP, presumably short for amplitude
 */
amplitudeFcnPtr readFromSharedLib(const std::string& library, const std::string& name = "AMP");

/*
 * Returns the amplitude associated with a decay according to the model provided in library
 *
 * Needs to know about the kinematic parameters of the decay and the location of a fcn for calculating amplitudes
 *
 * NB implementation assumes that true == 1 (which should always be true)
 */
std::complex<double> amplitude(const dDecay_t& dDecayParameters, const amplitudeFcnPtr amplitudeFcn);

#endif // AMPLITUDES_H
