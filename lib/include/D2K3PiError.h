/*
 * Analysis-specific exceptions
 */

#ifndef D2K3PIERROR_H
#define D2K3PIERROR_H

#include <exception>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnPrint.h"

struct D2K3PiException : public std::exception {
    const char* what() const throw() { return "Base D2K3Pi Exception"; }
};

struct BadFitException : public D2K3PiException {

    BadFitException(ROOT::Minuit2::FunctionMinimum min) : min(min) { ; }

    const char* what() const throw()
    {
        std::cerr << min << std::endl;
        return "Fit failed";
    }

    ROOT::Minuit2::FunctionMinimum min;
};

#endif // D2K3PIERROR_h
