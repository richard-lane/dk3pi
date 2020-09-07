#ifndef FIT_ERROR_H
#define FIT_ERROR_H

#include "D2K3PiError.h"

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>

struct BadFitException : public D2K3PiException {

    BadFitException(ROOT::Minuit2::FunctionMinimum min) : min(min) { ; }

    const char* what() const throw()
    {
        std::cerr << min << std::endl;
        return "Fit failed";
    }

    ROOT::Minuit2::FunctionMinimum min;
};

#endif // FIT_ERROR_H
