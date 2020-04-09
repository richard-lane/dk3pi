#include "BaseFitter.h"

BaseFitter::BaseFitter(const FitData_t& fitData)
{
    // Set our attributes to the right things
    // No need to perform consistency checks as they are performed by the FitData constructor
    _fitData = fitData;
}
