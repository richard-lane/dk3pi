#include <vector>

#ifndef TESTUTIL_H
#define TESTUTIL_H

/*
 * Find the expected ratio at a given time
 */
inline double ratio(const double a, const double b, const double c, const double time)
{
    return a + b * time + c * time * time;
}

/*
 * Create a vector of ratios with gaussian noise
 */
std::vector<double>
idealRatios(const std::vector<double>& times, const double error, const double a, const double b, const double c);

#endif // TESTUTIL_H
