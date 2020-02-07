#ifndef UTIL_H
#define UTIL_H

#include <boost/filesystem.hpp>

#include "TObject.h"

namespace util
{

/*
 * Combine a directory, filename and extension into a single boost:filesystem::path object.
 *
 * Effectively just does plotDir + plotName + fileExtension in an OS-intelligent way.
 *
 */
boost::filesystem::path concatPaths(std::string plotDir, std::string plotName, std::string fileExtension);

/*
 * Save a TObject to file
 * The extension in the specified path determines the format of the file.
 */
void saveToFile(TObject *myObject, const std::string &path, const std::string &drawOptions = "");

/*
 * Utility function to find bin limits given a dataset and the minimum number of points in each bin
 * Creates as many bins as possible containing exactly the minimum number of points; the last bin may contain more
 * points than the minimum.
 *
 * Note this may result in some very wide bins at the start/end if lowBin/highBin are specified as much lower or higher
 * than the extreme values in the dataset.
 *
 * TODO: this currently sucks for the above reason. so make it not suck
 *
 * Params:
 *   dataset         - data to find bins for. Should be sorted.
 *   minPointsPerBin - minimum number of points per bin.
 *   lowBin          - lower edge of the first bin.
 *   highBin         - upper edge of the last bin.
 */
std::vector<double> findBinLimits(const std::vector<double> &dataSet,
                                  const size_t               minPointsPerBin,
                                  const double               lowBin,
                                  const double               highBin);

} // namespace util

#endif // UTIL_H
