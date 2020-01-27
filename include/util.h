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

} // namespace util

#endif // UTIL_H
