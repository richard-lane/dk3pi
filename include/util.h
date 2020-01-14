#ifndef UTIL_H
#define UTIL_H

#include <boost/filesystem.hpp>

#include "TObject.h"

namespace util
{

boost::filesystem::path concatPaths(std::string plotDir, std::string plotName, std::string fileExtension);

void saveToFile(TObject *myObject, const std::string &path, const std::string &drawOptions="");

} // namespace util

#endif // UTIL_H
