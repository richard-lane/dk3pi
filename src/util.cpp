/*
 * Utility functions that will be useful in multiple places.
 */
#ifndef UTIL_CPP
#define UTIL_CPP

#include <boost/filesystem.hpp>

#include "TCanvas.h"

#include "../include/util.h"

namespace util
{

/*
 * Combine a directory, filename and extension into a single boost:filesystem::path object.
 *
 * Effectively just does plotDir + plotName + fileExtension in an OS-intelligent way.
 *
 */
boost::filesystem::path concatPaths(std::string plotDir, std::string plotName, std::string fileExtension)
{
    boost::filesystem::path dir(plotDir);
    boost::filesystem::path file(plotName + fileExtension);
    return dir / file;
}

/*
 * Save a TObject to file
 * The extension in the specified path determines the format of the file.
 */
void saveToFile(TObject *myObject, const std::string &path, const std::string &drawOptions)
{
    TCanvas *c = new TCanvas();
    c->cd();
    myObject->Draw(drawOptions.c_str());
    c->SaveAs(path.c_str());

    delete c;
}

} // namespace util

#endif // UTIL_CPP
