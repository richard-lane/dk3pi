#ifndef BDT_REWEIGHTING_H
#define BDT_REWEIGHTING_H

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <vector>

#include <Python.h>
#include <numpy/numpyconfig.h>

#include "D2K3PiError.h"
#include "efficiencyUtil.h"

/*
 * Generic python error
 */
struct PythonError : public D2K3PiException {
    const char* what() const throw()
    {
        PyErr_Print();
        return "Error with Python interface";
    }
};

/*
 * Error importing something from python
 */
struct PythonImportError : public D2K3PiException {
    PythonImportError(const std::string module) : module(module) {}
    const char* what() const throw() { return ("Error importing " + module).c_str(); }

  private:
    std::string module;
};

/*
 * Cast a PyObject to a char*, mostly for debugging and err messages
 *
 * Probably doesn't work with python strings
 */
inline char* py2chars(PyObject* obj)
{
    return PyBytes_AS_STRING(PyUnicode_AsEncodedString(PyObject_Repr(obj), "utf-8", "~E~"));
}

/*
 * Set up a hep_ml reweighting BDT using the MC and real data passed in
 *
 * Optionally pass in weights. i don't know why they're not const
 *
 */
PyObject* initBDT(const std::vector<PhspPoint>& mcData,
                  const std::vector<PhspPoint>& realData,
                  std::vector<double>*          mcWeights       = nullptr,
                  std::vector<double>*          realDataWeights = nullptr);

/*
 * Find the efficiency at a collection of points given a trained BDT
 *
 * weights           - initial weights
 * expectedNumPoints - the number of points expected in the distribution we're targeting
 */
std::vector<double> efficiency(PyObject*                     bdt,
                               const std::vector<PhspPoint>& points,
                               std::vector<double>*          weights,
                               const size_t                  expectedNumPoints);

#endif // BDT_REWEIGHTING_H
