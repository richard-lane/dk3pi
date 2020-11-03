#include <array>
#include <iostream>

#include <boost/filesystem/path.hpp>

#include "bdt_reweighting.h"
#include <numpy/arrayobject.h>

/*
 * Take a vector of PhspPoints (i.e. a vector of vectors) and create a numpy reference to it
 *
 * Does copy the data, which could be an issue
 */
static PyArrayObject* copyPhspPoints2Numpy(const std::vector<PhspPoint>& points)
{
    // Assuming that all our points are the same dimensionality. which is reasonable
    // These have to be numpy's intp type
    npy_intp dimensions[2]{static_cast<npy_intp>(points.size()), static_cast<npy_intp>(points[0].size())};

    // Create an empty Python array
    PyArrayObject* array    = (PyArrayObject*)PyArray_SimpleNew(2, dimensions, NPY_DOUBLE);
    double*        arrayPtr = (double*)PyArray_DATA(array);

    // Copy data in line by line
    // Numpy arrays are 1d
    for (npy_intp iRow = 0; iRow < dimensions[0]; ++iRow) {
        std::copy(points[iRow].begin(), points[iRow].end(), arrayPtr + iRow * dimensions[1]);
    }

    // Are we returning a local variable by reference here..?
    return array;
}

static PyObject* castVector2Array(std::vector<double>& points)
{
    npy_intp size = points.size();
    return PyArray_SimpleNewFromData(1, &size, NPY_DOUBLE, points.data());
}

PyObject* initBDT(const std::vector<PhspPoint>& mcData,
                  const std::vector<PhspPoint>& realData,
                  std::vector<double>*          mcWeights,
                  std::vector<double>*          realWeights)
{
    // Initialise Python interpreter
    Py_Initialize();

    // Initialise Numpy array API or something
    if (!PyArray_API) {
        import_array();
    }

    // Weird hack telling this script where the Python file to import is
    boost::filesystem::path currentDir = boost::filesystem::path(__FILE__).parent_path();
    PyRun_SimpleString(("import sys; sys.path.append(\"" + currentDir.string() + "/../\")").c_str());

    // Import the reweighting module
    const std::string moduleName{"reweighting"};
    PyObject*         name   = PyUnicode_DecodeFSDefault(moduleName.c_str());
    PyObject*         module = PyImport_Import(name);
    Py_DECREF(name); // No longer need to keep track of the module name now that we have imported it
    if (!module) {
        throw PythonImportError(moduleName);
    }

    // Import the BDT initialisation function, check if it is callable
    const std::string fcnName{"init"};
    PyObject*         initFunc = PyObject_GetAttrString(module, fcnName.c_str());
    if (!initFunc || !PyCallable_Check(initFunc)) {
        throw PythonImportError(fcnName.c_str());
    }

    // Convert data to numpy
    PyArrayObject* numpyMCData   = copyPhspPoints2Numpy(mcData);
    PyArrayObject* numpyRealData = copyPhspPoints2Numpy(realData);
    assert(numpyMCData);
    assert(numpyRealData);

    // Create + return the BDT
    // Four different cases requiring different types to be passed to PyObject_CallFunctionObjArgs;
    // Both/neither/one of mc and real weights specified
    // Kind of a mess
    PyObject* bdt = nullptr;

    // Both weights
    if (mcWeights && realWeights) {
        PyObject* mcWeightArray   = castVector2Array(*mcWeights);
        PyObject* realWeightArray = castVector2Array(*realWeights);
        bdt = PyObject_CallFunctionObjArgs(initFunc, numpyMCData, numpyRealData, mcWeightArray, realWeightArray, NULL);
    }

    // Only MC weights
    else if (mcWeights) {
        PyObject* mcWeightArray = castVector2Array(*mcWeights);
        bdt = PyObject_CallFunctionObjArgs(initFunc, numpyMCData, numpyRealData, mcWeightArray, NULL);
    }

    // Only real weights
    else if (realWeights) {
        PyObject* mcWeightArray   = Py_BuildValue(""); // python None
        PyObject* realWeightArray = castVector2Array(*realWeights);
        bdt = PyObject_CallFunctionObjArgs(initFunc, numpyMCData, numpyRealData, mcWeightArray, realWeightArray, NULL);
    }

    // Neither weights
    else {
        bdt = PyObject_CallFunctionObjArgs(initFunc, numpyMCData, numpyRealData, NULL);
    }

    if (!bdt) {
        throw PythonError();
    }
    return bdt;

    // Might have some issues with reference counts not being 0 when i return
    // But maybe it's fine
}

std::vector<double> efficiency(PyObject* bdt, const std::vector<PhspPoint>& points, const size_t expectedNumPoints)
{
    assert(bdt);
    // Find the Python function we'll be calling
    const std::string moduleName{"reweighting"};
    const std::string fcnName{"predicted_weights"};
    PyObject*         name   = PyUnicode_DecodeFSDefault(moduleName.c_str());
    PyObject*         module = PyImport_Import(name);
    PyObject*         func   = PyObject_GetAttrString(module, fcnName.c_str());
    assert(func);
    Py_DECREF(name);
    Py_DECREF(module);

    // Convert the number of points we expect
    PyObject* numPoints = PyLong_FromSize_t(expectedNumPoints);
    assert(numPoints);

    // Convert the points to numpy
    PyArrayObject* numpyData = copyPhspPoints2Numpy(points);
    assert(numpyData);

    // Find a numpy array of efficiency values
    PyObject* efficiencies = PyObject_CallFunctionObjArgs(func, bdt, numpyData, numPoints, NULL);
    assert(efficiencies);
    Py_DECREF(func);
    Py_DECREF(numpyData);
    Py_DECREF(numPoints);

    // Convert numpy array of efficiencies to C++ vector
    std::vector<double> data(points.size(), std::numeric_limits<double>::quiet_NaN());
    PyArrayIterObject*  iter = (PyArrayIterObject*)PyArray_IterNew(efficiencies);
    if (!iter) {
        throw PythonError();
    }
    double* dptr = static_cast<double*>(PyArray_DATA((PyArrayObject*)efficiencies));
    for (size_t i = 0; i < points.size(); ++i) {
        data[i] = *dptr;
        dptr++;
    }
    return data;
}
