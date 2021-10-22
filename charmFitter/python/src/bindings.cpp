#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "bindings.h"
#include "simulator.h"
#include "cleoFitter.h"
#include "ConstrainedFitter.h"

static void initSimulator(pybind11::module_& m) {
    m.def("simulate", &simulate, "decay simulator wrapper");
    m.def("expectedParams", &expectedParamsBinding, "expected abc params");
}

static void initCleoScan (pybind11::module_& m) {
    m.def("cleoZScan", &cleoZScan, "cleo scan wrapper");
}

static void initCleoFitter (pybind11::module_& m) {
    pybind11::class_<CharmFitter::CLEOCombinationFitter>(m, "CLEOCombinationFitter")
        .def(pybind11::init<const std::vector<double>&,
                            const std::array<double, 6>&,
                            const std::array<double, 6>&,
                            const int>())
        .def("addRSPoints", &CharmFitter::CLEOCombinationFitter::addRSPoints<std::vector<double>>)
        .def("addWSPoints", &CharmFitter::CLEOCombinationFitter::addWSPoints<std::vector<double>>)
        .def("addWSPoints", &CharmFitter::CLEOCombinationFitter::addWSPoints<std::vector<double>>)
        .def("fixParameter", &CharmFitter::CLEOCombinationFitter::fixParameter)
        .def("freeParameter", &CharmFitter::CLEOCombinationFitter::freeParameter)
        .def("fixParameters", &CharmFitter::CLEOCombinationFitter::fixParameters<std::vector<std::string>>)
        .def("freeParameters", &CharmFitter::CLEOCombinationFitter::freeParameters<std::vector<std::string>>)
        .def("setParameter", &CharmFitter::CLEOCombinationFitter::setParameter)
        .def("getBinCentres", &CharmFitter::CLEOCombinationFitter::getBinCentres)
        .def("getBinWidths", &CharmFitter::CLEOCombinationFitter::getBinWidths)
        .def("getRSBinContent", &CharmFitter::CLEOCombinationFitter::getRSBinContent)
        .def("getWSBinContent", &CharmFitter::CLEOCombinationFitter::getWSBinContent)
        .def("ratios", &CharmFitter::CLEOCombinationFitter::ratios)
        .def("errors", &CharmFitter::CLEOCombinationFitter::errors)
        .def("fit", &CharmFitter::CLEOCombinationFitter::fit);

}

static void initFitResults (pybind11::module_& m) {
    pybind11::class_<CharmFitter::FitResults_t>(m, "FitResults_t")
        .def_readwrite("fitStatistic", &CharmFitter::FitResults_t::fitStatistic)
        .def_readwrite("fitParams", &CharmFitter::FitResults_t::fitParams)
        .def_readwrite("fitParamErrors", &CharmFitter::FitResults_t::fitParamErrors);
}

static void initConstants (pybind11::module_& m) {
    m.attr("WORLD_AVERAGE_X")     = pybind11::float_(CharmFitter::WORLD_AVERAGE_X);
    m.attr("WORLD_AVERAGE_X_ERR") = pybind11::float_(CharmFitter::WORLD_AVERAGE_X_ERR);
    m.attr("WORLD_AVERAGE_Y")     = pybind11::float_(CharmFitter::WORLD_AVERAGE_Y);
    m.attr("WORLD_AVERAGE_Y_ERR") = pybind11::float_(CharmFitter::WORLD_AVERAGE_Y_ERR);
    m.attr("X_Y_CORRELATION")     = pybind11::float_(CharmFitter::X_Y_CORRELATION);
}

PYBIND11_MODULE(libcleoScan, m) {
    m.doc() = "Scan with CLEO fitter; also other stuff";

    initSimulator(m);
    initCleoScan(m);
    initCleoFitter(m);
    initFitResults(m);
    initConstants(m);
}

