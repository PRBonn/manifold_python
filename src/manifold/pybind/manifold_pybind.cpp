// @file      manifold_pybind.cpp
// @author    Ignacio Vizzo     [ivizzo@uni-bonn.de]
//
// Copyright (c) 2020 Ignacio Vizzo, all rights reserved
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <Eigen/Core>
#include <vector>

#include "manifold/ManifoldProcessor.h"

PYBIND11_MAKE_OPAQUE(std::vector<Eigen::Vector3d>);
PYBIND11_MAKE_OPAQUE(std::vector<Eigen::Vector3i>);

namespace py = pybind11;
using namespace py::literals;

namespace manifold {

PYBIND11_MODULE(manifold_pybind, m) {
    py::bind_vector<std::vector<Eigen::Vector3d>>(m, "_VectorEigen3d");
    py::bind_vector<std::vector<Eigen::Vector3i>>(m, "_VectorEigen3i");
    py::class_<ManifoldProcessor>(
            m, "_Processor",
            "This is the low level C++ bindings, all the methods and "
            "constructor defined within this module (starting with a ``_`` "
            "should not be used. Please reffer to the python Procesor class to "
            "check how to use the API")
            .def(py::init<const std::vector<Eigen::Vector3d> &,
                          const std::vector<Eigen::Vector3i> &>())
            .def("_get_manifold_mesh", &ManifoldProcessor::GetManifoldMesh)
            .def_static("_is_manifold", &ManifoldProcessor::is_manifold);
}
}  // namespace manifold
