#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "rnp.h"

namespace py = pybind11;

PYBIND11_MODULE(pyrnp, m) {
    m.doc() = "python binding for the RnP library"; // optional module docstring

    py::class_<RSDoublelinCameraPose>(m, "RSDoublelinCameraPose")
        .def(py::init<>())
        .def_readwrite("v", &RSDoublelinCameraPose::v)
        .def_readwrite("C", &RSDoublelinCameraPose::C)
        .def_readwrite("w", &RSDoublelinCameraPose::w)
        .def_readwrite("t", &RSDoublelinCameraPose::t)
        .def_readwrite("f", &RSDoublelinCameraPose::f)
        .def_readwrite("rd", &RSDoublelinCameraPose::rd);

    m.def(
        "r7pf",
        [](const Eigen::Matrix<double,3,7> &X,
            const Eigen::Matrix<double,2,7> &u,
            const Eigen::Vector3d &vk,
            double r0,
            int maxIter
            ) {
            RSDoublelinCameraPose model;
            int res =  iterativeRnP<RSDoublelinCameraPose, R7PfLin>(X.transpose(), u.transpose(), vk, 7, r0, maxIter, model);
            return model;
        },
        "R7Pf - RS absolute pose with unknown focal length from 7 correspondences");

}
