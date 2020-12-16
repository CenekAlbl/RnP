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
        "R7Pf - RS absolute pose with unknown focal length from 7 correspondences, both rotations linearized, R needs to be close to I");
    
    m.def(
        "r7pfr",
        [](const Eigen::Matrix<double,3,7> &X,
            const Eigen::Matrix<double,2,7> &u,  
            const Eigen::Vector3d &vk,
            double r0,
            int maxIter
            ) {
            RSDoublelinCameraPose model;
            int res =  iterativeRnP<RSDoublelinCameraPose, R7PfrLin>(X.transpose(), u.transpose(), vk, 7, r0, maxIter, model);
            return model;
            },
        "R7Pfr - RS absolute pose with unknown focal length and radial distortion from 7 correspondences, both rotations linearized, R needs to be close to I");

    m.def(
        "r6pDoubleLin",
        [](const Eigen::Matrix<double,3,6> &X,
            const Eigen::Matrix<double,2,6> &u, 
            int direction,
            double r0
            ) {
            RSDoublelinCameraPoseVector * results;
            int res =  r6pDoubleLin(X, u, direction, r0, results);
            return results;
            },
        "R6P double linearized - calibrated camera RS absolute pose from 6 correspondences, both rotations linearized, R needs to be close to I");

    m.def(
        "r6pSingleLin",
        [](const Eigen::Matrix<double,3,6> &X,
            const Eigen::Matrix<double,2,6> &u, 
            int direction,
            double r0
            ) {
            RSSinglelinCameraPoseVector * results;
            int res =  r6pSingleLin(X, u, direction, r0, 2, results);
            return results;
            },
        "R6P single linearized - calibrated camera RS absolute pose from 6 correspondences, only rotational velocity linearized, direct solution, no initialization needed");


}
