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

    py::class_<RSSinglelinCameraPose>(m, "RSSinglelinCameraPose")
        .def(py::init<>())
        .def_readwrite("v", &RSSinglelinCameraPose::v)
        .def_readwrite("C", &RSSinglelinCameraPose::C)
        .def_readwrite("w", &RSSinglelinCameraPose::w)
        .def_readwrite("t", &RSSinglelinCameraPose::t)
        .def_readwrite("f", &RSSinglelinCameraPose::f)
        .def_readwrite("rd", &RSSinglelinCameraPose::rd);

    m.def(
        "R7Pf",
        [](Eigen::MatrixXd X,
            Eigen::MatrixXd u,
            const Eigen::Vector3d &vk,
            double r0,
            int direction,
            int maxIter
            ) {
            RSDoublelinCameraPose model;
            int res =  iterativeRnP<RSDoublelinCameraPose, R7PfIter>(X, u, vk, 7, r0, direction, maxIter, model);
            return model;
            },
        "R7Pf - RS absolute pose with unknown focal length from 7 correspondences, both rotations linearized, R needs to be close to I");
    
    m.def(
        "R7Pfr",
        [](Eigen::MatrixXd X,
            Eigen::MatrixXd u,  
            const Eigen::Vector3d &vk,
            double r0,
            int direction,
            int maxIter
            ) {
            RSDoublelinCameraPose model;
            int res =  iterativeRnP<RSDoublelinCameraPose, R7PfrIter>(X, u, vk, 7, r0, direction, maxIter, model);
            return model;
            },
        "R7Pfr - RS absolute pose with unknown focal length and radial distortion from 7 correspondences, both rotations linearized, R needs to be close to I");

    m.def(
        "R6P2Lin",
        [](const Eigen::MatrixXd &X,
            const Eigen::MatrixXd &u, 
            double r0,
            int direction
            ) {
            RSDoublelinCameraPoseVector * results;
            int res =  R6P2Lin(X, u, direction, r0, results);
            return results;
            },
        "R6P double linearized - calibrated camera RS absolute pose from 6 correspondences, both rotations linearized, R needs to be close to I");

    m.def(
        "R6P1Lin",
        [](const Eigen::MatrixXd &X,
            const Eigen::MatrixXd &u, 
            double r0,
            int direction
            ) {
            RSSinglelinCameraPoseVector * results;
            int res =  R6P1Lin(X, u, direction, r0, 2, results);
            return results;
            },
        "R6P single linearized - calibrated camera RS absolute pose from 6 correspondences, only rotational velocity linearized, direct solution, no initialization needed");

    m.def(
        "R6PIter",
        [](Eigen::MatrixXd &X,
            Eigen::MatrixXd &u, 
            const Eigen::Vector3d &vk,
            double r0,
            int direction,
            int maxIter
            ) {
            RSDoublelinCameraPose model;
            int res =  iterativeRnP<RSDoublelinCameraPose, R6PIter>(X, u, vk, 6, r0, direction, maxIter, model);
            return model;
            },
        "R6P single linearized - calibrated camera RS absolute pose from 6 correspondences, only rotational velocity linearized, direct solution, no initialization needed");


}
