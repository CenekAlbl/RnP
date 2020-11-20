#pragma once

#include <Eigen/Eigen>
#include <vector>
#include "utils.h"

struct RSCameraPose{
    Eigen::Matrix3d R; // Camera orientation as rotation matrix
    Eigen::Vector3d C; // Camera center
    Eigen::Vector3d w; // Camera rotational velocity (rotation axis where ||w|| is the angular velocity)
    Eigen::Vector3d t; // Camera translational velocity
    double f; // focal length (only R7Pf and R7Pfr compute focal lenght, otherwise 1)
    double rd; // radial distortion coefficient (only R7Pfr computes it, otherwise 0)
};

struct RSDoublelinCameraPose{
    Eigen::Vector3d v; // Camera orientation as I + [v]_x
    Eigen::Vector3d C; // Camera center
    Eigen::Vector3d w; // Camera rotational velocity (rotation axis where ||w|| is the angular velocity)
    Eigen::Vector3d t; // Camera translational velocity
    double f; // focal length (only R7Pf and R7Pfr compute focal lenght, otherwise 1)
    double rd; // radial distortion coefficient (only R7Pfr computes it, otherwise 0)
};

struct RSSinglelinCameraPose{
    Eigen::Vector3d v; // Camera orientation in Cayley parameterization
    Eigen::Vector3d C; // Camera center
    Eigen::Vector3d w; // Camera rotational velocity (rotation axis where ||w|| is the angular velocity)
    Eigen::Vector3d t; // Camera translational velocity
    double f; // focal length (only R7Pf and R7Pfr compute focal lenght, otherwise 1)
    double rd; // radial distortion coefficient (only R7Pfr computes it, otherwise 0)
};

typedef std::vector<RSCameraPose> RSCameraPoseVector;
typedef std::vector<RSDoublelinCameraPose> RSDoublelinCameraPoseVector;
typedef std::vector<RSSinglelinCameraPose> RSSinglelinCameraPoseVector;

int R7PfLin(Eigen::Matrix<double,7,3> X, Eigen::Matrix<double,7,2> u, Eigen::Vector3d vk, double r0, RSDoublelinCameraPoseVector * results);

int R7PfrLin(Eigen::Matrix<double,7,3> X, Eigen::Matrix<double,7,2> u, Eigen::Vector3d vk, double r0, RSDoublelinCameraPoseVector * results);

template<typename Model, int (*Solver)(Eigen::Matrix<double,7,3> X, Eigen::Matrix<double,7,2> u, Eigen::Vector3d vk, double r0, std::vector<Model> * results)>
int iterativeRnP(Eigen::Matrix<double,7,3> X, Eigen::Matrix<double,7,2> u, Eigen::Vector3d vk, int sampleSize, double r0, int maxIter, Model & result);

int r6pSingleLin(double * X, double * u, int direction, double r0, int maxpow, RSSinglelinCameraPoseVector * output);

int r6pDoubleLin(double * X, double * u, int direction, double r0, RSDoublelinCameraPoseVector * output);

double calcErrAlgebraicRnPFocalRadialDoubleLin(Eigen::Vector3d vr, Eigen::Vector3d Cr, Eigen::Vector3d wr, Eigen::Vector3d tr, double f, double rd, double r0, Eigen::MatrixXd X, Eigen::MatrixXd u);

