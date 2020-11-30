#include "rnp.h"
#include <iostream>

double calcErrAlgebraicRnPFocalRadialDoubleLin(Eigen::Vector3d vr, Eigen::Vector3d Cr, Eigen::Vector3d wr, Eigen::Vector3d tr, double f, double rd, double r0, Eigen::MatrixXd X, Eigen::MatrixXd u, int direction){
    double err = 0;
    for (int i = 0; i < X.cols(); i++)
    {
        Eigen::Vector3d uh; 
        uh << u.col(i), 1 + rd*(u(0,i)*u(0,i) + u(1,i)*u(1,i));
        Eigen::Matrix3d K = Eigen::Matrix3d::Identity();
        K(2,2) = 1/f;
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        double rs = (u(direction,i)-r0);
        Eigen::Vector3d eq = X_(uh)*K*((I + rs*X_(wr))*(I + X_(vr))*X.col(i) + Cr + rs*tr);
        err += eq.cwiseAbs().sum();
    }
    return err;
}

double calcErrAlgebraicRnPFocalRadialSingleLin(Eigen::Vector3d vr, Eigen::Vector3d Cr, Eigen::Vector3d wr, Eigen::Vector3d tr, double f, double rd, double r0, Eigen::MatrixXd X, Eigen::MatrixXd u, int direction){
    double err = 0;
    for (int i = 0; i < X.cols(); i++)
    {
        Eigen::Vector3d uh; 
        uh << u.col(i), 1 + rd*(u(0,i)*u(0,i) + u(1,i)*u(1,i));
        Eigen::Matrix3d K = Eigen::Matrix3d::Identity();
        K(2,2) = 1/f;
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        double rs = (u(direction,i)-r0);
        Eigen::Matrix3d Rv;
        Rv = Eigen::AngleAxis<double>(vr.norm(), vr).toRotationMatrix();
        Eigen::Vector3d eq = X_(uh)*K*((I + rs*X_(wr))*Rv*X.col(i) + Cr + rs*tr);
        err += eq.cwiseAbs().sum();
    }
    return err;
}

