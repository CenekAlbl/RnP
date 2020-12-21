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
        double n = vr.norm();
        if(n < 1e-15){
            Rv = Eigen::AngleAxis<double>(0, (Eigen::Vector3d() << 1,0,0).finished()).toRotationMatrix();
        }else{
            Rv = Eigen::AngleAxis<double>(n, vr/n).toRotationMatrix();
        }
        
        Eigen::Vector3d eq = X_(uh)*K*((I + rs*X_(wr))*Rv*X.col(i) + Cr + rs*tr);
        err += eq.cwiseAbs().sum();
    }
    return err;
}


int rsDoubleLinProjection(const Eigen::Vector3d &X, Eigen::Vector2d &u, const Eigen::Vector3d &v, const Eigen::Vector3d &C, const Eigen::Vector3d &w, const Eigen::Vector3d &t, double f, double rd, double r0, int direction){
    // First initiate u with a global shutter projection
    Eigen::Matrix<double,3,4> P;
    Eigen::Vector3d uh;
    Eigen::Matrix3d K = Eigen::Matrix3d::Zero();
    K.diagonal() << f, f, 1;
    Eigen::Matrix3d Rv = Eigen::Matrix3d::Identity() + X_(v); 
    uh = K* ( Rv * X + C );
    u = uh.head(2)/uh(2);


    double diff = 1e15;

    int niter = 0;

    Eigen::Vector3d temp;

    while(diff > 1e-10 ){
        double rc2 = u(0)*u(0) + u(1)*u(1);
        temp = K * ((Eigen::Matrix3d::Identity() + X_(u(direction) * w)) * Rv * X + C + u(direction) * t);
        temp = temp/temp(2);
        temp.head(2) *= (1 + rd * rc2);
        diff = (u - temp.head(2)).norm();
        u = temp.head(2);
        if(niter > 100){
            return WARN_NO_CONVERGENCE;
        }
        niter++;
    }

    return 0;
    
}

int rsSingleLinProjection(const Eigen::Vector3d &X, Eigen::Vector2d &u, const Eigen::Vector3d &v, const Eigen::Vector3d &C, const Eigen::Vector3d &w, const Eigen::Vector3d &t, double f, double rd, double r0, int direction){
    Eigen::Matrix<double,3,4> P;
    Eigen::Vector3d uh;
    Eigen::Matrix3d K = Eigen::Matrix3d::Zero();
    K.diagonal() << f, f, 1;
    // Eigen::Matrix3d Rv = Eigen::Matrix3d::Identity() + X_(v); 
    Eigen::Matrix3d Rv;
    double n = v.norm();
    if(n < 1e-15){
        Rv = Eigen::AngleAxis<double>(0, (Eigen::Vector3d() << 1,0,0).finished()).toRotationMatrix();
    }else{
        Rv = Eigen::AngleAxis<double>(n, v/n).toRotationMatrix();
    }

    uh = K* ( Rv * X + C );
    u = uh.head(2)/uh(2);

    double diff = 1e15;

    int niter = 0;
    
    Eigen::Vector3d temp;

    while(diff > 1e-10 ){
        double rc2 = u(0)*u(0) + u(1)*u(1);
        temp = K * ((Eigen::Matrix3d::Identity() + X_(u(direction) * w)) * Rv * X + C + u(direction) * t);
        temp = temp/temp(2);
        temp.head(2) *= (1 + rd * rc2);
        diff = (u - temp.head(2)).norm();
        u = temp.head(2);
        if(niter > 100){
            return WARN_NO_CONVERGENCE;
        }
        niter++;
    }


    return 0;

}

template<typename T>
bool isPoseApproxEqual(T pose1, T pose2){
    double err = (pose1.C-pose2.C).norm();   
    err += (pose1.v-pose2.v).norm();
    err += (pose1.t-pose2.t).norm();
    err += (pose1.w-pose2.w).norm();
    err += std::abs(pose1.f-pose2.f)/std::max(pose1.f,1.0);
    err += std::abs(pose1.rd-pose2.rd);
    std::cout << "err: " << err << "\n";
    return err < 1e-6;
}

template bool isPoseApproxEqual(RSDoublelinCameraPose pose1, RSDoublelinCameraPose pose2);
template bool isPoseApproxEqual(RSSinglelinCameraPose pose1, RSSinglelinCameraPose pose2);

