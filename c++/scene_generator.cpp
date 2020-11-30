#include "scene_generator.h"
#include <iostream>
#include <Eigen/Geometry> 



Eigen::Vector2d rsDoubleLinProjection(const Eigen::Vector3d &X, Eigen::Vector2d &u, const Eigen::Vector3d &v, const Eigen::Vector3d &C, const Eigen::Vector3d &w, const Eigen::Vector3d &t, double f, double rd, double r0, int direction){
    // First initiate u with a global shutter projection
    Eigen::Matrix<double,3,4> P;
    Eigen::Vector3d uh;
    uh << u, 1;
    Eigen::Matrix3d K;
    K.diagonal() << f, f, 1;
    Eigen::Matrix3d Rv = Eigen::Matrix3d::Identity() + X_(v); 
    Eigen::Vector3d ut = K* ( Rv * X + C );
    ut = ut/ut(2);

    double rc2 = ut(0)*ut(0) + ut(1)*ut(1);

    double diff = 1e15;

    while(diff > 1e-10){
        Eigen::Vector3d temp = K * ((Eigen::Matrix3d::Identity() + X_(ut(direction) * w)) * Rv * X + C + ut(direction) * t);
        temp = temp/temp(2);
        diff = (ut - temp).norm();
        ut = temp;
        
        std::cout << diff << "\n";
    }

    u = ut.head(2);
    
}

Eigen::Vector2d rsSingleLinProjection(const Eigen::Vector3d &X, Eigen::Vector2d &u, const Eigen::Vector3d &v, const Eigen::Vector3d &C, const Eigen::Vector3d &w, const Eigen::Vector3d &t, double f, double rd, double r0, int direction){
    Eigen::Matrix<double,3,4> P;
    Eigen::Vector3d uh;
    uh << u, 1;
    Eigen::Matrix3d K;
    K.diagonal() << f, f, 1;
    // Eigen::Matrix3d Rv = Eigen::Matrix3d::Identity() + X_(v); 
    Eigen::Matrix3d Rv;
    Rv = Eigen::AngleAxis<double>(v.norm(), v).toRotationMatrix();
    
    Eigen::Vector3d ut = K* ( Rv * X + C );
    ut = ut/ut(2);

    double rc2 = ut(0)*ut(0) + ut(1)*ut(1);

    double diff = 1e15;

    while(diff > 1e-10){
        Eigen::Vector3d temp = K * ((Eigen::Matrix3d::Identity() + X_(ut(direction) * w)) * Rv * X + C + ut(direction) * t);
        temp = temp/temp(2);
        diff = (ut - temp).norm();
        ut = temp;
        
        std::cout << diff << "\n";
    }

    u = ut.head(2);

}

int main(int argc, char ** argv){

    Eigen::MatrixXd X(3,7);
    X << -0.004985829085851,   0.751978827846046,   0.703321607967492,   0.324836770098960,  -0.150331413612924,   0.362116095206949,   0.771109175334294,
  -0.050017318391826,   0.544760816964284,   0.455853510042298,   0.189075200710939,  -0.196453169329558,   0.269614531134675,   0.458682324161694,
  -1.109345725478448,   0.252824851009254,  -1.044324845933275,  -0.994706943626732,  -2.140253201676392,   0.295623763317194,  -2.126185255504355;

   Eigen::Vector3d v,w,C,t;
   v << 0,0,0;

   w << 0.334626153478911e-03,
  -0.116809228874349e-03,
  -0.185420984404916e-03;

   C << -0.738827831825120,
  -1.479842410229943,
  -3.846628250527604;

   t << 0.186789143458615e-03,
   0.019215829596831e-03,
  -0.068851781232604e-03;

   Eigen::Vector2d u;
   double f = 1500;
   double rd = -1.111111111111111e-07;

    for (int i = 0; i < X.cols(); i++)
    {
        rsSingleLinProjection(X.col(i), u, v, C, w, t, f, rd, 0, 0);
        std::cout << u << "\n";

        double err =  calcErrAlgebraicRnPFocalRadialSingleLin(v, C, w, t, f, rd, 0, X.col(i), u, 0);

        std::cout << "error: " <<err << "\n";
    }
    
   

   

}