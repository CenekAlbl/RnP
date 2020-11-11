#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <chrono>
#include "RnP.h"

int solver_R7Pfr(Eigen::MatrixXd data){

Eigen::MatrixXd C0 = Eigen::MatrixXd::Random(26,26);
Eigen::MatrixXd C1 = Eigen::MatrixXd::Random(26,10);

C0.lu().solve(C1);

Eigen::MatrixXd AM = Eigen::MatrixXd::Random(10,10);

Eigen::EigenSolver<Eigen::MatrixXd> eig(AM);
eig.eigenvalues();
eig.eigenvectors();

return 0;
}


int lin_w_t_v_C_focal_radial(Eigen::Matrix<double,7,3> X, Eigen::Matrix<double,7,2> u, Eigen::Vector3d vk, double r0, RSDoublelinCameraPoseVector * results){


Eigen::MatrixXd A = Eigen::MatrixXd::Zero(7,11);


// A = [ 
//     -X.col(2).array().*u(:,1), 
//     -X(:,3).*u(:,2), 
//     X(:,1).*u(:,1) + X(:,2).*u(:,2), 
//     u(:,1).*(X(:,3).*(r0 - u(:,1)) - X(:,1).*vk(2).*(r0 - u(:,1)) + X(:,2).*vk(0).*(r0 - u(:,1))), 
//     u(:,2).*(X(:,3).*(r0 - u(:,1)) - X(:,1).*vk(2).*(r0 - u(:,1)) + X(:,2).*vk(1).*(r0 - u(:,1))), 
//     - u(:,1).*(X(:,1).*(r0 - u(:,1)) - X(:,2).*vk(3).*(r0 - u(:,1)) + X(:,3).*vk(2).*(r0 - u(:,1))) - u(:,2).*(X(:,2).*(r0 - u(:,1)) + X(:,1).*vk(3).*(r0 - u(:,1)) - X(:,3).*vk(1).*(r0 - u(:,1))), 
//     -u(:,2), 
//     u(:,1), 
//     u(:,2).*(r0 - u(:,1)), 
//     -u(:,1).*(r0 - u(:,1)), 
//     X(:,2).*u(:,1) - X(:,1).*u(:,2)];

A.col(0) = -X.col(2).array() * u.col(0).array();
A.col(1) = -X.col(2).array() * u.col(1).array();
A.col(2) = X.col(0).array() * u.col(0).array() + X.col(1).array() * u.col(1).array();
A.col(3) = u.col(0).array() * (X.col(2).array() * (r0 - u.col(0).array()) - X.col(0).array() * vk(1) * (r0 - u.col(0).array()) + X.col(1).array() * vk(0) * (r0 - u.col(0).array()));
A.col(4) = u.col(1).array() * (X.col(2).array() * (r0 - u.col(0).array()) - X.col(0).array() * vk(1) * (r0 - u.col(0).array()) + X.col(1).array() * vk(0) * (r0 - u.col(0).array()));
A.col(5) = -u.col(0).array() * (X.col(0).array() * (r0 - u.col(0).array()) - X.col(1).array() * vk(2) * (r0 - u.col(0).array()) + X.col(2).array() * vk(1) * (r0 - u.col(0).array())) - u.col(1).array() * (X.col(1).array() * (r0 - u.col(0).array()) + X.col(0).array() * vk(2) * (r0 - u.col(0).array()) - X.col(2).array() * vk(0) * (r0 - u.col(0).array()));
A.col(6) = -u.col(1);
A.col(7) = u.col(0);
A.col(8) = u.col(1).array() * (r0 - u.col(0).array());
A.col(9) = -u.col(0).array() * (r0 - u.col(0).array());
A.col(10) = X.col(1).array() * u.col(0).array() - X.col(0).array() * u.col(1).array();



Eigen::MatrixXd nn = A.fullPivLu().kernel();

Eigen::MatrixXd n(10,4);

for ( int i = 0; i < 10; i++)
{
    n(i,0) = nn(i,0) - (nn(i,3) * nn(10,0)/nn(10,3));
    n(i,1) = nn(i,1) - (nn(i,3) * nn(10,1)/nn(10,3));
    n(i,2) = nn(i,2) - (nn(i,3) * nn(10,2)/nn(10,3));
    n(i,3) = nn(i,3)/nn(10,3);
}


n << -0.494882181407799,   0.411226263985881,   0.318345335859557,   0.927819398937113,
   0.993052896562505,   0.687977958425696,   0.768893995395170,  -1.793033698569230,
  -0.518424090047564,  -0.923289492859002,   0.472050129410013,   1.129772420556994,
  -0.000818047233864,  -0.000355381576460,   0.000238057785641,   0.000167456002760,
  -0.001271430972275,  -0.001007920402057,  -0.000066867094954,   0.001183167323652,
   0.000238330557519,   0.001527458270602,  -0.000768958234596,  -0.001320982654342,
   0.551245446952516,   0.337315602380641,   0.619130851038882,  -1.908967112082546,
   0.550113067243128,  -0.402009451439514,  -0.178094285891012,   0.809029947913194,
  -0.000061795565197,  -0.000822741118985,  -0.001197822291309,   0.000595842167120,
  -0.002044071545096,   0.000806256182446,   0.000662084265002,   0.000148745747078;


Eigen::VectorXd r2 = u.col(0).array().pow(2)+u.col(1).array().pow(2);


Eigen::MatrixXd AA(7,14);

AA.col(0) = -u.col(1).array()*(X.col(1).array()*(n(3,0)*(r0 - u.col(0).array()) - n(0,0) + n(4,0)*vk(2)*(r0 - u.col(0).array())) + X.col(0).array()*(n(1,0) - n(4,0)*(r0 - u.col(0).array()) + n(3,0)*vk(2)*(r0 - u.col(0).array())) - X.col(2).array()*(n(3,0)*vk(0)*(r0 - u.col(0).array()) + n(4,0)*vk(1)*(r0 - u.col(0).array())));
AA.col(1) = r2.array()*(X.col(0).array()*(n(5,0)*(r0 - u.col(0).array()) - n(2,0) + n(3,0)*vk(1)*(r0 - u.col(0).array())) - n(7,0) + n(9,0)*(r0 - u.col(0).array()) + X.col(2).array()*(n(0,0) - n(3,0)*(r0 - u.col(0).array()) + n(5,0)*vk(1)*(r0 - u.col(0).array())) - X.col(1).array()*(n(3,0)*vk(0)*(r0 - u.col(0).array()) + n(5,0)*vk(2)*(r0 - u.col(0).array())));
AA.col(2) = X.col(0).array()*(n(5,0)*(r0 - u.col(0).array()) - n(2,0) + n(3,0)*vk(1)*(r0 - u.col(0).array())) - n(7,0) + n(9,0)*(r0 - u.col(0).array()) + X.col(2).array()*(n(0,0) - n(3,0)*(r0 - u.col(0).array()) + n(5,0)*vk(1)*(r0 - u.col(0).array())) - X.col(1).array()*(n(3,0)*vk(0)*(r0 - u.col(0).array()) + n(5,0)*vk(2)*(r0 - u.col(0).array()));
AA.col(3) = -u.col(1).array()*(X.col(1).array()*(n(3,1)*(r0 - u.col(0).array()) - n(0,1) + n(4,1)*vk(2)*(r0 - u.col(0).array())) + X.col(0).array()*(n(1,1) - n(4,1)*(r0 - u.col(0).array()) + n(3,1)*vk(2)*(r0 - u.col(0).array())) - X.col(2).array()*(n(3,1)*vk(0)*(r0 - u.col(0).array()) + n(4,1)*vk(1)*(r0 - u.col(0).array()))); 
AA.col(4) = r2.array()*(X.col(0).array()*(n(5,1)*(r0 - u.col(0).array()) - n(2,1) + n(3,1)*vk(1)*(r0 - u.col(0).array())) - n(7,1) + n(9,1)*(r0 - u.col(0).array()) + X.col(2).array()*(n(0,1) - n(3,1)*(r0 - u.col(0).array()) + n(5,1)*vk(1)*(r0 - u.col(0).array())) - X.col(1).array()*(n(3,1)*vk(0)*(r0 - u.col(0).array()) + n(5,1)*vk(2)*(r0 - u.col(0).array())));
AA.col(7) = X.col(0).array()*(n(5,1)*(r0 - u.col(0).array()) - n(2,1) + n(3,1)*vk(1)*(r0 - u.col(0).array())) - n(7,1) + n(9,1)*(r0 - u.col(0).array()) + X.col(2).array()*(n(0,1) - n(3,1)*(r0 - u.col(0).array()) + n(5,1)*vk(1)*(r0 - u.col(0).array())) - X.col(1).array()*(n(3,1)*vk(0)*(r0 - u.col(0).array()) + n(5,1)*vk(2)*(r0 - u.col(0).array()));
AA.col(8) = -u.col(1).array()*(X.col(1).array()*(n(3,2)*(r0 - u.col(0).array()) - n(0,2) + n(4,2)*vk(2)*(r0 - u.col(0).array())) + X.col(0).array()*(n(1,2) - n(4,2)*(r0 - u.col(0).array()) + n(3,2)*vk(2)*(r0 - u.col(0).array())) - X.col(2).array()*(n(3,2)*vk(0)*(r0 - u.col(0).array()) + n(4,2)*vk(1)*(r0 - u.col(0).array())));
AA.col(9) = r2.array()*(X.col(0).array()*(n(5,2)*(r0 - u.col(0).array()) - n(2,2) + n(3,2)*vk(1)*(r0 - u.col(0).array())) - n(7,2) + n(9,2)*(r0 - u.col(0).array()) + X.col(2).array()*(n(0,2) - n(3,2)*(r0 - u.col(0).array()) + n(5,2)*vk(1)*(r0 - u.col(0).array())) - X.col(1).array()*(n(3,2)*vk(0)*(r0 - u.col(0).array()) + n(5,2)*vk(2)*(r0 - u.col(0).array())));
AA.col(10) = X.col(0).array()*(n(5,2)*(r0 - u.col(0).array()) - n(2,2) + n(3,2)*vk(1)*(r0 - u.col(0).array())) - n(7,2) + n(9,2)*(r0 - u.col(0).array()) + X.col(2).array()*(n(0,2) - n(3,2)*(r0 - u.col(0).array()) + n(5,2)*vk(1)*(r0 - u.col(0).array())) - X.col(1).array()*(n(3,2)*vk(0)*(r0 - u.col(0).array()) + n(5,2)*vk(2)*(r0 - u.col(0).array()));
AA.col(5) = u.col(1).array();
AA.col(6) = -u.col(1).array()*(r0 - u.col(0).array());
AA.col(11) = -u.col(1).array()*(X.col(1).array()*(n(3,3)*(r0 - u.col(0).array()) - n(0,3) + n(4,3)*vk(2)*(r0 - u.col(0).array())) + X.col(0).array()*(n(1,3) - n(4,3)*(r0 - u.col(0).array()) + n(3,3)*vk(2)*(r0 - u.col(0).array())) - X.col(2).array()*(n(3,3)*vk(0)*(r0 - u.col(0).array()) + n(4,3)*vk(1)*(r0 - u.col(0).array()) + 1));
AA.col(12) = r2.array()*(X.col(0).array()*(n(5,3)*(r0 - u.col(0).array()) - n(2,3) + n(3,3)*vk(1)*(r0 - u.col(0).array())) - n(7,3) + n(9,3)*(r0 - u.col(0).array()) + X.col(2).array()*(n(0,3) - n(3,3)*(r0 - u.col(0).array()) + n(5,3)*vk(1)*(r0 - u.col(0).array())) - X.col(1).array()*(n(3,3)*vk(0)*(r0 - u.col(0).array()) + n(5,3)*vk(2)*(r0 - u.col(0).array()) + 1));
AA.col(13) = X.col(0).array()*(n(5,3)*(r0 - u.col(0).array()) - n(2,3) + n(3,3)*vk(1)*(r0 - u.col(0).array())) - n(7,3) + n(9,3)*(r0 - u.col(0).array()) + X.col(2).array()*(n(0,3) - n(3,3)*(r0 - u.col(0).array()) + n(5,3)*vk(1)*(r0 - u.col(0).array())) - X.col(1).array()*(n(3,3)*vk(0)*(r0 - u.col(0).array()) + n(5,3)*vk(2)*(r0 - u.col(0).array()) + 1);

Eigen::MatrixXd AR = AA.leftCols(7).lu().solve(AA);

Eigen::VectorXd data << AA.col(0).tail(7), AA.col(1).tail(7), AA.col(2).tail(7), AA.col(3).tail(7), AA.col(4).tail(7);

solver_R7Pfr(data);

return 0;

}


static void compute_coeffs(const double data[35], double coeffs[34])
{
  coeffs[0] = -data[15];
  coeffs[1] = -data[16];
  coeffs[2] = -data[14];
  coeffs[3] = data[1] - data[17];
  coeffs[4] = -data[18];
  coeffs[5] = data[2];
  coeffs[6] = -data[19];
  coeffs[7] = data[0];
  coeffs[8] = data[3];
  coeffs[9] = data[4] - data[20];
  coeffs[10] = data[5];
  coeffs[11] = data[6];
  coeffs[12] = data[8];
  coeffs[13] = data[9] - data[17];
  coeffs[14] = data[7];
  coeffs[15] = data[10];
  coeffs[16] = data[11];
  coeffs[17] = data[12] - data[20];
  coeffs[18] = data[13];
  coeffs[19] = 1.0;
  coeffs[20] = data[22];
  coeffs[21] = data[23];
  coeffs[22] = data[21];
  coeffs[23] = data[24];
  coeffs[24] = data[25];
  coeffs[25] = data[26];
  coeffs[26] = data[27];
  coeffs[27] = data[29];
  coeffs[28] = data[30];
  coeffs[29] = data[28];
  coeffs[30] = data[31];
  coeffs[31] = data[32];
  coeffs[32] = data[33];
  coeffs[33] = data[34];
}

int main(int argc, char ** argv){

Eigen::MatrixXd X(3,7);
Eigen::MatrixXd u(2,7);

X << 0.930335618971609,  -0.412757128555297,  -0.300634768231611,  -0.472768038810491,   0.904183433628433,  -0.746899095560691,   0.782781495505683,
   0.784251388397344,   0.016563901098688,  -0.676295260848509,   0.944897105666936,   0.323358545689461,  -0.370185130712374,   0.202565077125717,
   0.881573031546951,   0.786836620540149,  -0.092437951112491,  -0.339174391869371,   0.821280370700545,   0.470064724843527,  -0.669632274378342;

u <<  858.8978263748063,   778.5936238794234,   304.1044457640757,   585.7861538645086,   674.3291756530226,   618.6345560318252,   33.8594439672178,
  -577.5015321917176,  -20.4871105409913,  -10.7388316467659,  -593.0554378121407,  -267.4839868749193,   43.4202366607074,  -555.0719801311780;

RSDoublelinCameraPoseVector results;

std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
for (size_t i = 0; i < 1; i++)
{
    lin_w_t_v_C_focal_radial(X.transpose(), u.transpose(), Eigen::Vector3d::Zero(), 0 , &results);
}

std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

std::cout << "Time difference = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/10000) << "[µs]" << std::endl;
    

}