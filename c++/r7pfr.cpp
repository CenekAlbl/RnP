#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <chrono>
#include "RnP.h"


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

static void setup_elimination_template(Eigen::Matrix<double,35,1> const input, Eigen::Matrix<double,26,26> &C0, Eigen::Matrix<double,26,10> &C1)
{
  double coeffs[34];
  int i0;
  static const short iv0[151] = { 0, 5, 27, 35, 54, 57, 60, 78, 79, 80, 81, 83,
    85, 87, 102, 105, 108, 113, 116, 130, 133, 153, 157, 162, 163, 179, 184, 187,
    195, 208, 211, 233, 234, 235, 236, 237, 240, 241, 247, 256, 257, 259, 260,
    262, 263, 264, 265, 270, 272, 280, 286, 289, 297, 313, 318, 319, 322, 323,
    324, 326, 335, 340, 341, 351, 363, 364, 366, 367, 374, 375, 376, 377, 378,
    380, 389, 394, 400, 401, 402, 404, 405, 408, 410, 426, 427, 428, 430, 452,
    453, 454, 456, 458, 459, 460, 482, 484, 485, 486, 496, 497, 499, 502, 503,
    507, 513, 515, 518, 519, 522, 525, 535, 554, 555, 565, 570, 573, 578, 579,
    580, 581, 591, 595, 596, 602, 606, 607, 613, 616, 617, 618, 622, 625, 628,
    630, 631, 632, 633, 634, 635, 636, 638, 639, 642, 643, 644, 647, 648, 654,
    665, 668, 670 };

  static const signed char iv1[151] = { 0, 0, 0, 0, 0, 12, 19, 3, 2, 1, 0, 13,
    19, 2, 19, 4, 0, 4, 0, 5, 1, 19, 5, 1, 21, 28, 2, 14, 19, 7, 2, 19, 8, 7, 5,
    13, 2, 22, 21, 19, 29, 28, 9, 6, 4, 2, 17, 19, 2, 19, 10, 6, 19, 10, 6, 25,
    21, 28, 5, 1, 32, 7, 14, 22, 29, 11, 10, 17, 22, 29, 7, 25, 2, 19, 32, 17,
    24, 31, 9, 4, 25, 6, 32, 25, 32, 10, 6, 26, 33, 11, 17, 32, 25, 10, 18, 33,
    26, 11, 3, 12, 15, 22, 14, 20, 2, 19, 29, 27, 4, 16, 19, 20, 12, 0, 27, 3, 0,
    20, 21, 13, 1, 27, 28, 12, 24, 16, 20, 0, 4, 27, 31, 9, 13, 4, 24, 25, 17,
    20, 27, 3, 0, 21, 1, 6, 28, 31, 32, 16, 24, 4, 31 };

  static const unsigned char uv0[87] = { 2U, 3U, 6U, 13U, 19U, 21U, 22U, 25U,
    28U, 29U, 30U, 31U, 39U, 41U, 43U, 44U, 46U, 51U, 54U, 55U, 65U, 66U, 68U,
    69U, 70U, 77U, 84U, 86U, 87U, 97U, 99U, 100U, 102U, 105U, 110U, 111U, 123U,
    125U, 126U, 127U, 136U, 149U, 151U, 152U, 160U, 162U, 164U, 165U, 170U, 171U,
    172U, 173U, 174U, 175U, 176U, 177U, 178U, 180U, 183U, 188U, 189U, 192U, 193U,
    194U, 196U, 198U, 199U, 200U, 201U, 203U, 204U, 205U, 214U, 222U, 224U, 225U,
    226U, 227U, 229U, 230U, 238U, 248U, 249U, 250U, 251U, 252U, 254U };

  static const signed char iv2[87] = { 8, 15, 14, 23, 7, 22, 29, 30, 9, 16, 14,
    18, 24, 22, 19, 2, 29, 31, 11, 18, 26, 14, 29, 22, 7, 33, 12, 23, 15, 3, 20,
    27, 30, 8, 13, 23, 5, 21, 28, 30, 15, 8, 23, 30, 15, 16, 26, 18, 12, 23, 27,
    20, 3, 9, 30, 24, 31, 33, 11, 17, 26, 23, 30, 8, 13, 28, 21, 5, 10, 25, 32,
    33, 18, 15, 30, 23, 8, 11, 26, 33, 18, 16, 26, 31, 24, 9, 33 };

  compute_coeffs(input.data(), coeffs);
  memset(&C0.data()[0], 0, 676U * sizeof(double));
  memset(&C1.data()[0], 0, 260U * sizeof(double));
  for (i0 = 0; i0 < 151; i0++) {
    C0.data()[iv0[i0]] = coeffs[iv1[i0]];
  }

  for (i0 = 0; i0 < 87; i0++) {
    C1.data()[uv0[i0]] = coeffs[iv2[i0]];
  }
}

int solver_R7Pfr(Eigen::Matrix<double,35,1> const data){

Eigen::Matrix<double, 26, 26> C0;
Eigen::Matrix<double, 26, 10> C1;

std::cout << "setting up elimination template \n";

setup_elimination_template(data, C0, C1);

std::cout << "elimination template set\n";


C1 = C0.lu().solve(C1);

std::cout << C1 << "\n";

Eigen::Matrix<double, 17, 10> RR;

RR.topRows(7) = -C1.bottomRows(7);
RR.bottomRows(10) = Eigen::MatrixXd::Identity(10,10);

std::vector<int> ind{1,2,9,3,4,11,5,6,14,7};

Eigen::Matrix<double, 10, 10> AM;




// AM = RR(ind,Eigen::all)

Eigen::EigenSolver<Eigen::MatrixXd> eig(AM);
std::cout << eig.eigenvalues() << "\n";
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


Eigen::Matrix<double,35,1> data; 
data << AR.row(0).tail(7).transpose(), AR.row(1).tail(7).transpose(), AR.row(2).tail(7).transpose(), AR.row(3).tail(7).transpose(), AR.row(4).tail(7).transpose();


solver_R7Pfr(data);

return 0;

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