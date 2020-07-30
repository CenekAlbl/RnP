#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <chrono>

int lin_w_t_v_C_focal_6(Eigen::Matrix<double,7,3> X, Eigen::Matrix<double,7,2> u, Eigen::Vector3d vk, double r0){

Eigen::MatrixXd A = Eigen::MatrixXd::Zero(7,11);

A.col(0) = -X.col(2).cwiseProduct(u.col(0));
A.col(1) = -X.col(2).cwiseProduct(u.col(1));
A.col(2) = X.col(0).cwiseProduct(u.col(0)) + X.col(1).cwiseProduct(u.col(1));
A.col(3) = u.col(0).array() * (X.col(2).array() * (r0 - u.col(0).array())) - X.col(0).array() * (vk(1) * (r0 - u.col(0).array())) + X.col(1).array() * (vk(0) * (r0 - u.col(0).array()));
A.col(4) = u.col(1).array() * (X.col(2).array() * (r0 - u.col(0).array())) - X.col(0).array() * (vk(1) * (r0 - u.col(0).array())) + X.col(1).array() * (vk(0) * (r0 - u.col(0).array()));
A.col(5) = -u.col(0).array() * (X.col(0).array() * (r0 - u.col(0).array())) - X.col(1).array() * (vk(2) * (r0 - u.col(0).array())) + X.col(2).array() * (vk(1) * (r0 - u.col(0).array())) - u.col(1).array() * (X.col(1).array() * (r0 - u.col(0).array())) + X.col(0).array() * (vk(2) * (r0 - u.col(0).array()))- X.col(2).array() * (vk(0) * (r0 - u.col(0).array()));
A.col(6) = -u.col(1);
A.col(7) = u.col(0);
A.col(8) = u.col(1).array() * (r0 - u.col(0).array());
A.col(9) = -u.col(0).array() * (r0 - u.col(0).array());
A.col(10) = -X.col(1).cwiseProduct(u.col(0)) - X.col(0).cwiseProduct(u.col(1));

Eigen::MatrixXd nn = A.fullPivLu().kernel();
// std::cout << nn << "\n";

Eigen::MatrixXd n(11,4);

for ( size_t i = 0; i < 10; i++)
{
    n(i,0) = nn(i,0) - (nn(i,3) * nn(10,0)/nn(10,3));
    n(i,1) = nn(i,1) - (nn(i,3) * nn(10,1)/nn(10,3));
    n(i,2) = nn(i,2) - (nn(i,3) * nn(10,2)/nn(10,3));
    n(i,3) = nn(i,3)/nn(10,3);
}

Eigen::MatrixXd A0 = Eigen::MatrixXd::Zero(7,6);

A0.col(0) = X.col(0).array() * (n(5,0) * (r0 - u.col(0).array()) - n(2,0) + n(3,0) * vk(1) * (r0 - u.col(0).array())) - n(7,0) + n(9,0) * (r0 - u.col(0).array()) + X.col(2).array() * (n(0,0) - n(3,0) * (r0 - u.col(0).array()) + n(5,0) * vk(1) * (r0 - u.col(0).array())) - X.col(1).array() * (n(3,0) * vk(0) * (r0 - u.col(0).array()) + n(5,0) * vk(2) * (r0 - u.col(0).array()));
A0.col(1) = X.col(0).array() * (n(5,1) * (r0 - u.col(0).array()) - n(2,1) + n(3,1) * vk(1) * (r0 - u.col(0).array())) - n(7,1) + n(9,1) * (r0 - u.col(0).array()) + X.col(2).array() * (n(0,1) - n(3,1) * (r0 - u.col(0).array()) + n(5,1) * vk(1) * (r0 - u.col(0).array())) - X.col(1).array() * (n(3,1) * vk(0) * (r0 - u.col(0).array()) + n(5,1) * vk(2) * (r0 - u.col(0).array()));
A0.col(2) = X.col(0).array() * (n(5,2) * (r0 - u.col(0).array()) - n(2,2) + n(3,2) * vk(1) * (r0 - u.col(0).array())) - n(7,2) + n(9,2) * (r0 - u.col(0).array()) + X.col(2).array() * (n(0,2) - n(3,2) * (r0 - u.col(0).array()) + n(5,2) * vk(1) * (r0 - u.col(0).array())) - X.col(1).array() * (n(3,2) * vk(0) * (r0 - u.col(0).array()) + n(5,2) * vk(2) * (r0 - u.col(0).array()));
A0.col(5) = X.col(0).array() * (n(5,3) * (r0 - u.col(0).array()) - n(2,3) + n(3,3) * vk(1) * (r0 - u.col(0).array())) - n(7,3) + n(9,3) * (r0 - u.col(0).array()) + X.col(2).array() * (n(0,3) - n(3,3) * (r0 - u.col(0).array()) + n(5,3) * vk(1) * (r0 - u.col(0).array())) - X.col(1).array() * (n(3,3) * vk(0) * (r0 - u.col(0).array()) + n(5,3) * vk(2) * (r0 - u.col(0).array())+1);

Eigen::MatrixXd A1 = Eigen::MatrixXd::Zero(7,6);

A1.col(0) = -u.col(1).array() * (X.col(1).array() * (n(3,0) * (r0 - u.col(0).array()) - n(0,0) + n(4,0) * vk(2) * (r0 - u.col(0).array())) + X.col(0).array() * (n(1,0) - n(4,0) * (r0 - u.col(0).array()) + n(3,0) * vk(2) * (r0-u.col(0).array())) - X.col(2).array() * (n(3,0) * vk(0) * (r0 - u.col(0).array()) + n(4,0) * vk(1) * (r0 - u.col(0).array())));
A1.col(1) = -u.col(1).array() * (X.col(1).array() * (n(3,1) * (r0 - u.col(0).array()) - n(0,1) + n(4,1) * vk(2) * (r0 - u.col(0).array())) + X.col(0).array() * (n(1,1) - n(4,1) * (r0 - u.col(0).array()) + n(3,1) * vk(2) * (r0-u.col(0).array())) - X.col(2).array() * (n(3,1) * vk(0) * (r0 - u.col(0).array()) + n(4,1) * vk(1) * (r0 - u.col(0).array())));
A1.col(2) = -u.col(1).array() * (X.col(1).array() * (n(3,2) * (r0 - u.col(0).array()) - n(0,2) + n(4,2) * vk(2) * (r0 - u.col(0).array())) + X.col(0).array() * (n(1,2) - n(4,2) * (r0 - u.col(0).array()) + n(3,2) * vk(2) * (r0-u.col(0).array())) - X.col(2).array() * (n(3,2) * vk(0) * (r0 - u.col(0).array()) + n(4,2) * vk(1) * (r0 - u.col(0).array())));
A1.col(5) = -u.col(1).array() * (X.col(1).array() * (n(3,3) * (r0 - u.col(0).array()) - n(0,3) + n(4,3) * vk(2) * (r0 - u.col(0).array())) + X.col(0).array() * (n(1,3) - n(4,3) * (r0 - u.col(0).array()) + n(3,3) * vk(2) * (r0-u.col(0).array())) - X.col(2).array() * (n(3,3) * vk(0) * (r0 - u.col(0).array()) + n(4,3) * vk(1) * (r0 - u.col(0).array())));
A1.col(3) = u.col(1);
A1.col(4) = -u.col(1).array() * (r0 - u.col(0).array());


Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
ges.compute(A0.topRows(6), A1.topRows(6),true);
Eigen::MatrixXcd f = ges.eigenvalues();
Eigen::MatrixXcd xx = ges.eigenvectors();
std::cout << f << "\n";
std::cout << xx << "\n";


// std::cout << A0 << "\n";
return 0;

}

int main(int argc, char ** argv){





// std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
// for (size_t i = 0; i < 10000; i++)
// {
    lin_w_t_v_C_focal_6(Eigen::Matrix<double,7,3>::Random(), Eigen::Matrix<double,7,2>::Random(), Eigen::Vector3d::Random(), 0 );
// }

// std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

// std::cout << "Time difference = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/10000) << "[µs]" << std::endl;
    

}