#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <chrono>

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


int lin_w_t_v_C_focal_radial(Eigen::Matrix<double,7,3> X, Eigen::Matrix<double,7,2> u, Eigen::Vector3d vk, double r0){


Eigen::MatrixXd A = Eigen::MatrixXd::Random(7,11);

Eigen::MatrixXd nn = A.fullPivLu().kernel();

Eigen::MatrixXd AA = Eigen::MatrixXd::Random(7,14);

Eigen::MatrixXd AR = AA.leftCols(7).lu().solve(AA);

Eigen::VectorXd data = Eigen::VectorXd::Random(35);

solver_R7Pfr(data);

return 0;

}

int main(int argc, char ** argv){





std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
for (size_t i = 0; i < 10000; i++)
{
    lin_w_t_v_C_focal_radial(Eigen::Matrix<double,7,3>::Random(), Eigen::Matrix<double,7,2>::Random(), Eigen::Vector3d::Random(), 0 );
}

std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

std::cout << "Time difference = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/10000) << "[µs]" << std::endl;
    

}