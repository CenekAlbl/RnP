#include <Eigen/Eigen>
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


Eigen::MatrixXd null = A.fullPivLu().kernel();

// A = [ -X(:,3).*u(:,1), 
// -X(:,3).*u(:,2), 
// X(:,1).*u(:,1) + X(:,2).*u(:,2), 
// u(:,1).*(X(:,3).*(r0 - u(:,1)) - X(:,1).*vk(2).*(r0 - u(:,1)) + X(:,2).*vk(1).*(r0 - u(:,1))), 
// u(:,2).*(X(:,3).*(r0 - u(:,1)) - X(:,1).*vk(2).*(r0 - u(:,1)) + X(:,2).*vk(1).*(r0 - u(:,1))), 
// - u(:,1).*(X(:,1).*(r0 - u(:,1)) - X(:,2).*vk(3).*(r0 - u(:,1)) + X(:,3).*vk(2).*(r0 - u(:,1))) - u(:,2).*(X(:,2).*(r0 - u(:,1)) + X(:,1).*vk(3).*(r0 - u(:,1)) - X(:,3).*vk(1).*(r0 - u(:,1))), 
// -u(:,2), 
// u(:,1), 
// u(:,2).*(r0 - u(:,1)), 
// -u(:,1).*(r0 - u(:,1)), 
// X(:,2).*u(:,1) - X(:,1).*u(:,2)];

}

int main(int argc, char ** argv){





std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
for (size_t i = 0; i < 10000; i++)
{
    lin_w_t_v_C_focal_6(Eigen::Matrix<double,7,3>::Random(), Eigen::Matrix<double,7,2>::Random(), Eigen::Vector3d::Random(), 0 );
}

std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

std::cout << "Time difference = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/10000) << "[µs]" << std::endl;
    

}