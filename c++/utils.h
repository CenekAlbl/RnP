#pragma once
#include <Eigen/Eigen>

// converts from unit quaternion to rotation matrix
void quat2R(Eigen::Vector4d const q, Eigen::Matrix3d & R) {
	double q1 = q(0);
	double q2 = q(1);
	double q3 = q(2);
	double q4 = q(3);

	R << q1*q1 + q2*q2 - q3*q3 - q4*q4, 2 * q2*q3 - 2 * q1*q4, 2 * q1*q3 + 2 * q2*q4,
		2 * q1*q4 + 2 * q2*q3, q1*q1 - q2*q2 + q3*q3 - q4*q4, 2 * q3*q4 - 2 * q1*q2,
		2 * q2*q4 - 2 * q1*q3, 2 * q1*q2 + 2 * q3*q4, q1*q1 - q2*q2 - q3*q3 + q4*q4;

}

// creates a skew symmetric matrix 
 Eigen::Matrix3d X_(Eigen::Vector3d const v) {
	 Eigen::Matrix3d S;
	 S << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
	 return S;
}