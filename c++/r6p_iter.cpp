#include "rnp.h"
#include <iostream>

// def linwtvC(X,u,vk):    
//     A = np.zeros((12,13))
//     A[:,[0]] = np.vstack((X[:,[2]] + X[:,[1]]*u[:,[1]], -X[:,[1]]*u[:,[0]]))

//     A[:,[1]] = np.vstack((-X[:,[0]]*u[:,[1]], X[:,[2]] + X[:,[0]]*u[:,[0]]))

//     A[:,[2]] = np.vstack((-X[:,[0]], -X[:,[1]]))

//     A[:,[3]] = np.vstack(
    // (X[:,[0]]*vk[1]*( - u[:,[0]]) - X[:,[2]]*( - u[:,[0]]) - u[:,[1]]*(X[:,[1]]*( - u[:,[0]]) + X[:,[0]]*vk[2]*( - u[:,[0]]) - X[:,[2]]*vk[0]*( - u[:,[0]])) - X[:,[1]]*vk[0]*( - u[:,[0]]), 
    // u[:,[0]]*(X[:,[1]]*( - u[:,[0]]) + X[:,[0]]*vk[2]*( - u[:,[0]]) - X[:,[2]]*vk[0]*( - u[:,[0]]))))

//     A[:,[4]] = np.vstack((
        // u[:,[1]]*(X[:,[0]]*(-u[:,[0]]) - X[:,[1]]*vk[2]*( - u[:,[0]]) + X[:,[2]]*vk[1]*( - u[:,[0]])), 
        // X[:,[0]]*vk[1]*( - u[:,[0]]) - X[:,[2]]*( - u[:,[0]]) - u[:,[0]]*(X[:,[0]]*( - u[:,[0]]) - X[:,[1]]*vk[2]*( - u[:,[0]]) + X[:,[2]]*vk[1]*( - u[:,[0]])) - X[:,[1]]*vk[0]*( - u[:,[0]])))

//     A[:,[5]] = np.vstack((
        // X[:,[0]]*( - u[:,[0]]) - X[:,[1]]*vk[2]*( - u[:,[0]]) + X[:,[2]]*vk[1]*( - u[:,[0]])
        // X[:,[1]]*( - u[:,[0]]) + X[:,[0]]*vk[2]*( - u[:,[0]]) - X[:,[2]]*vk[0]*( - u[:,[0]])))

//     A[:,[6]] = np.vstack((np.zeros((6,1)), np.ones((6,1))))
//     A[:,[7]] = np.vstack((-np.ones((6,1)), np.zeros((6,1))))
//     A[:,[8]] = np.vstack((u[:,[1]], -u[:,[0]]))
//     A[:,[9]] = np.vstack((np.zeros((6,1)), u[:,[0]]))
//     A[:,[10]] = np.vstack((-u[:,[0]], np.zeros((6,1))))
//     A[:,[11]] = np.vstack((-u[:,[1]]*( - u[:,[0]]), u[:,[0]]*( - u[:,[0]])))
//     A[:,[12]] = np.vstack((X[:,[2]]*u[:,[1]] - X[:,[1]], X[:,[0]] - X[:,[2]]*u[:,[0]]))
//     n = null_space(A)
//     s = n[12,-1]
//     v = n[0:3,[-1]]/s
//     w = n[3:6,[-1]]/s
//     C = n[6:9,[-1]]/s
//     t = n[9:12,[-1]]/s
//     return w,t,v,C

int R6PIter(const Eigen::MatrixXd & X, const Eigen::MatrixXd & u, const Eigen::Vector3d & vk, double r0, RSDoublelinCameraPoseVector * results){
    Eigen::Matrix<double,12,13> A = Eigen::Matrix<double,12,13>::Zero();
    A.col(0).head(6) << X.row(2).transpose().array() + X.row(1).transpose().array() * u.row(1).transpose().array();  
    A.col(0).tail(6) << -X.row(1).transpose().array() * u.row(0).transpose().array();
    A.col(1).head(6) << -X.row(0).transpose().array() * u.row(1).transpose().array();
    A.col(1).tail(6) << X.row(2).transpose().array() + X.row(0).transpose().array() * u.row(0).transpose().array();
    A.col(2).head(6) << - X.row(0).transpose().array();
    A.col(2).tail(6) << - X.row(1).transpose().array();
    A.col(3).head(6) << X.row(0).transpose().array() * vk(1) * (-u.row(0).transpose().array()) - X.row(2).transpose().array() * (- u.row(0).transpose().array()) - u.row(1).transpose().array() * (X.row(1).transpose().array() * ( -u.row(0).transpose().array()) + X.row(0).transpose().array() * vk(2) * (-u.row(0).transpose().array()) - X.row(2).transpose().array() * vk(0) * (-u.row(0).transpose().array())) - X.row(1).transpose().array() * vk(0) * (-u.row(0).transpose().array());
    A.col(3).tail(6) << u.row(0).transpose().array() * (X.row(1).transpose().array() * (-u.row(0).transpose().array()) + X.row(0).transpose().array() * vk(2) * (-u.row(0).transpose().array()) - X.row(2).transpose().array() * vk(0) * (-u.row(0).transpose().array()));
    A.col(4).head(6) << u.row(1).transpose().array() * (X.row(0).transpose().array() * (-u.row(0).transpose().array()) - X.row(1).transpose().array() * vk(2) * (-u.row(0).transpose().array()) + X.row(2).transpose().array() * vk(1) * (-u.row(0).transpose().array()));
    A.col(4).tail(6) << X.row(0).transpose().array() * vk(1) * (-u.row(0).transpose().array()) - X.row(2).transpose().array() * (-u.row(0).transpose().array()) - u.row(0).transpose().array() * (X.row(0).transpose().array() * (-u.row(0).transpose().array()) -X.row(1).transpose().array() * vk(2) * (-u.row(0).transpose().array()) + X.row(2).transpose().array() * vk(1) * (-u.row(0).transpose().array())) - X.row(1).transpose().array() * vk(0) * (-u.row(0).transpose().array());
    A.col(5).head(6) << X.row(0).transpose().array() * (-u.row(0).transpose().array()) - X.row(1).transpose().array() * vk(2) * (-u.row(0).transpose().array()) + X.row(2).transpose().array() * vk(1) * (-u.row(0).transpose().array());
    A.col(5).tail(6) << X.row(1).transpose().array() * (-u.row(0).transpose().array()) + X.row(0).transpose().array() * vk(2) * (-u.row(0).transpose().array()) - X.row(2).transpose().array() * vk(0) * (-u.row(0).transpose().array());
    A.col(6).tail(6) << Eigen::MatrixXd::Ones(6,1);
    A.col(7).head(6) << -Eigen::MatrixXd::Ones(6,1);
    A.col(8).head(6) << u.row(1).transpose().array();
    A.col(8).tail(6) << -u.row(0).transpose().array();
    A.col(9).tail(6) << u.row(0).transpose().array();
    A.col(10).head(6) << -u.row(0).transpose().array();
    A.col(11).head(6) << -u.row(1).transpose().array() * (-u.row(0).transpose().array());
    A.col(11).tail(6) << u.row(0).transpose().array() * (-u.row(0).transpose().array());
    A.col(12).head(6) << X.row(2).transpose().array() * u.row(1).transpose().array() - X.row(1).transpose().array();
    A.col(12).tail(6) << X.row(0).transpose().array() - X.row(2).transpose().array() * u.row(0).transpose().array();
    Eigen::MatrixXd n = A.fullPivLu().kernel();
    // std::cout << "A inverse done \n";
    int end = n.cols()-1;
    double s = n(12,end);
    // std::cout << "scalar done \n";
    Eigen::Vector3d v = n.col(end).head(3)/s; 
    // std::cout << "v: \n" << v << "\n";
    Eigen::Vector3d w = n.col(end).segment(3,3)/s;
    // std::cout << "w done \n" << w << "\n";
    Eigen::Vector3d C = n.col(end).segment(6,3)/s;
    // std::cout << "C done \n" << C << "\n";
    Eigen::Vector3d t = n.col(end).segment(9,3)/s;
    // std::cout << "t done \n" << t << "\n";

    results->push_back({v, C, w, t, 1, 0});

    return 0;

}