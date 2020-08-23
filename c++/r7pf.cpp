#include <iostream>
#include <chrono>
#include "RnP.h"
#include <iomanip>


template<typename Model, int (*Solver)(Eigen::Matrix<double,7,3> X, Eigen::Matrix<double,7,2> u, Eigen::Vector3d vk, double r0, std::vector<Model> * results)>
int iterativeRnP(Eigen::Matrix<double,7,3> X, Eigen::Matrix<double,7,2> u, Eigen::Vector3d vk, int sampleSize, double r0, int maxIter, Model & result){
    int notFound = 1;
    int k = 0;

    result.v = Eigen::Vector3d::Zero();
    result.C = Eigen::Vector3d::Zero();
    result.w = Eigen::Vector3d::Zero();
    result.t = Eigen::Vector3d::Zero();
    result.f = 1;
    result.rd = 0;

    double errPrev = 1e6;

    std::vector<Model> results;

    while(notFound && k < maxIter){
        Solver(X, u, result.v, 0,  &results); 
        // if the inner solver returned no solution
        if(!results.size()){
            return 0;
        }

        for(auto const& res: results){
            double errNew = calcErrAlgebraicRnPFocalRadialDoubleLin(res.v, res.C, res.w, res.t, res.f, res.rd,  r0,  X.transpose(),  u.transpose());
            if(errNew < errPrev){
                result = res;
                errPrev = errNew;
                std::cout << "errNew: " << errNew << "\n";
            }
        }

        k++;       
    }
    return 0;
}

int lin_w_t_v_C_focal_6(Eigen::Matrix<double,7,3> X, Eigen::Matrix<double,7,2> u, Eigen::Vector3d vk, double r0, RSDoublelinCameraPoseVector * results){

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
A.col(10) = X.col(1).cwiseProduct(u.col(0)) - X.col(0).cwiseProduct(u.col(1));

Eigen::MatrixXd nn = A.fullPivLu().kernel();

Eigen::MatrixXd n(10,4);


for ( int i = 0; i < 10; i++)
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
A1.col(5) = -u.col(1).array() * (X.col(1).array() * (n(3,3) * (r0 - u.col(0).array()) - n(0,3) + n(4,3) * vk(2) * (r0 - u.col(0).array())) + X.col(0).array() * (n(1,3) - n(4,3) * (r0 - u.col(0).array()) + n(3,3) * vk(2) * (r0-u.col(0).array())) - X.col(2).array() * (n(3,3) * vk(0) * (r0 - u.col(0).array()) + n(4,3) * vk(1) * (r0 - u.col(0).array())+1));
A1.col(3) = u.col(1);
A1.col(4) = -u.col(1).array() * (r0 - u.col(0).array());


Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
ges.compute(A0.topRows(6), A1.topRows(6),true);
Eigen::VectorXcd ff = ges.eigenvalues();
Eigen::MatrixXcd xx = ges.eigenvectors();

for (int i = 0; i < ff.rows(); i++)
{
    if(std::abs(ff(i).real()) > 10e-6 & std::abs(ff(i).imag()) < 10e-6){
        double f = 1/ff(i).real();
        double a, b, c;
        Eigen::Vector3d C,v,w,t;
        Eigen::VectorXd x = xx.col(i).real();
        a = x(0)/x(5);
        b = x(1)/x(5);
        c = x(2)/x(5);
        C(2) = x(3)/x(5);
        t(2) = x(4)/x(5);
        // TODO decrement indices
        v(0) = a*n(0,0)+b*n(0,1)+c*n(0,2)+n(0,3);
        v(1) = a*n(1,0)+b*n(1,1)+c*n(1,2)+n(1,3);
        v(2) = a*n(2,0)+b*n(2,1)+c*n(2,2)+n(2,3);
        w(0) = a*n(3,0)+b*n(3,1)+c*n(3,2)+n(3,3);
        w(1) = a*n(4,0)+b*n(4,1)+c*n(4,2)+n(4,3);
        w(2) = a*n(5,0)+b*n(5,1)+c*n(5,2)+n(5,3);
        C(0) = a*n(6,0)+b*n(6,1)+c*n(6,2)+n(6,3);
        C(1) = a*n(7,0)+b*n(7,1)+c*n(7,2)+n(7,3);
        t(0) = a*n(8,0)+b*n(8,1)+c*n(8,2)+n(8,3);
        t(1) = a*n(9,0)+b*n(9,1)+c*n(9,2)+n(9,3);
        std::cout << "C: " << C << "\n";
        results->push_back({v, C, w, t, f, 0});
    }
}


// std::cout << A0 << "\n";
return 0;

}

int main(int argc, char ** argv){

std::srand((unsigned int) std::time(0));

Eigen::MatrixXd X(3,7);
Eigen::MatrixXd u(2,7);

X << 0.537667139546100,   0.862173320368121,  -0.433592022305684,   2.769437029884877,   0.725404224946106,  -0.204966058299775,  1.409034489800479,
   1.833885014595086,   0.318765239858981,   0.342624466538650,  -1.349886940156521,  -0.063054873189656,  -0.124144348216312, 1.417192413429614,
  -2.258846861003648,  -1.307688296305273,   3.578396939725760,   3.034923466331855,   0.714742903826096,   1.489697607785465, 0.671497133608080;

u << -1.207486922685038,   1.630235289164729,   1.034693009917860,  -0.303440924786016,  -0.787282803758638,  -1.147070106969150,  -0.809498694424876,
   0.717238651328838,   0.488893770311789,   0.726885133383238,   0.293871467096658,   0.888395631757642,  -1.068870458168032,  -2.944284161994896;

// synthetic data generated with RS model

X << -0.503909660144734,  -0.298950855822539,   0.169495966151121,  -0.425456298098001,  -0.489857699307035,  -0.123874335620373,  -0.446421966521499,
   0.200846972671030,   0.593347927125033,  -0.482886654281585,   0.885086516289968,   0.806102240976620,   0.027467242888846,   0.501076519424143,
  -2.065389955548187,   0.054205758855624,  -0.425973877083998,   0.186843179728400,  -0.359508816247434,  -0.566837662670934,  -0.957820172224483;
  
u << -0.270836283773268,  -0.517671729904335,  -0.575222810376785,  -0.501436620776611,  -0.409834538142448,  -0.475797263728415,  -0.358882073270815,
  -0.050773801610145,  -0.222719459414455,   0.075702713083869,  -0.312858144184682,  -0.242367363513907,  -0.044286902417554,  -0.138074970066516;

RSDoublelinCameraPose result;

// int res =  iterativeRnP<RSDoublelinCameraPose,lin_w_t_v_C_focal_6>(Eigen::Matrix<double,7,3>::Random(), Eigen::Matrix<double,7,2>::Random(), Eigen::Vector3d::Random(), 7, 0.0, 5, result);

int res =  iterativeRnP<RSDoublelinCameraPose,lin_w_t_v_C_focal_6>(X.transpose(),u.transpose(), Eigen::Vector3d::Zero(), 7, 0.0, 5, result);


// std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
// for (size_t i = 0; i < 10000; i++)
// {
//     RSDoublelinCameraPose result;
//     int res =  iterativeRnP<RSDoublelinCameraPose,lin_w_t_v_C_focal_6>(Eigen::Matrix<double,7,3>::Random(), Eigen::Matrix<double,7,2>::Random(), Eigen::Vector3d::Random(), 7, 0.0, 1, result);
// }

// std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

// std::cout << "Time difference = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/10000) << "[µs]" << std::endl;
    

}