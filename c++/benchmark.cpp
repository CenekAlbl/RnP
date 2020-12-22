#include "rnp.h"
#include <chrono>
#include <iostream>

int main(int argc, char ** argv){
    Eigen::MatrixXd X(3,6);
    Eigen::MatrixXd u(2,6);
    
    double r0 = 0;
    int direction = 0;
    int n_points = 6;
    int maxpow = 2;
    int niters = 100;

    RSSinglelinCameraPose gt;

    gt.f = 1;
    gt.C << Eigen::Vector3d::Random();
    gt.v << Eigen::Vector3d::Random()/10;
    gt.w << Eigen::Vector3d::Random()/10;
    gt.t << Eigen::Vector3d::Random()/10;
    gt.rd = 0;
    X << Eigen::MatrixXd::Random(3,6) + Eigen::Vector3d(0,0,5).replicate(1,6);

    // create exact 2D projections
    for (int i = 0; i < n_points; i++)
    {
        Eigen::Vector2d temp;
        if(rs1LinProjection(X.col(i), temp, gt.v, gt.C, gt.w, gt.t,  gt.f,  gt.rd,  r0,  direction)){
            std::cout << "R6P singlelin, RS in x direction, test failed due to projection function\n";
            return 0;
        }
        u.col(i) = temp;
    }

    // BENCHMARK R6P1Lin

    RSSinglelinCameraPoseVector results1Lin;
    
    auto t1 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i< niters; i++){
        int res = R6P1Lin(X, u, direction, r0, maxpow, &results1Lin);
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count()/niters;

    std::cout << "R6P1Lin runtime is: " << duration << " microseconds \n";

    // BENCHMARK R6P2Lin

    RSDoublelinCameraPoseVector results2Lin;
    
    t1 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i< niters; i++){
        int res = R6P2Lin(X, u, direction, r0, &results2Lin);
    }
    t2 = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count()/niters;

    std::cout << "R6P2Lin runtime is: " << duration << " microseconds \n";


    // BENCHMARK R6PIter

    int maxIter = 5;

    RSDoublelinCameraPose model;
    
    t1 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i< niters; i++){
        int res =  iterativeRnP<RSDoublelinCameraPose, R6PIter>(X, u, gt.v, n_points, r0, direction, maxIter, model);
    }
    t2 = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count()/niters;

    std::cout << "R6PIter (5 iterations) runtime is: " << duration << " microseconds \n";

    // BENCHMARK R7Pfr

    X = Eigen::MatrixXd::Zero(3,7);
    u = Eigen::MatrixXd::Zero(2,7);
    
    n_points = 7;


    gt.f = 1500;
    gt.C << Eigen::Vector3d::Random();
    gt.v << Eigen::Vector3d::Random()/10;
    gt.w << Eigen::Vector3d::Random()/gt.f/10;
    gt.t << Eigen::Vector3d::Random()/gt.f/10;
    gt.rd = -1e-7;
    X << Eigen::MatrixXd::Random(3,7) + Eigen::Vector3d(0,0,5).replicate(1,7);

    // create exact 2D projections
    for (int i = 0; i < n_points; i++)
    {
        Eigen::Vector2d temp;
        if(rs1LinProjection(X.col(i), temp, gt.v, gt.C, gt.w, gt.t,  gt.f,  gt.rd,  r0,  direction)){
            std::cout << "R6P singlelin, RS in x direction, test failed due to projection function\n";
            return 0;
        }
        u.col(i) = temp;
    }

    
    t1 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i< niters; i++){
        int res =  iterativeRnP<RSDoublelinCameraPose, R7PfIter>(X, u, gt.v, n_points, r0, direction, maxIter, model);
    }
    t2 = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count()/niters;

    std::cout << "R7Pf (5 iterations) runtime is: " << duration << " microseconds \n";

    // BENCHMARK R7Pfr

    t1 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i< niters; i++){
        int res =  iterativeRnP<RSDoublelinCameraPose, R7PfrIter>(X, u, gt.v, n_points, r0, direction, maxIter, model);
    }
    t2 = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count()/niters;

    std::cout << "R7Pfr (5 iterations) runtime is: " << duration << " microseconds \n";


}