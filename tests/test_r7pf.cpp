#include "rnp.h"
#include "utils.h"
#include <math.h>
#include <iostream>
#include <random>

bool testR7Pf(){
    Eigen::MatrixXd X(3,7);
    Eigen::MatrixXd u(2,7);
    double r0 = 0;
    int direction = 0;
    int n_points = 7;
    int maxIter = 10;
    
    bool passed = true;

    std::default_random_engine random_engine;

    std::uniform_real_distribution<double> f_gen(1000, 1500);
    std::uniform_real_distribution<double> rd_gen(-0.5, 0);


    // // first test all zeros and f = 1 and X random
    RSDoublelinCameraPose gt;
    gt.v << 0,0,0;
    gt.w << 0,0,0;
    gt.C << 0,0,0;
    gt.t << 0,0,0;
    X = Eigen::MatrixXd::Random(3,7);
    gt.f = 1;
    gt.rd = 0;


    // create exact 2D projections
    for (int i = 0; i < n_points; i++)
    {
        Eigen::Vector2d temp;
        if(rs2LinProjection(X.col(i), temp, gt.v, gt.C, gt.w, gt.t,  gt.f,  gt.rd,  r0,  direction)){
            std::cout << "R7Pf test failed due to projection function\n";
            return 0;
        }
        u.col(i) = temp;
    }

    RSDoublelinCameraPose model;

    // run R7Pf
    int res =  iterativeRnP<RSDoublelinCameraPose, R7PfIter>(X, u, gt.v, n_points, r0, direction, maxIter, model);

    if(res == ERR_NO_SOLUTION){
        std::cout << "R7Pf for random X, direction 0, all params 0 and f = 1 returned no solution\n";
        passed = false;
    }else if(res == WARN_NO_CONVERGENCE){
        std::cout << "R7Pf for random X direction 0, all params 0 and f = 1 did not converge\n";
        passed = false;
    }else{
        if(!isPoseApproxEqual(gt,model)){
            std::cout << "R7Pf for random X direction 0, all params 0 and f = 1 returned bad pose\n";
        passed = false;
        }
    }

    // now test with random reasonable values

    gt.f = f_gen(random_engine);
    gt.C << Eigen::Vector3d::Random();
    gt.v << Eigen::Vector3d::Random()/10;
    gt.w << Eigen::Vector3d::Random()/gt.f/10;
    gt.t << Eigen::Vector3d::Random()/gt.f/10;
    gt.rd = 0;
    X << Eigen::Matrix<double,3,7>::Random() + Eigen::Vector3d(0,0,5).replicate(1,7);
    
    // create exact 2D projections
    for (int i = 0; i < n_points; i++)
    {
        Eigen::Vector2d temp;
        if(rs2LinProjection(X.col(i), temp, gt.v, gt.C, gt.w, gt.t,  gt.f,  gt.rd,  r0,  direction)){
            std::cout << "R7Pf test failed due to projection function\n";
            return 0;
        }
        u.col(i) = temp;
    }

    // run R7Pf
    res =  iterativeRnP<RSDoublelinCameraPose, R7PfIter>(X, u, gt.v, n_points, r0, direction, maxIter, model);

    if(res == ERR_NO_SOLUTION){
        std::cout << "R7Pf for random X, direction 0, random params and f = 1 returned no solution\n";
        passed = false;
    }else if(res == WARN_NO_CONVERGENCE){
        std::cout << "R7Pf for random X , direction 0, random params and f = 1 did not converge\n";
        passed = false;
    }else{
        if(!isPoseApproxEqual(gt,model)){
            std::cout << "R7Pf for random X , direction 0, random params and f = 1 returned bad pose\n";
        passed = false;
        }
    }

    // Now try other direction
    direction = 1;

    // create exact 2D projections
    for (int i = 0; i < n_points; i++)
    {
        Eigen::Vector2d temp;
        if(rs2LinProjection(X.col(i), temp, gt.v, gt.C, gt.w, gt.t,  gt.f,  gt.rd,  r0,  direction)){
            std::cout << "R7Pf test failed due to projection function\n";
            return 0;
        }
        u.col(i) = temp;
    }

    std::cout << "gt: \n";
    std::cout << "C: \n" << gt.C << "\n"; 
    std::cout << "t: \n" << gt.t << "\n"; 
    std::cout << "v: \n" << gt.v << "\n"; 
    std::cout << "w: \n" << gt.w << "\n"; 

    // run R7Pf
    res =  iterativeRnP<RSDoublelinCameraPose, R7PfIter>(X, u, gt.v, n_points, r0, direction, maxIter, model);

    std::cout << "res: \n";
    std::cout << "C: \n" << model.C << "\n"; 
    std::cout << "t: \n" << model.t << "\n"; 
    std::cout << "v: \n" << model.v << "\n"; 
    std::cout << "w: \n" << model.w << "\n"; 

    if(res == ERR_NO_SOLUTION){
        std::cout << "R7Pf for random X, direction 1, random params and f = 1 returned no solution\n";
        passed = false;
    }else if(res == WARN_NO_CONVERGENCE){
        std::cout << "R7Pf for random X , direction 1, random params and f = 1 did not converge\n";
        passed = false;
    }else{
        if(!isPoseApproxEqual(gt,model)){
            std::cout << "R7Pf for random X , direction 1, random params and f = 1 returned bad pose\n";
        passed = false;
        }
    }

    return passed;

}

int main(int argc, char ** argv){


if(testR7Pf()){
    std::cout << "R7Pf test passed\n";
}else{
    std::cerr << "R7Pf test did not pass\n";
    return 1;
};

return 0;
}