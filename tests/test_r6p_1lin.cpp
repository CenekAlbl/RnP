#include "rnp.h"
#include "utils.h"
#include <math.h>
#include <iostream>
#include <random>

bool testR6PSingleLin(){
    std::srand(2);
    Eigen::MatrixXd X(3,6);
    Eigen::MatrixXd u(2,6);
    
    double r0 = 0;
    int direction = 0;
    int n_points = 6;
    int maxpow = 2;
    
    bool passed = true;



    // // first test all zeros and f = 1 and X random
    RSSinglelinCameraPose gt;
    gt.v << 0,0,0;
    gt.w << 0,0,0;
    gt.C << 0,0,0;
    gt.t << 0,0,0;
    X = Eigen::Matrix<double,3,6>::Random();
    gt.f = 1;
    gt.rd = 0;


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

    RSSinglelinCameraPoseVector results;

    int res = R6P1Lin(X, u, direction, r0, maxpow, &results);
    if(res == ERR_NO_SOLUTION){
        std::cout << "R6P singlelin, RS in x direction, for random X all params 0 and f = 1 returned no solution\n";
        passed = false;
    }else if(res == WARN_NO_CONVERGENCE){
        std::cout << "R6P singlelin, RS in x direction, for random X all params 0 and f = 1 did not converge\n";
        passed = false;
    }else{
        bool ident = false;
        for (auto const& model: results)
        {
            if(isPoseApproxEqual(gt,model)){
                ident = true;
            }
        }
        if(!ident){
            std::cout << "R6P singlelin, RS in x direction, for random X all params 0 and f = 1 returned bad pose\n";
            passed = false;
        }
    }

    // now test with random reasonable values

    gt.f = 1;
    gt.C << Eigen::Vector3d::Random();
    gt.v << Eigen::Vector3d::Random()/10;
    gt.w << Eigen::Vector3d::Random()/gt.f/10;
    gt.t << Eigen::Vector3d::Random()/gt.f/10;
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


    res = R6P1Lin(X, u, direction, r0, maxpow, &results);

    if(res == ERR_NO_SOLUTION){
        std::cout << "R6P singlelin, RS in x direction, for random values returned no solution\n";
        passed = false;
    }else if(res == WARN_NO_CONVERGENCE){
        std::cout << "R6P singlelin, RS in x direction, for random values did not converge\n";
        passed = false;
    }else{
        bool ident = false;
        for (auto const& model: results)
        {
            if(isPoseApproxEqual(gt,model)){
                ident = true;
            }
        }
        if(!ident){
            std::cout << "R6P singlelin, RS in x direction, for random values returned bad pose\n";
            passed = false;
        }
    }

    
    

    // Now test with the other RS direction
    direction = 1;
    // create exact 2D projections
    for (int i = 0; i < n_points; i++)
    {
        Eigen::Vector2d temp;
        if(rs1LinProjection(X.col(i), temp, gt.v, gt.C, gt.w, gt.t,  gt.f,  gt.rd,  r0,  direction)){
            std::cout << "R6P singlelin, RS in y direction, test failed due to projection function\n";
            return 0;
        }
        u.col(i) = temp;
    }

    res = R6P1Lin(X, u, direction, r0, maxpow, &results);

    if(res == ERR_NO_SOLUTION){
        std::cout << "R6P singlelin, RS in y direction, for random values returned no solution\n";
        passed = false;
    }else if(res == WARN_NO_CONVERGENCE){
        std::cout << "R6P singlelin, RS in y direction, for random values did not converge\n";
        passed = false;
    }else{
        bool ident = false;
        for (auto const& model: results)
        {
            if(isPoseApproxEqual(gt,model)){
                ident = true;
            }
        }
        if(!ident){
            std::cout << "R6P singlelin, RS in y direction, for random values returned bad pose\n";
            passed = false;
        }
    }

    return passed;

}

int main(int argc, char ** argv){
    if(testR6PSingleLin()){
        std::cout << "R6P singlelin test passed\n";
    }else{
        std::cerr << "R6P singlelin test did not pass\n";
        return 1;
    };

    return 0;
}