#include "rnp.h"
#include "utils.h"
#include <math.h>
#include <iostream>
#include <random>

bool testR7Pfr(){
    Eigen::Matrix<double,3,7> X;
    Eigen::Matrix<double,2,7> u;
    double r0 = 0;
    int direction = 0;
    int n_points = 7;
    int maxIter = 10;
    
    bool passed = true;


    std::default_random_engine random_engine;

    std::uniform_real_distribution<double> f_gen(1000, 1500);
    std::uniform_real_distribution<double> rd_gen(-0.5, 0);


    //  test with random reasonable values

    RSDoublelinCameraPose gt, model;

    gt.f = f_gen(random_engine);
    gt.C << Eigen::Vector3d::Random();
    gt.v << Eigen::Vector3d::Random()/10;
    gt.w << Eigen::Vector3d::Random()/gt.f/10;
    gt.t << Eigen::Vector3d::Random()/gt.f/10;
    gt.rd = rd_gen(random_engine)/gt.f/gt.f;
    X << Eigen::Matrix<double,3,7>::Random() + Eigen::Vector3d(0,0,5).replicate(1,7);
    
    // create exact 2D projections
    for (int i = 0; i < n_points; i++)
    {
        Eigen::Vector2d temp;
        if(rsDoubleLinProjection(X.col(i), temp, gt.v, gt.C, gt.w, gt.t,  gt.f,  gt.rd,  r0,  direction)){
            std::cout << "R7Pfr test failed due to projection function\n";
            return 0;
        }
        u.col(i) = temp;
    }

    // run R7Pf
    int res =  iterativeRnP<RSDoublelinCameraPose, R7PfrLin>(X.transpose(), u.transpose(), gt.v, n_points, r0, maxIter, model);

    

    if(res == ERR_NO_SOLUTION){
        std::cout << "R7Pfr returned no solution\n";
        passed = false;
    }else if(res == WARN_NO_CONVERGENCE){
        std::cout << "R7Pfr  did not converge\n";
        passed = false;
    }else{
        if(!isPoseApproxEqual(gt,model)){
            std::cout << "R7Pfr returned bad pose\n";
        passed = false;
        }
    }

    return passed;

}

int main(int argc, char ** argv){


if(testR7Pfr()){
    std::cout << "R7Pfr test passed\n";
}else{
    std::cerr << "R7Pfr test did not pass\n";
    return 1;
};

return 0;
}