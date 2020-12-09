#include "rnp.h"
#include "utils.h"
#include <math.h>
#include "scene_generator.h"
#include <iostream>
#include <random>
#include <chrono>

template<typename T>
bool isPoseApproxEqual(T pose1, T pose2){
    double err = (pose1.C-pose2.C).norm();   
    err += (pose1.v-pose2.v).norm();
    err += (pose1.t-pose2.t).norm();
    err += (pose1.w-pose2.w).norm();
    err += std::abs(pose1.f-pose2.f)/std::max(pose1.f,1.0);
    err += std::abs(pose1.rd-pose2.rd);

    return err < 1e-6;
}


bool testRSDoubleLinProjection(){
    Eigen::Vector3d v,w,C,t,X;
    double r0 = 0;
    int direction = 0;
    double rd = 0;
    double f = 0;
    
    bool passed = true;

    v << 0,0,0;
    w << 0,0,0;
    C << 0,0,0;
    t << 0,0,0;
    X << 0,0,0;

    Eigen::Vector2d u;

    std::default_random_engine random_engine;

    int res;

    // test all zero input, should return NaN
    res = rsDoubleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
    if(!std::isnan(u(0)) || !std::isnan(u(1))){
        passed = false;
        std::cout << "Error: when passing all 0 should return NaN\n";
    }

    // now test with focal length 1 and non-zero X, should return a number
    f = 1;
    X << 1,1,1;
    res = rsDoubleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
    if(std::isnan(u(0)) || std::isnan(u(1))){
        passed = false;
        std::cout << "Error: when passing f = 1 should not return NaN\n";
    }

    // now test with random X and C 
    for (int i = 0; i < 1000; i++)
    {
        X = Eigen::Vector3d::Random();
        C = Eigen::Vector3d::Random();
        res = rsDoubleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
        if(std::isnan(u(0)) || std::isnan(u(1)) || res){
            passed = false;
            std::cout << "Error: projection failed for random X and C\n";
        }
    }

    // now test with random  v
    for (int i = 0; i < 1000; i++)
    {
        X << 0, 0, 2;
        C << 0,0,0;
        v = Eigen::Vector3d::Random();
        res = rsDoubleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
        if(std::isnan(u(0)) || std::isnan(u(1)) || res){
            passed = false;
            std::cout << "Error: projection failed for random v\n";
        }
    }

    // now test with random w
    for (int i = 0; i < 1000; i++)
    {
        X << 0, 0, 2;
        C << 0,0,0;
        v << 0,0,0;
        w = Eigen::Vector3d::Random();
        res = rsDoubleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
        if(std::isnan(u(0)) || std::isnan(u(1)) || res){
            passed = false;
            std::cout << "Error: projection failed for random w\n";
        }
    }

    // now test with random rd
    std::uniform_real_distribution<double> rd_gen(-0.5, 0);
    for (int i = 0; i < 1000; i++)
    {
        X << 0, 0, 2;
        C << 0,0,0;
        v << 0,0,0;
        w << 0,0,0;
        rd = rd_gen(random_engine);
        res = rsDoubleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
        if(std::isnan(u(0)) || std::isnan(u(1)) || res){
            passed = false;
            std::cout << "Error: projection failed for random rd\n";
        }
    }

    // now test reasonable values of X,C,t,v,w,f,rd
    std::uniform_real_distribution<double> f_gen(1000, 1500);
    for (int i = 0; i < 1000; i++)
    {
        f = f_gen(random_engine);
        X << Eigen::Vector3d::Random() + Eigen::Vector3d(0,0,5);
        C << Eigen::Vector3d::Random();
        v << Eigen::Vector3d::Random()/10;
        w << Eigen::Vector3d::Random()/f/10;
        t << Eigen::Vector3d::Random()/f/10;
        rd = rd_gen(random_engine)/f/f;
        res = rsDoubleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
        if(std::isnan(u(0)) || std::isnan(u(1)) || res){
            passed = false;
            std::cout << "Error: projection failed for random inputs\n";
        }
    }

    //now verify that the solution is accurate
    for (int i = 0; i < 1000; i++)
    {
        f = f_gen(random_engine);
        X << Eigen::Vector3d::Random() + Eigen::Vector3d(0,0,5);
        C << Eigen::Vector3d::Random();
        v << Eigen::Vector3d::Random()/10;
        w << Eigen::Vector3d::Random()/f/10;
        t << Eigen::Vector3d::Random()/f/10;
        rd = rd_gen(random_engine)/f/f;
        res = rsDoubleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
        double err =  calcErrAlgebraicRnPFocalRadialDoubleLin(v, C, w, t, f, rd, 0, X, u, direction);
        if(err > 1e-8){
            passed = false;
            std::cout << "Error: projection failed the accuracy test\n";
        }


    }

    //now verify that the solution is accurate for other RS direction
    direction = 1;
    for (int i = 0; i < 1000; i++)
    {
        f = f_gen(random_engine);
        X << Eigen::Vector3d::Random() + Eigen::Vector3d(0,0,5);
        C << Eigen::Vector3d::Random();
        v << Eigen::Vector3d::Random()/10;
        w << Eigen::Vector3d::Random()/f/10;
        t << Eigen::Vector3d::Random()/f/10;
        rd = rd_gen(random_engine)/f/f;
        res = rsDoubleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
        double err =  calcErrAlgebraicRnPFocalRadialDoubleLin(v, C, w, t, f, rd, 0, X, u, direction);
        if(err > 1e-8){
            passed = false;
            std::cout << "Error: projection failed the accuracy test for vertical RS direction\n";
        }


    }

    return passed;
}

bool testRSSingleLinProjection(){
    Eigen::Vector3d v,w,C,t,X;
    double r0 = 0;
    int direction = 0;
    double rd = 0;
    double f = 0;
    
    bool passed = true;

    v << 0,0,0;
    w << 0,0,0;
    C << 0,0,0;
    t << 0,0,0;
    X << 0,0,0;

    Eigen::Vector2d u;

    std::default_random_engine random_engine;

    int res;

    // test all zero input, should return NaN
    res = rsSingleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
    if(!std::isnan(u(0)) || !std::isnan(u(1))){
        passed = false;
        std::cout << "Error: when passing all 0 should return NaN\n";
    }

    // now test with focal length 1 and non-zero X, should return a number
    f = 1;
    X << 1,1,1;
    res = rsSingleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
    if(std::isnan(u(0)) || std::isnan(u(1))){
        passed = false;
        std::cout << "Error: when passing f = 1 should not return NaN\n";
    }

    // now test with random X and C 
    for (int i = 0; i < 1000; i++)
    {
        X = Eigen::Vector3d::Random();
        C = Eigen::Vector3d::Random();
        res = rsSingleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
        if(std::isnan(u(0)) || std::isnan(u(1)) || res){
            passed = false;
            std::cout << "Error: projection failed for random X and C\n";
        }
    }

    // now test with random  v
    for (int i = 0; i < 1000; i++)
    {
        X << 0, 0, 2;
        C << 0,0,0;
        v = Eigen::Vector3d::Random();
        res = rsSingleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
        if(std::isnan(u(0)) || std::isnan(u(1)) || res){
            passed = false;
            std::cout << "Error: projection failed for random v\n";
        }
    }

    // now test with random w
    for (int i = 0; i < 1000; i++)
    {
        X << 0, 0, 2;
        C << 0,0,0;
        v << 0,0,0;
        w = Eigen::Vector3d::Random();
        res = rsSingleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
        if(std::isnan(u(0)) || std::isnan(u(1)) || res){
            passed = false;
            std::cout << "Error: projection failed for random w\n";
        }
    }

    // now test with random rd
    std::uniform_real_distribution<double> rd_gen(-0.5, 0);
    for (int i = 0; i < 1000; i++)
    {
        X << 0, 0, 2;
        C << 0,0,0;
        v << 0,0,0;
        w << 0,0,0;
        rd = rd_gen(random_engine);
        res = rsSingleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
        if(std::isnan(u(0)) || std::isnan(u(1)) || res){
            passed = false;
            std::cout << "Error: projection failed for random rd\n";
        }
    }

    // now test reasonable values of X,C,t,v,w,f,rd
    std::uniform_real_distribution<double> f_gen(1000, 1500);
    for (int i = 0; i < 1000; i++)
    {
        f = f_gen(random_engine);
        X << Eigen::Vector3d::Random() + Eigen::Vector3d(0,0,5);
        C << Eigen::Vector3d::Random();
        v << Eigen::Vector3d::Random()/10;
        w << Eigen::Vector3d::Random()/f/10;
        t << Eigen::Vector3d::Random()/f/10;
        rd = rd_gen(random_engine)/f/f;
        res = rsSingleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
        if(std::isnan(u(0)) || std::isnan(u(1)) || res){
            passed = false;
            std::cout << "Error: projection failed for random inputs\n";
        }
    }

    //now verify that the solution is accurate
    for (int i = 0; i < 1000; i++)
    {
        f = f_gen(random_engine);
        X << Eigen::Vector3d::Random() + Eigen::Vector3d(0,0,5);
        C << Eigen::Vector3d::Random();
        v << Eigen::Vector3d::Random()/10;
        w << Eigen::Vector3d::Random()/f/10;
        t << Eigen::Vector3d::Random()/f/10;
        rd = rd_gen(random_engine)/f/f;
        res = rsSingleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
        double err =  calcErrAlgebraicRnPFocalRadialSingleLin(v, C, w, t, f, rd, 0, X, u, direction);
        if(err > 1e-8){
            passed = false;
            std::cout << "Error: projection failed the accuracy test\n";
        }


    }

    //now verify that the solution is accurate for other RS direction
    direction = 1;
    for (int i = 0; i < 1000; i++)
    {
        f = f_gen(random_engine);
        X << Eigen::Vector3d::Random() + Eigen::Vector3d(0,0,5);
        C << Eigen::Vector3d::Random();
        v << Eigen::Vector3d::Random()/10;
        w << Eigen::Vector3d::Random()/f/10;
        t << Eigen::Vector3d::Random()/f/10;
        rd = rd_gen(random_engine)/f/f;
        res = rsSingleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
        double err =  calcErrAlgebraicRnPFocalRadialSingleLin(v, C, w, t, f, rd, 0, X, u, direction);
        if(err > 1e-8){
            passed = false;
            std::cout << "Error: projection failed the accuracy test for vertical RS direction\n";
        }


    }

    return passed;
}

bool testR7Pf(){
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


    // // first test all zeros and f = 1 and X random
    RSDoublelinCameraPose gt;
    gt.v << 0,0,0;
    gt.w << 0,0,0;
    gt.C << 0,0,0;
    gt.t << 0,0,0;
    X = Eigen::Matrix<double,3,7>::Random();
    gt.f = 1;
    gt.rd = 0;


    // create exact 2D projections
    for (int i = 0; i < n_points; i++)
    {
        Eigen::Vector2d temp;
        if(rsDoubleLinProjection(X.col(i), temp, gt.v, gt.C, gt.w, gt.t,  gt.f,  gt.rd,  r0,  direction)){
            std::cout << "R7Pf test failed due to projection function\n";
            return 0;
        }
        u.col(i) = temp;
    }

    RSDoublelinCameraPose model;

    // run R7Pf
    int res =  iterativeRnP<RSDoublelinCameraPose, R7PfLin>(X.transpose(), u.transpose(), gt.v, n_points, r0, maxIter, model);

    if(res == ERR_NO_SOLUTION){
        std::cout << "R7Pf for random X all params 0 and f = 1 returned no solution\n";
        passed = false;
    }else if(res == WARN_NO_CONVERGENCE){
        std::cout << "R7Pf for random X all params 0 and f = 1 did not converge\n";
        passed = false;
    }else{
        if(!isPoseApproxEqual(gt,model)){
            std::cout << "R7Pf for random X all params 0 and f = 1 returned bad pose\n";
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
        if(rsDoubleLinProjection(X.col(i), temp, gt.v, gt.C, gt.w, gt.t,  gt.f,  gt.rd,  r0,  direction)){
            std::cout << "R7Pf test failed due to projection function\n";
            return 0;
        }
        u.col(i) = temp;
    }

    // run R7Pf
    res =  iterativeRnP<RSDoublelinCameraPose, R7PfLin>(X.transpose(), u.transpose(), gt.v, n_points, r0, maxIter, model);

    if(res == ERR_NO_SOLUTION){
        std::cout << "R7Pf for random X, random params and f = 1 returned no solution\n";
        passed = false;
    }else if(res == WARN_NO_CONVERGENCE){
        std::cout << "R7Pf for random X , random params and f = 1 did not converge\n";
        passed = false;
    }else{
        if(!isPoseApproxEqual(gt,model)){
            std::cout << "R7Pf for random X , random params and f = 1 returned bad pose\n";
        passed = false;
        }
    }

    return passed;

}

bool testR7Pfr(){
    Eigen::Matrix<double,3,7> X;
    Eigen::Matrix<double,2,7> u;
    double r0 = 0;
    int direction = 0;
    int n_points = 7;
    int maxIter = 10;
    
    bool passed = true;

    typedef std::chrono::high_resolution_clock myclock;

    std::default_random_engine random_engine(myclock::now().time_since_epoch().count());

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

if(testRSDoubleLinProjection()){
    std::cout << "rsDoubleLinProjection test passed\n";
}else{
    std::cerr << "rsDoubleLinProjection test did not pass\n";
};

if(testRSSingleLinProjection()){
    std::cout << "rsSingleLinProjection test passed\n";
}else{
    std::cerr << "rsSingleLinProjection test did not pass\n";
};

if(testR7Pf()){
    std::cout << "R7Pf test passed\n";
}else{
    std::cerr << "R7Pf test did not pass\n";
};

if(testR7Pfr()){
    std::cout << "R7Pfr test passed\n";
}else{
    std::cerr << "R7Pfr test did not pass\n";
};

}