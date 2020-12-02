#include "rnp.h"
#include "utils.h"
#include <math.h>
#include "scene_generator.h"
#include <iostream>
#include <random>


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

    // test all zero input, should return NaN
    bool res = rsDoubleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
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
        if(std::isnan(u(0)) || std::isnan(u(1)) || !res){
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
        if(std::isnan(u(0)) || std::isnan(u(1)) || !res){
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
        if(std::isnan(u(0)) || std::isnan(u(1)) || !res){
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
        if(std::isnan(u(0)) || std::isnan(u(1)) || !res){
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
        if(std::isnan(u(0)) || std::isnan(u(1)) || !res){
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

    // test all zero input, should return NaN
    bool res = rsSingleLinProjection(X, u, v, C, w, t,  f,  rd,  r0,  direction);
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
        if(std::isnan(u(0)) || std::isnan(u(1)) || !res){
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
        if(std::isnan(u(0)) || std::isnan(u(1)) || !res){
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
        if(std::isnan(u(0)) || std::isnan(u(1)) || !res){
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
        if(std::isnan(u(0)) || std::isnan(u(1)) || !res){
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
        if(std::isnan(u(0)) || std::isnan(u(1)) || !res){
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

}