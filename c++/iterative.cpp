#include "rnp.h"
#include <iostream>

template<typename Model, int (*Solver)(const Eigen::Matrix<double,3,7> &X, const Eigen::Matrix<double,2,7> &u, const Eigen::Vector3d &vk, double r0, std::vector<Model> * results)>
int iterativeRnP(Eigen::Matrix<double,3,7> &X, Eigen::Matrix<double,2,7> &u, Eigen::Vector3d vk, int sampleSize, double r0, int direction, int maxIter, Model & result){
    int notFound = 1;
    int k = 0;

    if(direction){
        inputSwitchDirection(X, u);
    }

    result.v = Eigen::Vector3d::Zero();
    result.C = Eigen::Vector3d::Zero();
    result.w = Eigen::Vector3d::Zero();
    result.t = Eigen::Vector3d::Zero();
    result.f = 1;
    result.rd = 0;

    while(notFound && k < maxIter){
        std::vector<Model> results;
        double errPrev = 1e15;
        Solver(X, u, result.v, 0,  &results); 
        // if the inner solver returned no solution
        if(!results.size()){
            std::cout << "Solver returned no solution\n";
            return 1;
        }

        for(auto const& res: results){
            double errNew = calcErrAlgebraicRnPFocalRadialDoubleLin(res.v, res.C, res.w, res.t, res.f, res.rd,  r0,  X,  u, 0);
            if(errNew < errPrev){
                result = res;
                errPrev = errNew;
            }
            if(errNew < 1e-10){
                notFound = 0;
            }
        }
        k++;       
    }

    if(direction){
        result = outputSwitchDirection(result);
    }

    return 0;
}


template int iterativeRnP<RSDoublelinCameraPose,R7PfrLin>(Eigen::Matrix<double,3,7> &X, Eigen::Matrix<double,2,7> &u, Eigen::Vector3d vk, int sampleSize, double r0, int direction, int maxIter, RSDoublelinCameraPose & result);
template int iterativeRnP<RSDoublelinCameraPose,R7PfLin>(Eigen::Matrix<double,3,7> &X, Eigen::Matrix<double,2,7> &u, Eigen::Vector3d vk, int sampleSize, double r0, int direction, int maxIter, RSDoublelinCameraPose & result);