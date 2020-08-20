#include "rnp.h"

template<typename Model, void (*Solver)(Eigen::Matrix<double,7,3> X, Eigen::Matrix<double,7,2> u, Eigen::Vector3d vk, double r0, std::vector<Model> * results)>
int iterativeRnP(Eigen::Matrix<double,7,3> X, Eigen::Matrix<double,7,2> u, Eigen::Vector3d vk, int sampleSize, double r0, int maxIter, Model * result){
    int notFound = 1;
    int k = 0;

    std::vector<Model> results;

    result.v = Eigen::Vector3d::Zeros();
    result.C = Eigen::Vector3d::Zeros();
    result.w = Eigen::Vector3d::Zeros();
    result.t = Eigen::Vector3d::Zeros();
    result.f = 1;
    result.rd = 0;

    double errPrev = 1e6;

    while(notFound && k < maxIter){
        Solver(X, u, vk, &results); 
        // if the inner solver returned no solution
        if(!results.size()){
            return 0;
        }

        double errPrev = 1e6;

        for(auto const& res: results){
            std::cout << res.C(0) << "\n";
        }

        k++;       
        
    }
    return 0;
}