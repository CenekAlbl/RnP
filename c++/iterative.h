#include "rnp.h"

template<typename Model, void (*Solver)(Eigen::Matrix<double,7,3> X, Eigen::Matrix<double,7,2> u, Eigen::Vector3d vk, double r0, std::vector<Model> * results)>
int iterativeRnP(Eigen::Matrix<double,7,3> X, Eigen::Matrix<double,7,2> u, Eigen::Vector3d vk, int sampleSize, double r0, int maxIter){
    int notFound = 1;
    int k = 0;

    std::vector<Model> results;

    while(notFound && k < maxIter){
        Model model = Solver(X, u, vk, &results); 
        k++;       
        std::cout << results[0].C(0) << "\n";
    }
    return 0;
}