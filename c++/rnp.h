#include <Eigen/Eigen>
#include <vector>

struct RSCameraPose{
    Eigen::Matrix3d R; // Camera orientation
    Eigen::Vector3d C; // Camera center
    Eigen::Vector3d w; // Camera rotational velocity (rotation axis where ||w|| is the angular velocity)
    Eigen::Vector3d t; // Camera translational velocity
    double f; // focal length (only R7Pf and R7Pfr compute focal lenght, otherwise 1)
    double rd; // radial distortion coefficient (only R7Pfr computes it, otherwise 0)
};

typedef std::vector<RSCameraPose> RSCameraPoseVector;



int r6pSingleLin(double * X, double * u, int direction, double r0, int maxpow, RSCameraPoseVector * output);

int r6pDoubleLin(double * X, double * u, int direction, double r0, RSCameraPoseVector * output);

