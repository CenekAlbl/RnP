# RnP
Algorithms for computing the rolling shutter absolute camera pose. 

![Rolling shutter absolute pose](absolute_pose_rs.png)

More specifically:
* **R6P** - computes calibrated RS camera pose and velocity from 6 3D<->2D correspondences
* **R7Pf** - computes RS camera pose and focal length from 7 3D<->2D correspondences
* **R7Pfr** - computes RS camera pose, focal length and radial distortion from 7 3D<->2D correspondences

For more details please see below.

The algorithms are available in the following languages:

* **C++** (all algorithms)
* **Python** (all algorithms via C++ bindings, R6PIter also in pure Python)
* **Matlab** (R6P1Lin via MEX, R7Pf and R7Pfr also in pure Matlab)

## Installation
The library can be compiled with CMake. This will produce the static C++ library and if Python and/or Matlab is installed on the machine it will try to compile the respective bindings. The Python bindings are produced thanks to the amazing library [Pybind](https://github.com/pybind/pybind11) which is included as a submodule. 

Follow the classic procedure:  

    git clone https://github.com/CenekAlbl/RnP.git
    cd RnP
    git submodule update --init --recursive
    mkdir build
    cd build
    cmake ..
    make

Since the compilation with optimizations can take very long, we recommend trying DEBUG version first:

    cmake .. -DCMAKE_BUILD_TYPE=DEBUG

After a successful compilation, we recommend to run the tests:

    make test

You can also run the benchmark to see the runtime of each algorithm on your machine:

    ./benchmark

Naturally, you should do it after compilation in the default RELEASE mode, not DEBUG.

## Usage

For an example how to use the algorithms take a look in the folder tests or in c++/benchmark.cpp.

## R7P
This is the 7 point version for cameras with unknown intrinsics. There are 2 algorithms:
* R7Pf - computes the camera pose + velocity + the focal length
* R7Pfr - computes the camera pose + velocity + the focal length and radial distortion  

The algorithms were presented in:  

Z. Kukelova, C. Albl, A. Sugimoto, K. Schindler, T. Pajdla, "Minimal Rolling Shutter Absolute Pose with Unknown Focal Length and Radial", *European Conference on Computer Vision (ECCV)* 2020

## R6P1lin
This is the single-linearized standalone RS absolute pose solver presented in:

C. Albl, Z. Kukelova, V. Larsson and T. Pajdla, "Rolling Shutter Camera Absolute Pose," in *IEEE Transactions on Pattern Analysis and Machine Intelligence*, 2019.

## R6P2lin
This is the double-linearized RS absolute pose solver that needs an initial camera orientation estimate also presented in:

C. Albl, Z. Kukelova, V. Larsson and T. Pajdla, "Rolling Shutter Camera Absolute Pose," in *IEEE Transactions on Pattern Analysis and Machine Intelligence*, 2019.

which solves the same equations, but faster than the original:

C. Albl, Z. Kukelova and T. Pajdla, "R6P - Rolling shutter absolute pose problem," *2015 IEEE Conference on Computer Vision and Pattern Recognition (CVPR)*, 2015.

## R6PIter
This is the iterative version of R6P-2lin that needs an initial camera orientation estimate and was presented in:

Kukelova Z., Albl C., Sugimoto A., Pajdla T. Linear Solution to the Minimal Absolute Pose Rolling Shutter Problem. *Asian Conference on Computer Vision*, 2018.








