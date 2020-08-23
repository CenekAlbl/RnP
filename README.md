# RnP
Algorithms for computing the rolling shutter absolute camera pose. They are available in various languages (C++,Python,Matlab), but sometimes not in all of them, see below.

![Rolling shutter absolute pose][absolute_pose_rs.png]

*This repository is in the process of being completed, some parts still missing*

## R7P
This is the 7 point version for cameras with unknown intrinsics. There are 2 algorithms:
* R7Pf - computes the camera pose + velocity + the focal length
* R7Pfr - computes the camera pose + velocity + the focal length and radial distortion  

These algorithms are currently only available as research code in Matlab. C++ and Python implementation plus documentation and tests are underway. If you are interested in efficient implementation right now, please contact me and I can make it bigger priority.

## R6P-1lin
This is the single-linearized standalone RS absolute pose solver presented in:

C. Albl, Z. Kukelova, V. Larsson and T. Pajdla, "Rolling Shutter Camera Absolute Pose," in *IEEE Transactions on Pattern Analysis and Machine Intelligence*, 2019.

## R6P-2lin
This is the double-linearized RS absolute pose solver that needs an initial camera orientation estimate also presented in:

C. Albl, Z. Kukelova, V. Larsson and T. Pajdla, "Rolling Shutter Camera Absolute Pose," in *IEEE Transactions on Pattern Analysis and Machine Intelligence*, 2019.

which solves the same equations, but faster than the original:

C. Albl, Z. Kukelova and T. Pajdla, "R6P - Rolling shutter absolute pose problem," *2015 IEEE Conference on Computer Vision and Pattern Recognition (CVPR)*, 2015.

## R6P-2lin-iter
This is the iterative version of R6P-2lin that needs an initial camera orientation estimate and was presented in:

Kukelova Z., Albl C., Sugimoto A., Pajdla T. Linear Solution to the Minimal Absolute Pose Rolling Shutter Problem. *Asian Conference on Computer Vision*, 2018.








