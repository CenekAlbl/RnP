import numpy as np
import rnp
import ransac 


# test a single run of R6PLin on perfect perspective data
X = np.array([[-0.735203763665459, -0.496627573213431, -0.872484275284515, -0.427093140715349, -0.834433230284650,  -0.510357709134182], [0.238783627918357, -0.331469664458513, -0.774067735937978, -0.330087062598541, -0.532562868622104, -0.175836814526873], [0.348232082844991, 0.496786171364075, 0.842812754497372, 0.900023147138860, 0.062454513759772,0.920453051685240]])
u = np.array([[0.087573419907225, 0.234591535755890, -0.008470881920492,0.221656649106303, 0.020112962609198, 0.176817312337601], [0.377812297395972, -0.030443875680965, -0.260227638405494, -0.023394282943956, -0.225285886765297,  0.055631494170819]])

data = np.vstack((X,u))

res = rnp.r6pLin(data)
print("Camera center is:")
print(res[0][9:12]) 
print("Camera orientation in angle-axis representation is:")
print(res[0][6:9]) 
print("Camera rotation velocity is:")
print(res[0][0:3]) 
print("Camera translation velocity is:")
print(res[0][3:6]) 

# test R6PLin in LO-RANSAC on perfect perspective data
threshold = 1
maxiter = 10
result = ransac.loRansacSimple(rnp.r6pLin, rnp.calcR6PEAXErrFast, data, 6, threshold, 10, optimizeFn = rnp.calcR6PEAXErrFastForLsq, verbose = 1)
        