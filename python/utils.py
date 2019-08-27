import numpy as np

def x_(v):
    M = np.zeros((3,3))
    M[0,1] = -v[2]
    M[0,2] = v[1]
    M[1,0] = v[2]
    M[1,2] = -v[0]
    M[2,0] = -v[1]
    M[2,1] = v[0]
    return M

def h2a(v):
    return v[:-1,:]/v[[-1],:]

def a2h(v):
    temp = np.ones((v.shape[0]+1,v.shape[1]))
    temp[:-1,:] = v
    return temp