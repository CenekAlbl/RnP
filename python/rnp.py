import numpy as np
import math as m
import utils as ut
from scipy.linalg import null_space


# angle-axis to rotation matrix
def eax2R(t):
    """Converts angle-axis vector to rotation matrix"""
    theta = np.linalg.norm(t)
    if theta < np.finfo(float).eps:   # If the rotation is very small...
        return np.array((
            (1 , -t[2, 0], t[1, 0]),
            (t[2, 0], 1, -t[0, 0]),
            (-t[1, 0], t[0, 0], 1))
            )
    
    # Otherwise set up standard matrix, first setting up some convenience
    # variables
    t = t/theta;  x = t[0, 0]; y = t[1, 0]; z = t[2, 0]
    
    c = m.cos(theta); s = m.sin(theta); C = 1-c
    xs = x*s;   ys = y*s;   zs = z*s
    xC = x*C;   yC = y*C;   zC = z*C
    xyC = x*yC; yzC = y*zC; zxC = z*xC

    return np.array((
        (x*xC+c,xyC-zs,zxC+ys),
        (xyC+zs, y*yC+c, yzC-xs),
        (zxC-ys, yzC+xs, z*zC+c))
        )

# angle-axis to rotation matrix fast version for RS projection
def eax2RforRS(axis,angles):
    if(len(angles.shape)>1):
        n = angles.shape[1]
    else:
        n = angles.shape[0]
    theta = np.linalg.norm(axis)
    if theta < np.finfo(float).eps:   # If the rotation is very small...
        return np.tile(np.array((
        (1 , -axis[2, 0], axis[1, 0]),
        (axis[2, 0], 1, -axis[0, 0]),
        (-axis[1, 0], axis[0, 0], 1))
        ), (1,n))

    # Otherwise set up standard matrix, first setting up some convenience
    # variables
    axis = axis/theta
    
    vx = np.tile(ut.x_(axis),(n,1))
    vx2 = np.tile(ut.x_(axis).dot(ut.x_(axis)),(n,1))
    s = np.tile(np.sin(angles),(3,1)).transpose()
    c = np.tile(np.cos(angles),(3,1)).transpose()
    sM = np.tile(s.ravel().reshape(n*3,1),(1,3))
    cM = np.tile(c.ravel().reshape(n*3,1),(1,3))
    I = np.tile(np.eye(3),(n,1))

    Rw = I + vx*sM + (1-cM)*vx2
    return Rw

def calcR6PEAXErrFast(M,data):
        n = data.shape[1]
        w = M[0:3].reshape(3,1)
        t = M[3:6].reshape(3,1)
        v = M[6:9].reshape(3,1)
        C = M[9:12].reshape(3,1)
        u = data[3:5,:]
        X = data[0:3,:]
        col = u[0,:]
        r = np.tile(col,(3,1))
        Rw = eax2RforRS(w,col*np.linalg.norm(w))
        Rv = eax2R(v)
        Rt = Rw.dot(Rv)
        Xrep = np.tile(X,(3,1))
        RX = np.sum(Rt.transpose()*Xrep.ravel(order='F').reshape((3,3*n),order='F'),axis=0).reshape((3,n),order='F')
        Z = RX + np.tile(C,(1,n)) + r*np.tile(t,(1,n))
        u_rs = ut.h2a(Z)
        err  = np.sqrt(np.sum((u_rs-u)**2,axis=0))
        return err

def calcR6PEAXErrFastForLsq(M,data):
    n = data.shape[1]
    w = M[0:3].reshape(3,1)
    t = M[3:6].reshape(3,1)
    v = M[6:9].reshape(3,1)
    C = M[9:12].reshape(3,1)
    u = data[3:5,:]
    X = data[0:3,:]
    col = u[0,:]
    r = np.tile(col,(3,1))
    Rw = eax2RforRS(w,col*np.linalg.norm(w))
    Rv = eax2R(v)
    Rt = Rw.dot(Rv)
    Xrep = np.tile(X,(3,1))
    RX = np.sum(Rt.transpose()*Xrep.ravel(order='F').reshape((3,3*n),order='F'),axis=0).reshape((3,n),order='F')
    Z = RX + np.tile(C,(1,n)) + r*np.tile(t,(1,n))
    u_rs = ut.h2a(Z)
    err  = (u_rs-u).ravel()
    return err

def linwtvC(X,u,vk):    
    A = np.zeros((12,13))
    A[:,[0]] = np.vstack((X[:,[2]] + X[:,[1]]*u[:,[1]], -X[:,[1]]*u[:,[0]]))
    A[:,[1]] = np.vstack((-X[:,[0]]*u[:,[1]], X[:,[2]] + X[:,[0]]*u[:,[0]]))
    A[:,[2]] = np.vstack((-X[:,[0]], -X[:,[1]]))
    A[:,[3]] = np.vstack((X[:,[0]]*vk[1]*( - u[:,[0]]) - X[:,[2]]*( - u[:,[0]]) - u[:,[1]]*(X[:,[1]]*( - u[:,[0]]) + X[:,[0]]*vk[2]*( - u[:,[0]]) - X[:,[2]]*vk[0]*( - u[:,[0]])) - X[:,[1]]*vk[0]*( - u[:,[0]]), u[:,[0]]*(X[:,[1]]*( - u[:,[0]]) + X[:,[0]]*vk[2]*( - u[:,[0]]) - X[:,[2]]*vk[0]*( - u[:,[0]]))))
    A[:,[4]] = np.vstack((u[:,[1]]*(X[:,[0]]*(-u[:,[0]]) - X[:,[1]]*vk[2]*( - u[:,[0]]) + X[:,[2]]*vk[1]*( - u[:,[0]])), X[:,[0]]*vk[1]*( - u[:,[0]]) - X[:,[2]]*( - u[:,[0]]) - u[:,[0]]*(X[:,[0]]*( - u[:,[0]]) - X[:,[1]]*vk[2]*( - u[:,[0]]) + X[:,[2]]*vk[1]*( - u[:,[0]])) - X[:,[1]]*vk[0]*( - u[:,[0]])))
    A[:,[5]] = np.vstack((X[:,[0]]*( - u[:,[0]]) - X[:,[1]]*vk[2]*( - u[:,[0]]) + X[:,[2]]*vk[1]*( - u[:,[0]]), X[:,[1]]*( - u[:,[0]]) + X[:,[0]]*vk[2]*( - u[:,[0]]) - X[:,[2]]*vk[0]*( - u[:,[0]])))
    A[:,[6]] = np.vstack((np.zeros((6,1)), np.ones((6,1))))
    A[:,[7]] = np.vstack((-np.ones((6,1)), np.zeros((6,1))))
    A[:,[8]] = np.vstack((u[:,[1]], -u[:,[0]]))
    A[:,[9]] = np.vstack((np.zeros((6,1)), u[:,[0]]))
    A[:,[10]] = np.vstack((-u[:,[0]], np.zeros((6,1))))
    A[:,[11]] = np.vstack((-u[:,[1]]*( - u[:,[0]]), u[:,[0]]*( - u[:,[0]])))
    A[:,[12]] = np.vstack((X[:,[2]]*u[:,[1]] - X[:,[1]], X[:,[0]] - X[:,[2]]*u[:,[0]]))
    n = null_space(A)
    s = n[12,-1]
    v = n[0:3,[-1]]/s
    w = n[3:6,[-1]]/s
    C = n[6:9,[-1]]/s
    t = n[9:12,[-1]]/s
    return w,t,v,C

def calcR6P2linErrEq(M,data):
    w = M[0:3]
    t = M[3:6]
    v = M[6:9]
    C = M[9:12]
    err = np.zeros(data.shape[1])
    i = 0
    for temp in data.transpose():
        X = temp[0:3,None]
        u = temp[3:5,None]
        u = ut.a2h(u)
        eq = np.matmul(ut.x_(u),np.matmul((np.eye(3)+u[0]*ut.x_(w)),np.matmul((np.eye(3)+ut.x_(v)),X))+C+u[0]*t)
        err[i] = np.sum(np.absolute(eq))
        i+=1
    return err

def r6pLin(data,maxiter=5,eps=1e-10):
    """An implementation of ACCV'18 linear R6P iterative algorithm 

    Computes absolute pose of a RS camera and the motion parameters.

    Parameters
    ----------
    data : 2D numpy array
        a 5x6 array containing 6 3D-2D correspondences such that first 3 rows are the 3D points coordinates and next 2 rows are the 2D coordinates, i.e. 3D points X are data[0:3,0:6] and the correspondning 2D points are data[3:5,0:6].
    maxiter : int, optional
        maximum number of iterations, default is 5
    eps : float, optional
        the algebraic error under which convergence is declared, default is 1e-10
    Returns
    -------
    [M]
        a list with one element (this is for the purposes of RANSAC routines)which is a 12 element numpy array containing the pose and motion parameters. Camera center is C = M[9:12], angle-axis representation of camera orientation is v = M[6:9], rotation velocity is w = M[0:3] and translation velocity is t = M[3:6].

    """
    X = data[0:3,0:6]
    u = data[3:5,0:6]
    v = np.zeros((3,1))
    maxiter = 5
    k = 0
    while k<maxiter:
        w,t,v,C = linwtvC(X.transpose(),u.transpose(),v)
        M = np.vstack((w,t,v,C))
        err = np.sum(calcR6P2linErrEq(M,data))
        if err < eps:
            return [M]
        else:
            k += 1
    return []