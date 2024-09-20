import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.stats import norm
import pdb
import warnings
import time
import pickle
import SolveLinSys
# from numba import njit

def finiteDiff(data, dim, order, dlt, cap = None):  
    """
    Computes finite difference approximations of derivatives for multidimensional data arrays.
    Supports data arrays of dimensions 2, 3, or 4 and can compute first or second-order derivatives
    along a specified dimension.

    Parameters:
    - data: np.ndarray
        The input data array (2D, 3D, or 4D).
    - dim: int
        The dimension along which to compute the derivative (0-based index).
    - order: int
        The order of the derivative (1 for first-order, 2 for second-order).
    - dlt: float
        The grid spacing (delta) in the specified dimension.
    - cap: float, optional
        If provided, any values in the result less than 'cap' are set to 'cap'.

    Returns:
    - res: np.ndarray
        An array of the same shape as 'data' containing the computed derivatives.
    """
    res = np.zeros(data.shape)
    l = len(data.shape)
    if l == 3:
        if order == 1:                    # first order derivatives
            
            if dim == 0:                  # to first dimension

                res[1:-1,:,:] = (1 / (2 * dlt)) * (data[2:,:,:] - data[:-2,:,:])
                res[-1,:,:] = (1 / dlt) * (data[-1,:,:] - data[-2,:,:])
                res[0,:,:] = (1 / dlt) * (data[1,:,:] - data[0,:,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1,:] = (1 / (2 * dlt)) * (data[:,2:,:] - data[:,:-2,:])
                res[:,-1,:] = (1 / dlt) * (data[:,-1,:] - data[:,-2,:])
                res[:,0,:] = (1 / dlt) * (data[:,1,:] - data[:,0,:])

            elif dim == 2:                # to third dimension

                res[:,:,1:-1] = (1 / (2 * dlt)) * (data[:,:,2:] - data[:,:,:-2])
                res[:,:,-1] = (1 / dlt) * (data[:,:,-1] - data[:,:,-2])
                res[:,:,0] = (1 / dlt) * (data[:,:,1] - data[:,:,0])

            else:
                raise ValueError('wrong dim')
                
        elif order == 2:
            
            if dim == 0:                  # to first dimension

                res[1:-1,:,:] = (1 / dlt ** 2) * (data[2:,:,:] + data[:-2,:,:] - 2 * data[1:-1,:,:])
                res[-1,:,:] = (1 / dlt ** 2) * (data[-1,:,:] + data[-3,:,:] - 2 * data[-2,:,:])
                res[0,:,:] = (1 / dlt ** 2) * (data[2,:,:] + data[0,:,:] - 2 * data[1,:,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1,:] = (1 / dlt ** 2) * (data[:,2:,:] + data[:,:-2,:] - 2 * data[:,1:-1,:])
                res[:,-1,:] = (1 / dlt ** 2) * (data[:,-1,:] + data[:,-3,:] - 2 * data[:,-2,:])
                res[:,0,:] = (1 / dlt ** 2) * (data[:,2,:] + data[:,0,:] - 2 * data[:,1,:])

            elif dim == 2:                # to third dimension

                res[:,:,1:-1] = (1 / dlt ** 2) * (data[:,:,2:] + data[:,:,:-2] - 2 * data[:,:,1:-1])
                res[:,:,-1] = (1 / dlt ** 2) * (data[:,:,-1] + data[:,:,-3] - 2 * data[:,:,-2])
                res[:,:,0] = (1 / dlt ** 2) * (data[:,:,2] + data[:,:,0] - 2 * data[:,:,1])

            else:
                raise ValueError('wrong dim')
            
        else:
            raise ValueError('wrong order')
    elif l == 2:
        if order == 1:                    # first order derivatives
            
            if dim == 0:                  # to first dimension

                res[1:-1,:] = (1 / (2 * dlt)) * (data[2:,:] - data[:-2,:])
                res[-1,:] = (1 / dlt) * (data[-1,:] - data[-2,:])
                res[0,:] = (1 / dlt) * (data[1,:] - data[0,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1] = (1 / (2 * dlt)) * (data[:,2:] - data[:,:-2])
                res[:,-1] = (1 / dlt) * (data[:,-1] - data[:,-2])
                res[:,0] = (1 / dlt) * (data[:,1] - data[:,0])

            else:
                raise ValueError('wrong dim')
                
        elif order == 2:
            
            if dim == 0:                  # to first dimension

                res[1:-1,:] = (1 / dlt ** 2) * (data[2:,:] + data[:-2,:] - 2 * data[1:-1,:])
                res[-1,:] = (1 / dlt ** 2) * (data[-1,:] + data[-3,:] - 2 * data[-2,:])
                res[0,:] = (1 / dlt ** 2) * (data[2,:] + data[0,:] - 2 * data[1,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1] = (1 / dlt ** 2) * (data[:,2:] + data[:,:-2] - 2 * data[:,1:-1])
                res[:,-1] = (1 / dlt ** 2) * (data[:,-1] + data[:,-3] - 2 * data[:,-2])
                res[:,0] = (1 / dlt ** 2) * (data[:,2] + data[:,0] - 2 * data[:,1])

            else:
                raise ValueError('wrong dim')
            
        else:
            raise ValueError('wrong order')

    elif l == 4:
        if order == 1:                    # first order derivatives
            
            if dim == 0:                  # to first dimension

                res[1:-1,:,:,:] = (1 / (2 * dlt)) * (data[2:,:,:,:] - data[:-2,:,:,:])
                res[-1,:,:,:] = (1 / dlt) * (data[-1,:,:,:] - data[-2,:,:,:])
                res[0,:,:,:] = (1 / dlt) * (data[1,:,:,:] - data[0,:,:,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1,:,:] = (1 / (2 * dlt)) * (data[:,2:,:,:] - data[:,:-2,:,:])
                res[:,-1,:,:] = (1 / dlt) * (data[:,-1,:,:] - data[:,-2,:,:])
                res[:,0,:,:] = (1 / dlt) * (data[:,1,:,:] - data[:,0,:,:])

            elif dim == 2:                # to third dimension

                res[:,:,1:-1,:] = (1 / (2 * dlt)) * (data[:,:,2:,:] - data[:,:,:-2,:])
                res[:,:,-1,:] = (1 / dlt) * (data[:,:,-1,:] - data[:,:,-2,:])
                res[:,:,0,:] = (1 / dlt) * (data[:,:,1,:] - data[:,:,0,:])
            
            elif dim == 3:

                res[:,:,:, 1:-1] = (1 / (2 * dlt)) * (data[:,:,:,2:] - data[:,:,:,:-2])
                res[:,:,:,-1] = (1 / dlt) * (data[:,:,:,-1]- data[:, :, :, -2])
                res[:,:,:,0] = (1 / dlt) * (data[:,:,:,1] - data[:,:,:,0])
                

            else:
                raise ValueError('wrong dim')
                
        elif order == 2:
            
            if dim == 0:                  # to first dimension

                res[1:-1,:,:,:] = (1 / dlt ** 2) * (data[2:,:,:,:] + data[:-2,:,:,:] - 2 * data[1:-1,:,:,:])
                res[-1,:,:,:] = (1 / dlt ** 2) * (data[-1,:,:,:] + data[-3,:,:, :] - 2 * data[-2,:,:,:])
                res[0,:,:,:] = (1 / dlt ** 2) * (data[2,:,:,:] + data[0,:,:,:] - 2 * data[1,:,:,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1,:,:] = (1 / dlt ** 2) * (data[:,2:,:,:] + data[:,:-2,:,:] - 2 * data[:,1:-1,:,:])
                res[:,-1,:,:] = (1 / dlt ** 2) * (data[:,-1,:,:] + data[:,-3,:,:] - 2 * data[:,-2,:,:])
                res[:,0,:,:] = (1 / dlt ** 2) * (data[:,2,:,:] + data[:,0,:,:] - 2 * data[:,1,:,:])

            elif dim == 2:                # to third dimension

                res[:,:,1:-1,:] = (1 / dlt ** 2) * (data[:,:,2:,:] + data[:,:,:-2,:] - 2 * data[:,:,1:-1,:])
                res[:,:,-1,:] = (1 / dlt ** 2) * (data[:,:,-1,:] + data[:,:,-3,:] - 2 * data[:,:,-2,:])
                res[:,:,0,:] = (1 / dlt ** 2) * (data[:,:,2,:] + data[:,:,0,:] - 2 * data[:,:,1,:])

            elif dim == 3:                # to third dimension

                res[:,:,:,1:-1] = (1 / dlt ** 2) * (data[:,:,:,2:] + data[:,:,:,:-2] - 2 * data[:,:,:,1:-1])
                res[:,:,:,-1] = (1 / dlt ** 2) * (data[:,:,:,-1] + data[:,:,:,-3] - 2 * data[:,:,:,-2])
                res[:,:,:,0] = (1 / dlt ** 2) * (data[:,:,:,2] + data[:,:,:,0] - 2 * data[:,:,:,1])
            else:
                raise ValueError('wrong dim')
            
        else:
            raise ValueError('wrong order')
    else:
        raise ValueError("Dimension NOT supported")
        
    if cap is not None:
        res[res < cap] = cap
    return res

def finiteDiff_3D(data, dim, order, dlt, cap = None):  
    """
    Specialized function for computing finite difference approximations for 3D data arrays.
    Computes first or second-order derivatives along a specified dimension.

    Parameters:
    - data: np.ndarray (3D)
        The input data array.
    - dim: int
        The dimension along which to compute the derivative (0, 1, or 2).
    - order: int
        The order of the derivative (1 or 2).
    - dlt: float
        The grid spacing (delta) in the specified dimension.
    - cap: float, optional
        If provided, any values in the result less than 'cap' are set to 'cap'.

    Returns:
    - res: np.ndarray (3D)
        An array containing the computed derivatives.
    """
    res = np.zeros(data.shape)
    l = len(data.shape)
    if l == 3:
        if order == 1:                    # first order derivatives
            
            if dim == 0:                  # to first dimension

                res[1:-1,:,:] = (1 / (2 * dlt)) * (data[2:,:,:] - data[:-2,:,:])
                res[-1,:,:] = (1 / dlt) * (data[-1,:,:] - data[-2,:,:])
                res[0,:,:] = (1 / dlt) * (data[1,:,:] - data[0,:,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1,:] = (1 / (2 * dlt)) * (data[:,2:,:] - data[:,:-2,:])
                res[:,-1,:] = (1 / dlt) * (data[:,-1,:] - data[:,-2,:])
                res[:,0,:] = (1 / dlt) * (data[:,1,:] - data[:,0,:])

            elif dim == 2:                # to third dimension

                res[:,:,1:-1] = (1 / (2 * dlt)) * (data[:,:,2:] - data[:,:,:-2])
                res[:,:,-1] = (1 / dlt) * (data[:,:,-1] - data[:,:,-2])
                res[:,:,0] = (1 / dlt) * (data[:,:,1] - data[:,:,0])

            else:
                raise ValueError('wrong dim')
                
        elif order == 2:
            
            if dim == 0:                  # to first dimension

                res[1:-1,:,:] = (1 / dlt ** 2) * (data[2:,:,:] + data[:-2,:,:] - 2 * data[1:-1,:,:])
                res[-1,:,:] = (1 / dlt ** 2) * (data[-1,:,:] + data[-3,:,:] - 2 * data[-2,:,:])
                res[0,:,:] = (1 / dlt ** 2) * (data[2,:,:] + data[0,:,:] - 2 * data[1,:,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1,:] = (1 / dlt ** 2) * (data[:,2:,:] + data[:,:-2,:] - 2 * data[:,1:-1,:])
                res[:,-1,:] = (1 / dlt ** 2) * (data[:,-1,:] + data[:,-3,:] - 2 * data[:,-2,:])
                res[:,0,:] = (1 / dlt ** 2) * (data[:,2,:] + data[:,0,:] - 2 * data[:,1,:])

            elif dim == 2:                # to third dimension

                res[:,:,1:-1] = (1 / dlt ** 2) * (data[:,:,2:] + data[:,:,:-2] - 2 * data[:,:,1:-1])
                res[:,:,-1] = (1 / dlt ** 2) * (data[:,:,-1] + data[:,:,-3] - 2 * data[:,:,-2])
                res[:,:,0] = (1 / dlt ** 2) * (data[:,:,2] + data[:,:,0] - 2 * data[:,:,1])

            else:
                raise ValueError('wrong dim')
            
        else:
            raise ValueError('wrong order')
    elif l == 2:
        if order == 1:                    # first order derivatives
            
            if dim == 0:                  # to first dimension

                res[1:-1,:] = (1 / (2 * dlt)) * (data[2:,:] - data[:-2,:])
                res[-1,:] = (1 / dlt) * (data[-1,:] - data[-2,:])
                res[0,:] = (1 / dlt) * (data[1,:] - data[0,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1] = (1 / (2 * dlt)) * (data[:,2:] - data[:,:-2])
                res[:,-1] = (1 / dlt) * (data[:,-1] - data[:,-2])
                res[:,0] = (1 / dlt) * (data[:,1] - data[:,0])

            else:
                raise ValueError('wrong dim')
                
        elif order == 2:
            
            if dim == 0:                  # to first dimension

                res[1:-1,:] = (1 / dlt ** 2) * (data[2:,:] + data[:-2,:] - 2 * data[1:-1,:])
                res[-1,:] = (1 / dlt ** 2) * (data[-1,:] + data[-3,:] - 2 * data[-2,:])
                res[0,:] = (1 / dlt ** 2) * (data[2,:] + data[0,:] - 2 * data[1,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1] = (1 / dlt ** 2) * (data[:,2:] + data[:,:-2] - 2 * data[:,1:-1])
                res[:,-1] = (1 / dlt ** 2) * (data[:,-1] + data[:,-3] - 2 * data[:,-2])
                res[:,0] = (1 / dlt ** 2) * (data[:,2] + data[:,0] - 2 * data[:,1])

            else:
                raise ValueError('wrong dim')
            
        else:
            raise ValueError('wrong order')

            
    else:
        raise ValueError("Dimension NOT supported")
        
    if cap is not None:
        res[res < cap] = cap
    return res


def DiffOne(data, dim, order, dlt, cap = None):  
    """
    Computes finite differences using a forward difference scheme for first-order derivatives.
    Supports data arrays of dimensions 2, 3, or 4.

    Parameters:
    - data: np.ndarray
        The input data array (2D, 3D, or 4D).
    - dim: int
        The dimension along which to compute the derivative (0-based index).
    - order: int
        The order of the derivative (1 or 2).
    - dlt: float
        The grid spacing in the specified dimension.
    - cap: float, optional
        If provided, any values in the result less than 'cap' are set to 'cap'.

    Returns:
    - res: np.ndarray
        An array containing the computed derivatives.
    """
    res = np.zeros(data.shape)
    l = len(data.shape)
    if l == 3:
        if order == 1:                    # first order derivatives
            
            if dim == 0:                  # to first dimension

                res[1:-1,:,:] = (1 / dlt) * (data[1:,:,:] - data[:-1,:,:])
                res[-1,:,:] = (1 / dlt) * (data[-1,:,:] - data[-2,:,:])
                res[0,:,:] = (1 / dlt) * (data[1,:,:] - data[0,:,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1,:] = (1 / dlt) * (data[:,1:,:] - data[:,:-1,:])
                res[:,-1,:] = (1 / dlt) * (data[:,-1,:] - data[:,-2,:])
                res[:,0,:] = (1 / dlt) * (data[:,1,:] - data[:,0,:])

            elif dim == 2:                # to third dimension

                res[:,:,1:-1] = (1 / dlt) * (data[:,:,1:] - data[:,:,:-1])
                res[:,:,-1] = (1 / dlt) * (data[:,:,-1] - data[:,:,-2])
                res[:,:,0] = (1 / dlt) * (data[:,:,1] - data[:,:,0])

            else:
                raise ValueError('wrong dim')
                
        elif order == 2:
            
            if dim == 0:                  # to first dimension

                res[1:-1,:,:] = (1 / dlt ** 2) * (data[2:,:,:] + data[:-2,:,:] - 2 * data[1:-1,:,:])
                res[-1,:,:] = (1 / dlt ** 2) * (data[-1,:,:] + data[-3,:,:] - 2 * data[-2,:,:])
                res[0,:,:] = (1 / dlt ** 2) * (data[2,:,:] + data[0,:,:] - 2 * data[1,:,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1,:] = (1 / dlt ** 2) * (data[:,2:,:] + data[:,:-2,:] - 2 * data[:,1:-1,:])
                res[:,-1,:] = (1 / dlt ** 2) * (data[:,-1,:] + data[:,-3,:] - 2 * data[:,-2,:])
                res[:,0,:] = (1 / dlt ** 2) * (data[:,2,:] + data[:,0,:] - 2 * data[:,1,:])

            elif dim == 2:                # to third dimension

                res[:,:,1:-1] = (1 / dlt ** 2) * (data[:,:,2:] + data[:,:,:-2] - 2 * data[:,:,1:-1])
                res[:,:,-1] = (1 / dlt ** 2) * (data[:,:,-1] + data[:,:,-3] - 2 * data[:,:,-2])
                res[:,:,0] = (1 / dlt ** 2) * (data[:,:,2] + data[:,:,0] - 2 * data[:,:,1])

            else:
                raise ValueError('wrong dim')
            
        else:
            raise ValueError('wrong order')
    elif l == 2:
        if order == 1:                    # first order derivatives
            
            if dim == 0:                  # to first dimension

                res[1:-1,:] = (1 / dlt) * (data[1:,:] - data[:-1,:])
                res[-1,:] = (1 / dlt) * (data[-1,:] - data[-2,:])
                res[0,:] = (1 / dlt) * (data[1,:] - data[0,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1] = (1 / dlt) * (data[:,1:] - data[:,:-1])
                res[:,-1] = (1 / dlt) * (data[:,-1] - data[:,-2])
                res[:,0] = (1 / dlt) * (data[:,1] - data[:,0])

            else:
                raise ValueError('wrong dim')
                
        elif order == 2:
            
            if dim == 0:                  # to first dimension

                res[1:-1,:] = (1 / dlt ** 2) * (data[2:,:] + data[:-2,:] - 2 * data[1:-1,:])
                res[-1,:] = (1 / dlt ** 2) * (data[-1,:] + data[-3,:] - 2 * data[-2,:])
                res[0,:] = (1 / dlt ** 2) * (data[2,:] + data[0,:] - 2 * data[1,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1] = (1 / dlt ** 2) * (data[:,2:] + data[:,:-2] - 2 * data[:,1:-1])
                res[:,-1] = (1 / dlt ** 2) * (data[:,-1] + data[:,-3] - 2 * data[:,-2])
                res[:,0] = (1 / dlt ** 2) * (data[:,2] + data[:,0] - 2 * data[:,1])

            else:
                raise ValueError('wrong dim')
            
        else:
            raise ValueError('wrong order')

    elif l == 4:
        if order == 1:                    # first order derivatives
            
            if dim == 0:                  # to first dimension

                res[1:-1,:,:,:] = (1 / dlt) * (data[1:-1,:,:,:] - data[:-2,:,:,:])
                res[-1,:,:,:] = (1 / dlt) * (data[-1,:,:,:] - data[-2,:,:,:])
                res[0,:,:,:] = (1 / dlt) * (data[1,:,:,:] - data[0,:,:,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1,:,:] = (1 / dlt) * (data[:,1:-1,:,:] - data[:,:-2,:,:])
                res[:,-1,:,:] = (1 / dlt) * (data[:,-1,:,:] - data[:,-2,:,:])
                res[:,0,:,:] = (1 / dlt) * (data[:,1,:,:] - data[:,0,:,:])

            elif dim == 2:                # to third dimension

                res[:,:,1:-1,:] = (1 /dlt) * (data[:,:,1:-1,:] - data[:,:,:-2,:])
                res[:,:,-1,:] = (1 / dlt) * (data[:,:,-1,:] - data[:,:,-2,:])
                res[:,:,0,:] = (1 / dlt) * (data[:,:,1,:] - data[:,:,0,:])
            
            elif dim == 3:

                res[:,:,:, 1:-1] = (1 / dlt) * (data[:,:,:,1:-1] - data[:,:,:,:-2])
                res[:,:,:,-1] = (1 / dlt) * (data[:,:,:,-1]- data[:, :, :, -2])
                res[:,:,:,0] = (1 / dlt) * (data[:,:,:,1] - data[:,:,:,0])
                

            else:
                raise ValueError('wrong dim')
                
        elif order == 2:
            
            if dim == 0:                  # to first dimension

                res[1:-1,:,:,:] = (1 / dlt ** 2) * (data[2:,:,:,:] + data[:-2,:,:,:] - 2 * data[1:-1,:,:,:])
                res[-1,:,:,:] = (1 / dlt ** 2) * (data[-1,:,:,:] + data[-3,:,:, :] - 2 * data[-2,:,:,:])
                res[0,:,:,:] = (1 / dlt ** 2) * (data[2,:,:,:] + data[0,:,:,:] - 2 * data[1,:,:,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1,:,:] = (1 / dlt ** 2) * (data[:,2:,:,:] + data[:,:-2,:,:] - 2 * data[:,1:-1,:,:])
                res[:,-1,:,:] = (1 / dlt ** 2) * (data[:,-1,:,:] + data[:,-3,:,:] - 2 * data[:,-2,:,:])
                res[:,0,:,:] = (1 / dlt ** 2) * (data[:,2,:,:] + data[:,0,:,:] - 2 * data[:,1,:,:])

            elif dim == 2:                # to third dimension

                res[:,:,1:-1,:] = (1 / dlt ** 2) * (data[:,:,2:,:] + data[:,:,:-2,:] - 2 * data[:,:,1:-1,:])
                res[:,:,-1,:] = (1 / dlt ** 2) * (data[:,:,-1,:] + data[:,:,-3,:] - 2 * data[:,:,-2,:])
                res[:,:,0,:] = (1 / dlt ** 2) * (data[:,:,2,:] + data[:,:,0,:] - 2 * data[:,:,1,:])

            elif dim == 3:                # to third dimension

                res[:,:,:,1:-1] = (1 / dlt ** 2) * (data[:,:,:,2:] + data[:,:,:,:-2] - 2 * data[:,:,:,1:-1])
                res[:,:,:,-1] = (1 / dlt ** 2) * (data[:,:,:,-1] + data[:,:,:,-3] - 2 * data[:,:,:,-2])
                res[:,:,:,0] = (1 / dlt ** 2) * (data[:,:,:,2] + data[:,:,:,0] - 2 * data[:,:,:,1])
            else:
                raise ValueError('wrong dim')
            
        else:
            raise ValueError('wrong order')
    else:
        raise ValueError("Dimension NOT supported")
        
    if cap is not None:
        res[res < cap] = cap
    return res




def PDESolver(stateSpace, A, B_r, B_f, B_k, C_rr, C_ff, C_kk, D, v0, ε = 1, tol = -10, smartguess = False, solverType = 'False Transient'):
   """
    Solves partial differential equations (PDEs) using either the False Transient method or the
    Feynman-Kac method for three-dimensional problems.

    Parameters:
    - stateSpace: tuple
        The grid representing the state space over which the PDE is defined.
    - A: np.ndarray
        Zeroth-order term (constant term) in the PDE.
    - B_r, B_f, B_k: np.ndarray
        Coefficients of the first-order derivative terms in the PDE along different dimensions.
    - C_rr, C_ff, C_kk: np.ndarray
        Coefficients of the second-order derivative terms (diffusion terms) in the PDE.
    - D: np.ndarray
        Source term or forcing term in the PDE.
    - v0: np.ndarray
        Initial guess for the solution of the PDE.
    - ε: float, default 1
        Controls the convergence speed in the False Transient method.
    - tol: float, default -10
        Tolerance level for convergence (e.g., -10 corresponds to 1e-10).
    - smartguess: bool, default False
        If True, uses a smart initial guess with fewer iterations.
    - solverType: str, default 'False Transient'
        Solver method to use ('False Transient' or 'Feyman Kac').

    Returns:
    - out: np.ndarray
        The computed solution to the PDE.
    """
    if solverType == 'False Transient':
        A = A.reshape(-1,1,order = 'F')
        B = np.hstack([B_r.reshape(-1,1,order = 'F'),B_f.reshape(-1,1,order = 'F'),B_k.reshape(-1,1,order = 'F')])
        C = np.hstack([C_rr.reshape(-1,1,order = 'F'), C_ff.reshape(-1,1,order = 'F'), C_kk.reshape(-1,1,order = 'F')])
        D = D.reshape(-1,1,order = 'F')
        v0 = v0.reshape(-1,1,order = 'F')
        out = SolveLinSys.solveFT(stateSpace, A, B, C, D, v0, ε, tol)

        return out

    elif solverType == 'Feyman Kac':
        
        if smartguess:
            iters = 1
        else:
            iters = 400000
            
        A = A.reshape(-1, 1, order='F')
        B = np.hstack([B_r.reshape(-1, 1, order='F'), B_f.reshape(-1, 1, order='F'), B_k.reshape(-1, 1, order='F')])
        C = np.hstack([C_rr.reshape(-1, 1, order='F'), C_ff.reshape(-1, 1, order='F'), C_kk.reshape(-1, 1, order='F')])
        D = D.reshape(-1, 1, order='F')
        v0 = v0.reshape(-1, 1, order='F')
        out = SolveLinSys.solveFK(stateSpace, A, B, C, D, v0, iters)

        return out

def PDESolver_4D(stateSpace, A, B_1, B_2, B_3, B_4, C_11, C_22, C_33, C_44,  D, v0, ε = 1, tol = -10, smartguess = False, solverType = 'False Transient'):
    """
    Extends PDESolver to handle four-dimensional PDEs, suitable for problems involving four state variables.

    Parameters:
    - stateSpace: tuple
        The grid representing the four-dimensional state space over which the PDE is defined.
    - A: np.ndarray
        Zeroth-order term in the PDE.
    - B_1, B_2, B_3, B_4: np.ndarray
        Coefficients of the first-order derivative terms along the four dimensions.
    - C_11, C_22, C_33, C_44: np.ndarray
        Coefficients of the second-order derivative terms along the four dimensions.
    - D: np.ndarray
        Source term in the PDE.
    - v0: np.ndarray
        Initial guess for the solution of the PDE.
    - ε: float, default 1
        Controls the convergence speed in the False Transient method.
    - tol: float, default -10
        Tolerance level for convergence.
    - smartguess: bool, default False
        If True, uses a smart initial guess.
    - solverType: str, default 'False Transient'
        Solver method to use ('False Transient' or 'Feyman Kac').

    Returns:
    - out: np.ndarray
        The computed solution to the four-dimensional PDE.
    """
    if solverType == 'False Transient':
        A = A.reshape(-1,1,order = 'F')
        B = np.hstack([B_1.reshape(-1,1,order = 'F'),B_2.reshape(-1,1,order = 'F'), B_3.reshape(-1,1,order = 'F'), B_4.reshape(-1,1,order='F')])
        C = np.hstack([C_11.reshape(-1,1,order = 'F'), C_22.reshape(-1,1,order = 'F'), C_33.reshape(-1,1,order = 'F'), C_44.reshape(-1,1,order='F')])
        D = D.reshape(-1,1,order = 'F')
        v0 = v0.reshape(-1,1,order = 'F')
        out = SolveLinSys.solveFT(stateSpace, A, B, C, D, v0, ε, tol)

        return out

    elif solverType == 'Feyman Kac':
        
        if smartguess:
            iters = 1
        else:
            iters = 400000
            
        A = A.reshape(-1, 1, order='F')
        B = np.hstack([B_1.reshape(-1, 1, order='F'), B_2.reshape(-1, 1, order='F'), B_3.reshape(-1, 1, order='F'), B_4.reshape(-1,1,order='F')])
        C = np.hstack([C_11.reshape(-1, 1, order='F'), C_22.reshape(-1, 1, order='F'), C_33.reshape(-1, 1, order='F'), C_44.reshape(-1,1,order='F')])
        D = D.reshape(-1, 1, order='F')
        v0 = v0.reshape(-1, 1, order='F')
        out = SolveLinSys.solveFK(stateSpace, A, B, C, D, v0, iters)

        return out
