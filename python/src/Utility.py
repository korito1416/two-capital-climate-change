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


def finite_1D(data, dlt, cap = None):  
    '''
    This function finite_1D computes the finite difference approximation of the first derivative 
        for a one-dimensional array data. The central difference method is used for the interior points, 
    and forward/backward difference methods are used for the boundary points.

    Input:
    -----------------------------------
    data: a one-dimensional array of values.
    dlt: the grid spacing (the difference between consecutive points in the array).
    cap (optional): An optional threshold, where results lower than cap are set to cap.

    Return: 
    -----------------------------------
    res: containing the approximations of the first derivative for each point in the input array data.
    '''
    res = np.zeros(data.shape)

    res[1:-1] = (1 / (2 * dlt)) * (data[2:] - data[:-2])
    res[-1] = (1 / dlt) * (data[-1] - data[-2])
    res[0] = (1 / dlt) * (data[1] - data[0])
    if cap is not None:
        res[res < cap] = cap
    return res


def finiteDiff(data, dim, order, dlt, cap = None):  
    """
    Compute finite difference approximations of derivatives for multi-dimensional arrays.

    Parameters:
    - data: n-dimensional NumPy array containing function values.
    - dim: the dimension along which to compute the derivative (0-based indexing).
    - order: the order of the derivative (1 for first-order, 2 for second-order).
    - dlt: the grid spacing (Δx).
    - cap (optional): a lower bound to cap the results (if provided).

    The function supports data arrays of dimensions 2, 3, and 4.

    Returns:
    - res: array of the same shape as data, containing the finite difference approximations.
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
    Compute finite difference approximations of derivatives for multi-dimensional arrays (2D and 3D).

    Parameters:
    - data (np.ndarray): 2D or 3D NumPy array containing function values at discrete grid points.
    - dim (int): The dimension along which to compute the derivative (0-based indexing).
                 For 3D arrays:
                     0 - first dimension (e.g., x-axis)
                     1 - second dimension (e.g., y-axis)
                     2 - third dimension (e.g., z-axis)
                 For 2D arrays:
                     0 - first dimension (e.g., x-axis)
                     1 - second dimension (e.g., y-axis)
    - order (int): The order of the derivative to compute.
                   1 for first-order derivatives,
                   2 for second-order derivatives.
    - dlt (float): The grid spacing (Δx), representing the distance between adjacent grid points.
    - cap (float, optional): A lower bound to cap the results. If provided, any values in the result
                             less than `cap` are set to `cap`.

    Returns:
    - res (np.ndarray): An array of the same shape as `data`, containing the finite difference
                        approximations of the derivative.

    Raises:
    - ValueError: If an unsupported number of dimensions, dimension index, or order is provided.
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

def finiteDiff_4D(data, dim, order, dlt, cap=None):
    """
    Compute finite difference approximations of derivatives for multi-dimensional arrays (2D, 3D, and 4D).

    Parameters:
    - data (np.ndarray): 2D, 3D, or 4D NumPy array containing function values at discrete grid points.
    - dim (int): The dimension along which to compute the derivative (0-based indexing).
                For 4D arrays:
                    0 - first dimension (e.g., x-axis)
                    1 - second dimension (e.g., y-axis)
                    2 - third dimension (e.g., z-axis)
                    3 - fourth dimension (e.g., w-axis)
                For 3D arrays:
                    0 - first dimension (e.g., x-axis)
                    1 - second dimension (e.g., y-axis)
                    2 - third dimension (e.g., z-axis)
                For 2D arrays:
                    0 - first dimension (e.g., x-axis)
                    1 - second dimension (e.g., y-axis)
    - order (int): The order of the derivative to compute.
                  1 for first-order derivatives,
                  2 for second-order derivatives.
    - dlt (float): The grid spacing (Δx), representing the distance between adjacent grid points.
    - cap (float, optional): A lower bound to cap the results. If provided, any values in the result
                             less than `cap` are set to `cap`.

    Returns:
    - res (np.ndarray): An array of the same shape as `data`, containing the finite difference
                        approximations of the derivative.

    Raises:
    - ValueError: If an unsupported number of dimensions, dimension index, or order is provided.
    """
    res = np.zeros(data.shape)
    l = len(data.shape)
    if l == 4:
        if order == 1:  # first order derivatives
            if dim == 0:  # to first dimension
                res[1:-1, :, :, :] = (1 / (2 * dlt)) * (data[2:, :, :, :] - data[:-2, :, :, :])
                res[-1, :, :, :] = (1 / dlt) * (data[-1, :, :, :] - data[-2, :, :, :])
                res[0, :, :, :] = (1 / dlt) * (data[1, :, :, :] - data[0, :, :, :])

            elif dim == 1:  # to second dimension
                res[:, 1:-1, :, :] = (1 / (2 * dlt)) * (data[:, 2:, :, :] - data[:, :-2, :, :])
                res[:, -1, :, :] = (1 / dlt) * (data[:, -1, :, :] - data[:, -2, :, :])
                res[:, 0, :, :] = (1 / dlt) * (data[:, 1, :, :] - data[:, 0, :, :])

            elif dim == 2:  # to third dimension
                res[:, :, 1:-1, :] = (1 / (2 * dlt)) * (data[:, :, 2:, :] - data[:, :, :-2, :])
                res[:, :, -1, :] = (1 / dlt) * (data[:, :, -1, :] - data[:, :, -2, :])
                res[:, :, 0, :] = (1 / dlt) * (data[:, :, 1, :] - data[:, :, 0, :])

            elif dim == 3:  # to fourth dimension
                res[:, :, :, 1:-1] = (1 / (2 * dlt)) * (data[:, :, :, 2:] - data[:, :, :, :-2])
                res[:, :, :, -1] = (1 / dlt) * (data[:, :, :, -1] - data[:, :, :, -2])
                res[:, :, :, 0] = (1 / dlt) * (data[:, :, :, 1] - data[:, :, :, 0])

            else:
                raise ValueError('wrong dim')

        elif order == 2:
            if dim == 0:
                res[1:-1, :, :, :] = (1 / dlt ** 2) * (data[2:, :, :, :] + data[:-2, :, :, :] - 2 * data[1:-1, :, :, :])
                res[-1, :, :, :] = (1 / dlt ** 2) * (data[-1, :, :, :] + data[-3, :, :, :] - 2 * data[-2, :, :, :])
                res[0, :, :, :] = (1 / dlt ** 2) * (data[2, :, :, :] + data[0, :, :, :] - 2 * data[1, :, :, :])

            elif dim == 1:
                res[:, 1:-1, :, :] = (1 / dlt ** 2) * (data[:, 2:, :, :] + data[:, :-2, :, :] - 2 * data[:, 1:-1, :, :])
                res[:, -1, :, :] = (1 / dlt ** 2) * (data[:, -1, :, :] + data[:, -3, :, :] - 2 * data[:, -2, :, :])
                res[:, 0, :, :] = (1 / dlt ** 2) * (data[:, 2, :, :] + data[:, 0, :, :] - 2 * data[:, 1, :, :])

            elif dim == 2:
                res[:, :, 1:-1, :] = (1 / dlt ** 2) * (data[:, :, 2:, :] + data[:, :, :-2, :] - 2 * data[:, :, 1:-1, :])
                res[:, :, -1, :] = (1 / dlt ** 2) * (data[:, :, -1, :] + data[:, :, -3, :] - 2 * data[:, :, -2, :])
                res[:, :, 0, :] = (1 / dlt ** 2) * (data[:, :, 2, :] + data[:, :, 0, :] - 2 * data[:, :, 1, :])

            elif dim == 3:
                res[:, :, :, 1:-1] = (1 / dlt ** 2) * (data[:, :, :, 2:] + data[:, :, :, :-2] - 2 * data[:, :, :, 1:-1])
                res[:, :, :, -1] = (1 / dlt ** 2) * (data[:, :, :, -1] + data[:, :, :, -3] - 2 * data[:, :, :, -2])
                res[:, :, :, 0] = (1 / dlt ** 2) * (data[:, :, :, 2] + data[:, :, :, 0] - 2 * data[:, :, :, 1])

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
    Compute finite difference derivatives for multi-dimensional data.

    This function calculates the first or second order finite difference 
    derivatives of the input `data` array along a specified dimension `dim`. 
    It supports 2D, 3D, and 4D numpy arrays and handles boundary conditions 
    using forward and backward differences.

    Parameters
    ----------
    data : numpy.ndarray
        The input data array. Supported dimensions are 2, 3, or 4.
    dim : int
        The dimension along which to compute the derivative. 
        For a 3D array, valid values are 0, 1, or 2.
        For a 4D array, valid values are 0, 1, 2, or 3.
    order : int
        The order of the derivative to compute.
        - `1` for first-order derivatives.
        - `2` for second-order derivatives.
    dlt : float
        The grid spacing (Δx) in the specified dimension.
    cap : float, optional
        If provided, any derivative values below this threshold are set to `cap`.

    Returns
    -------
    numpy.ndarray
        An array of the same shape as `data`, containing the computed derivatives.

    Raises
    ------
    ValueError
        If an unsupported dimension (`dim`) or order (`order`) is specified, 
        or if the input array has a dimension other than 2, 3, or 4.

    Notes
    -----
    - **First-Order Derivative (order=1):**
        - **Interior Points:** Uses forward difference.
        - **Boundary Points:** Uses forward difference at the first point and 
          backward difference at the last point.
    - **Second-Order Derivative (order=2):**
        - **Interior Points:** Uses central difference.
        - **Boundary Points:** Uses one-sided differences to approximate the second derivative.
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
    Solve Partial Differential Equations (PDEs) using specified numerical methods.

    Parameters
    ----------
    stateSpace : numpy.ndarray
        Represents the discretized state space for the PDE. The structure and dimensionality depend on the specific problem being solved.
    
    A : numpy.ndarray
        Coefficient matrix associated with the PDE. Its shape and content are determined by the PDE formulation.
    
    B_r, B_f, B_k : numpy.ndarray
        Additional coefficient matrices corresponding to different terms in the PDE. These matrices are horizontally stacked to form the combined coefficient matrix `B`.
    
    C_rr, C_ff, C_kk : numpy.ndarray
        Coefficient matrices related to nonlinear or coupling terms in the PDE. These are horizontally stacked to form the combined coefficient matrix `C`.
    
    D : numpy.ndarray
        Source or forcing term in the PDE.
    
    v0 : numpy.ndarray
        Initial condition for the PDE solver. It represents the starting state from which the solution evolves.
    
    ε : float, optional (default=1)
        A parameter that may represent a scaling factor, regularization term, or other problem-specific coefficient.
    
    tol : float, optional (default=-10)
        Tolerance level for the solver convergence. The interpretation of this value depends on the solver implementation.
    
    smartguess : bool, optional (default=False)
        Determines the number of iterations for the Feynman-Kac solver. If `True`, a minimal number of iterations (`iters=1`) is used, suitable for scenarios where an initial guess is already close to the solution. If `False`, a larger number of iterations (`iters=400000`) is performed to ensure convergence.
    
    solverType : str, optional (default='False Transient')
        Specifies the numerical method to use for solving the PDE. Supported values:
        - `'False Transient'`: Utilizes the False Transient method for steady-state solutions.
        - `'Feyman Kac'`: Employs the Feynman-Kac probabilistic approach to solve the PDE.
    
    Returns
    -------
    out : numpy.ndarray
        The solution to the PDE, structured in the same shape as the input `stateSpace`. The exact nature of the output depends on the solver used.
    
    Raises
    ------
    ValueError
        If an unsupported `solverType` is specified.
    
    Notes
    -----
    - **False Transient Method**:
        - Transforms the steady-state PDE into a time-dependent problem by introducing a pseudo-time derivative.
        - Iteratively advances the solution in pseudo-time until steady-state is achieved.
    
    - **Feynman-Kac Method**:
        - Solves the PDE using stochastic processes and probabilistic representations.
        - Suitable for high-dimensional PDEs where traditional grid-based methods are computationally expensive.
    
    - **Data Reshaping**:
        - All input matrices (`A`, `B_r`, `B_f`, `B_k`, `C_rr`, `C_ff`, `C_kk`, `D`, `v0`) are reshaped into column vectors with Fortran (column-major) order to ensure compatibility with the solver functions in `SolveLinSys`.
    
    - **Clipping and Regularization**:
        - The `ε` and `tol` parameters allow for control over the solver's regularization and convergence criteria, respectively.
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
    Solve Partial Differential Equations (PDEs) in a 4D state space using specified numerical methods.

    Parameters
    ----------
    stateSpace : numpy.ndarray
        A multidimensional array representing the discretized state space for the PDE. For a 4D PDE, this typically includes dimensions such as time, space, and other relevant variables.
    
    A : numpy.ndarray
        The primary coefficient matrix associated with the PDE. Its shape and content are determined by the specific PDE formulation.
    
    B_1, B_2, B_3, B_4 : numpy.ndarray
        Additional coefficient matrices corresponding to different terms or dimensions in the PDE. These matrices are horizontally stacked to form the combined coefficient matrix `B`.
    
    C_11, C_22, C_33, C_44 : numpy.ndarray
        Coefficient matrices related to nonlinear or coupling terms in the PDE. These are horizontally stacked to form the combined coefficient matrix `C`.
    
    D : numpy.ndarray
        The source or forcing term in the PDE, represented as a multidimensional array matching the shape of `stateSpace`.
    
    v0 : numpy.ndarray
        The initial condition for the PDE solver. It represents the starting state from which the solution evolves and should have the same shape as `stateSpace`.
    
    ε : float, optional (default=1)
        A scaling factor or regularization parameter that influences the solver's behavior. The exact interpretation depends on the numerical method used.
    
    tol : float, optional (default=-10)
        The tolerance level for solver convergence. This parameter determines when the iterative solver should terminate based on the residual's magnitude.
    
    smartguess : bool, optional (default=False)
        Determines the number of iterations for the Feynman-Kac solver. If set to `True`, the solver performs a minimal number of iterations (`iters=1`), assuming that the initial guess is already close to the solution. If `False`, a larger number of iterations (`iters=400000`) is used to ensure convergence.
    
    solverType : str, optional (default='False Transient')
        Specifies the numerical method to use for solving the PDE. Supported values:
        - `'False Transient'`: Utilizes the False Transient Method for steady-state solutions.
        - `'Feyman Kac'`: Employs the Feynman-Kac Method, which leverages stochastic processes to solve the PDE.

    Returns
    -------
    numpy.ndarray
        An array containing the solution to the PDE, structured in the same shape as the input `stateSpace`. The exact nature of the output depends on the solver used.

    Raises
    ------
    ValueError
        If an unsupported `solverType` is specified. The function only supports `'False Transient'` and `'Feyman Kac'` as valid solver types.

    Notes
    -----
    - **False Transient Method**:
        - Transforms a steady-state PDE into a time-dependent problem by introducing a pseudo-time derivative.
        - Iteratively advances the solution in pseudo-time until a steady-state is achieved based on the specified tolerance.
    
    - **Feynman-Kac Method**:
        - Solves the PDE using stochastic processes and probabilistic representations.
        - Particularly useful for high-dimensional PDEs where traditional grid-based methods are computationally expensive.
    
    - **Data Reshaping**:
        - All input matrices (`A`, `B_1`, `B_2`, `B_3`, `B_4`, `C_11`, `C_22`, `C_33`, `C_44`, `D`, `v0`) are reshaped into column vectors with Fortran (column-major) order to ensure compatibility with the solver functions in `SolveLinSys`.
        - The reshaping process flattens the multidimensional arrays while preserving their data order, which is crucial for accurate numerical computations.
    
    - **Clipping and Regularization**:
        - The `ε` and `tol` parameters allow for control over the solver's regularization and convergence criteria, respectively.
        - Proper tuning of these parameters can enhance solver stability and performance.
    
    - **Solver Configuration**:
        - The `smartguess` parameter is particularly relevant for the Feynman-Kac method, where a good initial guess can significantly reduce the number of required iterations.
 
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
