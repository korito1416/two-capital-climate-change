import numpy as np

def finiteDiff_4D(data, dim, order, dlt, cap=None):
    """
    Computes the central difference derivatives for a 4D array along a specified dimension.

    Parameters:
    - data: np.ndarray (4D)
        Input data array.
    - dim: int
        Dimension along which the difference is computed (0, 1, 2, or 3).
    - order: int
        Order of the derivative (1 for first-order, 2 for second-order).
    - dlt: float
        Step size for the finite difference.
    - cap: float, optional
        If provided, the result will be capped at this value.

    Returns:
    - res: np.ndarray (4D)
        The computed central difference along the specified dimension.
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

def finiteDiff_3D(data, dim, order, dlt, cap = None):  
    """
    Computes the central difference derivatives for a 3D array along a specified dimension.

    Parameters:
    - data: np.ndarray (3D)
        Input data array.
    - dim: int
        Dimension along which the difference is computed (0, 1, or 2).
    - order: int
        Order of the derivative (1 for first-order, 2 for second-order).
    - dlt: float
        Step size for the finite difference.
    - cap: float, optional
        If provided, the result will be capped at this value.

    Returns:
    - res: np.ndarray (3D)
        The computed central difference along the specified dimension.
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

def finiteDiff_1D_first(data, dim, dlt, cap = None):  
    """
    Computes the first-order central difference derivative for a 1D array.

    Parameters:
    - data: np.ndarray (1D)
        Input data array.
    - dlt: float
        Step size for the finite difference.
    - cap: float, optional
        If provided, the result will be capped at this value.

    Returns:
    - res: np.ndarray (1D)
        The computed first-order derivative.
    """
    res = np.zeros(data.shape)
    res[1:-1] = (1 / (2 * dlt)) * (data[2:] - data[:-2])
    res[0] = (1 / (dlt)) * (data[1] - data[0])
    res[-1] = (1 / (dlt)) * (data[-1] - data[-2])
    return res

def finiteDiff_1D_second(data, dim, dlt, cap = None):  
    """
    Computes the second-order central difference derivative for a 1D array.

    Parameters:
    - data: np.ndarray (1D)
        Input data array.
    - dlt: float
        Step size for the finite difference.
    - cap: float, optional
        If provided, the result will be capped at this value.

    Returns:
    - res: np.ndarray (1D)
        The computed second-order derivative.
    """
    res = np.zeros(data.shape)
    res[1:-1] = (1 / dlt ** 2) * (data[2:] + data[:-2] - 2 * data[1:-1])
    res[-1] = (1 / dlt ** 2) * (data[-1] + data[-2] - 2 * data[-1])
    res[0] = (1 / dlt ** 2) * (data[1] + data[0] - 2 * data[0])
    return res

def finiteDiff_2D_first(data, dim, dlt, cap = None):  
    """
    Computes the first-order central difference derivative for a 2D array.

    Parameters:
    - data: np.ndarray (2D)
        Input data array.
    - dim: int
        Dimension along which the difference is computed (0 or 1).
    - dlt: float
        Step size for the finite difference.
    - cap: float, optional
        If provided, the result will be capped at this value.

    Returns:
    - res: np.ndarray (2D)
        The computed first-order derivative along the specified dimension.
    """
    res = np.zeros(data.shape)
    if dim == 0:                  # to first dimension

        res[1:-1,:] = (1 / (2 * dlt)) * (data[2:,:] - data[:-2,:])

    elif dim == 1:                # to second dimension

        res[:,1:-1] = (1 / (2 * dlt)) * (data[:,2:] - data[:,:-2])
    return res

def finiteDiff_3D_second(data, dim, dlt, cap = None):  
    """
    Computes the second-order central difference derivative for a 3D array.

    Parameters:
    - data: np.ndarray (3D)
        Input data array.
    - dim: int
        Dimension along which the difference is computed (0, 1, or 2).
    - dlt: float
        Step size for the finite difference.
    - cap: float, optional
        If provided, the result will be capped at this value.

    Returns:
    - res: np.ndarray (3D)
        The computed second-order derivative along the specified dimension.
    """
    res = np.zeros(data.shape)
    if dim == 0:                  # to first dimension
        res[1:-1,:,:] = (1 / dlt ** 2) * (data[2:,:,:] + data[:-2,:,:] - 2 * data[1:-1,:,:])
        res[-1,:,:] = (1 / dlt ** 2) * (data[-1,:,:] + data[-2,:,:] - 2 * data[-1,:,:])
        res[0,:,:] = (1 / dlt ** 2) * (data[1,:,:] + data[0,:,:] - 2 * data[0,:,:])

    elif dim == 1:                # to second dimension

        res[:,1:-1,:] = (1 / dlt ** 2) * (data[:,2:,:] + data[:,:-2,:] - 2 * data[:,1:-1,:])
        res[:,-1,:] = (1 / dlt ** 2) * (data[:,-1,:] + data[:,-2,:] - 2 * data[:,-1,:])
        res[:,0,:] = (1 / dlt ** 2) * (data[:,1,:] + data[:,0,:] - 2 * data[:,0,:])

    elif dim == 2:                # to third dimension

        res[:,:,1:-1] = (1 / dlt ** 2) * (data[:,:,2:] + data[:,:,:-2] - 2 * data[:,:,1:-1])
        res[:,:,-1] = (1 / dlt ** 2) * (data[:,:,-1] + data[:,:,-2] - 2 * data[:,:,-1])
        res[:,:,0] = (1 / dlt ** 2) * (data[:,:,1] + data[:,:,0] - 2 * data[:,:,0])
    
    return res

def finiteDiff_2D_second(data, dim, dlt, cap = None):  
    """
    Computes the second-order central difference derivative for a 2D array.

    Parameters:
    - data: np.ndarray (2D)
        Input data array.
    - dim: int
        Dimension along which the difference is computed (0 or 1).
    - dlt: float
        Step size for the finite difference.
    - cap: float, optional
        If provided, the result will be capped at this value.

    Returns:
    - res: np.ndarray (2D)
        The computed second-order derivative along the specified dimension.
    """
    res = np.zeros(data.shape)
    if dim == 0:                  # to first dimension
        res[1:-1,:] = (1 / dlt ** 2) * (data[2:,:] + data[:-2,:] - 2 * data[1:-1,:])
        res[-1,:] = (1 / dlt ** 2) * (data[-1,:] + data[-2,:] - 2 * data[-1,:])
        res[0,:] = (1 / dlt ** 2) * (data[1,:] + data[0,:] - 2 * data[0,:])

    elif dim == 1:                # to second dimension

        res[:,1:-1] = (1 / dlt ** 2) * (data[:,2:] + data[:,:-2] - 2 * data[:,1:-1])
        res[:,-1] = (1 / dlt ** 2) * (data[:,-1] + data[:,-2] - 2 * data[:,-1])
        res[:,0] = (1 / dlt ** 2) * (data[:,1] + data[:,0] - 2 * data[:,0])
    
    return res

def finiteDiff_2D_cross(data, dlt1, dlt2, cap = None):
    """
    Computes the mixed second-order central difference derivative for a 2D array.

    Parameters:
    - data: np.ndarray (2D)
        Input data array.
    - dlt1: float
        Step size for the first dimension.
    - dlt2: float
        Step size for the second dimension.
    - cap: float, optional
        If provided, the result will be capped at this value.

    Returns:
    - res: np.ndarray (2D)
        The computed mixed derivative between the two dimensions.
    """
    res = np.zeros(data.shape)
    res[1:-1,1:-1] = (1 / (2* dlt1 * dlt2)) * (2*data[1:-1,1:-1] +  data[2:,:-2] + data[:-2,2:] - data[2:, 1:-1] - data[:-2, 1:-1] - data[1:-1, 2:] - data[1:-1, :-2])
    res[-1,1:-1] = (1 / (2* dlt1 * dlt2)) * (2*data[-1,1:-1] + data[-1,:-2] + data[-2,2:] - data[-1, 1:-1] - data[-2, 1:-1] - data[-1, 2:] - data[-1, :-2])
    res[0,1:-1] = (1 / (2* dlt1 * dlt2)) * (2*data[0,1:-1] + data[1,:-2] + data[0,2:] - data[1, 1:-1] - data[0, 1:-1] - data[0, 2:] - data[0, :-2])
    res[1:-1,-1] = (1 / (2* dlt1 * dlt2)) * (2*data[1:-1,-1] + data[2:,-2] + data[:-2,-1] - data[2:, -1] - data[:-2, -1] - data[1:-1, -2] - data[1:-1, -1])
    res[1:-1,0] = (1 / (2* dlt1 * dlt2)) * (2*data[1:-1,0] + data[2:,0] + data[:-2,1] - data[2:, 0] - data[:-2, 0] - data[1:-1, 1] - data[1:-1, 0])

    res[-1,-1] = (1 / (2* dlt1 * dlt2)) * (2*data[-1,-1] + data[-1,-2] + data[-2,-1] - data[-1, -1] - data[-2, -1] - data[-1, -2] - data[-1, -1])
    res[-1,0] = (1 / (2* dlt1 * dlt2)) * (2*data[-1,0] + data[-1,0] + data[-2,1] - data[-1, 0] - data[-2, 0] - data[-1, 1] - data[-1, 0])
    res[0,-1] = (1 / (2* dlt1 * dlt2)) * (2*data[0,-1] + data[1,-2] + data[0,-1] - data[1, -1] - data[0, -1] - data[0, -2] - data[0, -1])
    res[0,0] = (1 / (2* dlt1 * dlt2)) * (2*data[0,0] + data[1,0] + data[0,1] - data[1, 0] - data[0, 0] - data[0, 1] - data[0, 0])

    return res

def finiteDiff_3D_cross(data, dim1, dim2, dlt1, dlt2, cap = None):  
    """
    Computes the mixed second-order central difference derivative for a 3D array.

    Parameters:
    - data: np.ndarray (3D)
        Input data array.
    - dim1: int
        First dimension for the mixed derivative (0, 1, or 2).
    - dim2: int
        Second dimension for the mixed derivative (0, 1, or 2).
    - dlt1: float
        Step size for the first dimension.
    - dlt2: float
        Step size for the second dimension.
    - cap: float, optional
        If provided, the result will be capped at this value.

    Returns:
    - res: np.ndarray (3D)
        The computed mixed derivative between the two dimensions.
    """
    res = np.zeros(data.shape)

    if (dim1 == 0) and (dim2 == 1):

        res[1:-1,1:-1,:] = (1 / (2* dlt1 * dlt2)) * (2*data[1:-1,1:-1,:] +  data[2:,:-2,:] + data[:-2,2:,:] - data[2:, 1:-1, :] - data[:-2, 1:-1, :] - data[1:-1, 2:, :] - data[1:-1, :-2, :])

        res[-1,1:-1,:] = (1 / (2* dlt1 * dlt2)) * (2*data[-1,1:-1,:] + data[-1,:-2,:] + data[-2,2:,:] - data[-1, 1:-1, :] - data[-2, 1:-1, :] - data[-1, 2:, :] - data[-1, :-2, :])
        res[0,1:-1,:] = (1 / (2* dlt1 * dlt2)) * (2*data[0,1:-1,:] + data[1,:-2,:] + data[0,2:,:] - data[1, 1:-1, :] - data[0, 1:-1, :] - data[0, 2:, :] - data[0, :-2, :])
        res[1:-1,-1,:] = (1 / (2* dlt1 * dlt2)) * (2*data[1:-1,-1,:] + data[2:,-2,:] + data[:-2,-1,:] - data[2:, -1, :] - data[:-2, -1, :] - data[1:-1, -2, :] - data[1:-1, -1, :])
        res[1:-1,0,:] = (1 / (2* dlt1 * dlt2)) * (2*data[1:-1,0,:] + data[2:,0,:] + data[:-2,1,:] - data[2:, 0, :] - data[:-2, 0, :] - data[1:-1, 1, :] - data[1:-1, 0, :])

        res[-1,-1,:] = (1 / (2* dlt1 * dlt2)) * (2*data[-1,-1,:] + data[-1,-2,:] + data[-2,-1,:] - data[-1, -1, :] - data[-2, -1, :] - data[-1, -2, :] - data[-1, -1, :])
        res[-1,0,:] = (1 / (2* dlt1 * dlt2)) * (2*data[-1,0,:] + data[-1,0,:] + data[-2,1,:] - data[-1, 0, :] - data[-2, 0, :] - data[-1, 1, :] - data[-1, 0, :])
        res[0,-1,:] = (1 / (2* dlt1 * dlt2)) * (2*data[0,-1,:] + data[1,-2,:] + data[0,-1,:] - data[1, -1, :] - data[0, -1, :] - data[0, -2, :] - data[0, -1, :])
        res[0,0,:] = (1 / (2* dlt1 * dlt2)) * (2*data[0,0,:] + data[1,0,:] + data[0,1,:] - data[1, 0, :] - data[0, 0, :] - data[0, 1, :] - data[0, 0, :])

    elif (dim1 == 0) and (dim2 == 2):

        res[1:-1,:,1:-1] = (1 / (2* dlt1 * dlt2)) * (2*data[1:-1,:,1:-1] + data[2:,:,:-2] + data[:-2,:,2:] - data[2:,:, 1:-1] - data[:-2,:, 1:-1] - data[1:-1,:, 2:] - data[1:-1,:, :-2])

        res[-1,:,1:-1] = (1 / (2* dlt1 * dlt2)) * (2*data[-1,:,1:-1] + data[-1,:,:-2] + data[-2,:,2:] - data[-1,:, 1:-1] - data[-2,:, 1:-1] - data[-1,:, 2:] - data[-1,:, :-2])
        res[0,:,1:-1] = (1 / (2* dlt1 * dlt2)) * (2*data[0,:,1:-1] + data[1,:,:-2] + data[0,:,2:] - data[1,:, 1:-1] - data[0,:, 1:-1] - data[0,:, 2:] - data[0,:, :-2])
        res[1:-1,:,-1] = (1 / (2* dlt1 * dlt2)) * (2*data[1:-1,:,-1] + data[2:,:,-2] + data[:-2,:,-1] - data[2:,:, -1] - data[:-2,:, -1] - data[1:-1,:, -2] - data[1:-1,:, -1])
        res[1:-1,:,0] = (1 / (2* dlt1 * dlt2)) * (2*data[1:-1,:,0] + data[2:,:,0] + data[:-2,:,1] - data[2:,:, 0] - data[:-2,:, 0] - data[1:-1,:, 1] - data[1:-1,:, 0])
        
        res[-1,:,-1] = (1 / (2* dlt1 * dlt2)) * (2*data[-1,:,-1] + data[-1,:,-2] + data[-2,:,-1] - data[-1,:, -1] - data[-2,:, -1] - data[-1,:, -2] - data[-1,:, -1])
        res[-1,:,0] = (1 / (2* dlt1 * dlt2)) * (2*data[-1,:,0] + data[-1,:,0] + data[-2,:,1] - data[-1,:, 0] - data[-2,:, 0] - data[-1,:, 1] - data[-1,:, 0])
        res[0,:,-1] = (1 / (2* dlt1 * dlt2)) * (2*data[0,:,-1] + data[1,:,-2] + data[0,:,-1] - data[1,:, -1] - data[0,:, -1] - data[0,:, -2] - data[0,:, -1])
        res[0,:,0] = (1 / (2* dlt1 * dlt2)) * (2*data[0,:,0] + data[1,:,0] + data[0,:,1] - data[1,:, 0] - data[0,:, 0] - data[0,:, 1] - data[0,:, 0])
    
    elif (dim1 == 1) and (dim2 == 2):
        
        res[:,1:-1,1:-1] = (1 / (2* dlt1 * dlt2)) * (2*data[:,1:-1,1:-1] + data[:,2:,:-2] + data[:,:-2,2:] - data[:,2:, 1:-1] - data[:,:-2, 1:-1] - data[:,1:-1, 2:] - data[:,1:-1, :-2])

        res[:,-1,1:-1] = (1 / (2* dlt1 * dlt2)) * (2*data[:,-1,1:-1] + data[:,-1,:-2] + data[:,-2,2:] - data[:,-1, 1:-1] - data[:,-2, 1:-1] - data[:,-1, 2:] - data[:,-1, :-2])
        res[:,0,1:-1] = (1 / (2* dlt1 * dlt2)) * (2*data[:,0,1:-1] + data[:,1,:-2] + data[:,0,2:] - data[:,1, 1:-1] - data[:,0, 1:-1] - data[:,0, 2:] - data[:,0, :-2])
        res[:,1:-1,-1] = (1 / (2* dlt1 * dlt2)) * (2*data[:,1:-1,-1] + data[:,2:,-2] + data[:,:-2,-1] - data[:,2:, -1] - data[:,:-2, -1] - data[:,1:-1, -2] - data[:,1:-1, -1])
        res[:,1:-1,0] = (1 / (2* dlt1 * dlt2)) * (2*data[:,1:-1,0] + data[:,2:,0] + data[:,:-2,1] - data[:,2:, 0] - data[:,:-2, 0] - data[:,1:-1, 1] - data[:,1:-1, 0])

        res[:,-1,-1] = (1 / (2* dlt1 * dlt2)) * (2*data[:,-1,-1] + data[:,-1,-2] + data[:,-2,-1] - data[:,-1, -1] - data[:,-2, -1] - data[:,-1, -2] - data[:,-1, -1])
        res[:,-1,0] = (1 / (2* dlt1 * dlt2)) * (2*data[:,-1,0] + data[:,-1,0] + data[:,-2,1] - data[:,-1, 0] - data[:,-2, 0] - data[:,-1, 1] - data[:,-1, 0])
        res[:,0,-1] = (1 / (2* dlt1 * dlt2)) * (2*data[:,0,-1] + data[:,1,-2] + data[:,0,-1] - data[:,1, -1] - data[:,0, -1] - data[:,0, -2] - data[:,0, -1])
        res[:,0,0] = (1 / (2* dlt1 * dlt2)) * (2*data[:,0,0] + data[:,1,0] + data[:,0,1] - data[:,1, 0] - data[:,0, 0] - data[:,0, 1] - data[:,0, 0])

    return res

def finiteDiff_3D2(data, dim, order, dlt, cap = None):  
    """
    Computes the central difference derivatives (first or second order) for a 3D array along a specified dimension.

    Parameters:
    - data: np.ndarray (3D)
        Input data array.
    - dim: int
        Dimension along which the difference is computed (0, 1, or 2).
    - order: int
        Order of the derivative (1 for first-order, 2 for second-order).
    - dlt: float
        Step size for the finite difference.
    - cap: float, optional
        If provided, the result will be capped at this value.

    Returns:
    - res: np.ndarray (3D)
        The computed central difference along the specified dimension.
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
                res[-1,:,:] = (1 / dlt ** 2) * (data[-2,:,:] + 0 - data[-1,:,:])
                res[0,:,:] = (1 / dlt ** 2) * (data[1,:,:] + 0 - data[0,:,:])

            elif dim == 1:                # to second dimension

                res[:,1:-1,:] = (1 / dlt ** 2) * (data[:,2:,:] + data[:,:-2,:] - 2 * data[:,1:-1,:])
                res[:,-1,:] = (1 / dlt ** 2) * (data[:,-2,:] + 0 -data[:,-1,:])
                res[:,0,:] = (1 / dlt ** 2) * (data[:,1,:] + 0 - data[:,0,:])

            elif dim == 2:                # to third dimension

                res[:,:,1:-1] = (1 / dlt ** 2) * (data[:,:,2:] + data[:,:,:-2] - 2 * data[:,:,1:-1])
                res[:,:,-1] = (1 / dlt ** 2) * (data[:,:,-2] + 0 - data[:,:,-1])
                res[:,:,0] = (1 / dlt ** 2) * (data[:,:,1] + 0 - data[:,:,0])

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


