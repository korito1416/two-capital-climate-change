"""
file for post HJB with k and y
"""
import os
import sys
sys.path.append("../src/")
import numpy as np
import pandas as pd
# from numba import njit
from src.supportfunctions import finiteDiff
import SolveLinSys



"""
solver.py
For 3D abatement solver
"""
import os
import sys
sys.path.append("../../src/")
from src.Utility import finiteDiff_3D
import SolveLinSys
import numpy as np
import petsc4py
from petsc4py import PETSc
import petsclinearsystem
import time
from datetime import datetime


def pde_one_interation(ksp, petsc_mat, X1_mat_1d, X2_mat_1d, X3_mat_1d, lowerLims, upperLims, dVec, increVec, v0, A, B_1, B_2, B_3, C_1, C_2, C_3, D, tol, epsilon):
    """
    Perform one iteration of solving the PDE using PETSc linear solvers.

    This function constructs the linear system representing the PDE,
    solves it using the provided Krylov subspace solver (`ksp`), and
    updates the value function.

    Parameters
    ----------
    ksp : PETSc.KSP
        Krylov subspace solver object from PETSc.
    petsc_mat : PETSc.Mat
        The PETSc matrix representing the linear system.
    X1_mat_1d, X2_mat_1d, X3_mat_1d : np.ndarray
        1D representations of the grid matrices for each state variable.
    lowerLims, upperLims : np.ndarray
        Lower and upper limits of the state variables.
    dVec : np.ndarray
        Vector of grid step sizes.
    increVec : np.ndarray
        Vector indicating grid increments.
    v0 : np.ndarray
        Current value function.
    A, B_1, B_2, B_3, C_1, C_2, C_3, D : np.ndarray
        Coefficients in the PDE.
    tol : float
        Tolerance for the linear solver convergence.
    epsilon : float
        Step size for updating the value function.

    Returns
    -------
    out_comp : np.ndarray
        Updated value function after solving the linear system.
    end_ksp : float
        Time taken by the PETSc solver.
    bpoint1 : float
        Start time of the iteration.
    """
    bpoint1 = time.time()
    A_1d   = A.ravel(order = 'F')
    C_1_1d = C_1.ravel(order = 'F')
    C_2_1d = C_2.ravel(order = 'F')
    C_3_1d = C_3.ravel(order = 'F')
    B_1_1d = B_1.ravel(order = 'F')
    B_2_1d = B_2.ravel(order = 'F')
    B_3_1d = B_3.ravel(order = 'F')
    D_1d   = D.ravel(order = 'F')
    v0_1d  = v0.ravel(order = 'F')
    petsclinearsystem.formLinearSystem(X1_mat_1d, X2_mat_1d, X3_mat_1d, A_1d, B_1_1d, B_2_1d, B_3_1d, C_1_1d, C_2_1d, C_3_1d, epsilon, lowerLims, upperLims, dVec, increVec, petsc_mat)
    b = v0_1d + D_1d * epsilon
    petsc_rhs = PETSc.Vec().createWithArray(b)
    x = petsc_mat.createVecRight()


    # create linear solver
    start_ksp = time.time()
    ksp.setOperators(petsc_mat)
    ksp.setTolerances(rtol=tol)
    ksp.solve(petsc_rhs, x)
    petsc_rhs.destroy()
    x.destroy()
    out_comp = np.array(ksp.getSolution()).reshape(A.shape,order = "F")
    end_ksp = time.time()
    num_iter = ksp.getIterationNumber()
    return out_comp,end_ksp,bpoint1

def _FOC_update(v0, steps= (), states = (), args=(), controls=(), fraction=0.5):
    """
    Update the First Order Conditions (FOC) based on the current value function.

    This function computes the derivatives of the value function, updates the
    control variables (i), and calculates various coefficients used in the PDE.

    Parameters
    ----------
    v0 : np.ndarray
        Current value function.
    steps : tuple
        Step sizes (hX1, hX2, hX3) for each state variable.
    states : tuple
        Current state matrices (K_mat, Y_mat, L_mat).
    args : tuple
        Model parameters.
    controls : tuple
        Current optimal control (i_star).
    fraction : float
        Fraction used for updating controls to ensure stability.

    Returns
    -------
    A, B_1, B_2, B_3, C_1, C_2, C_3, D : np.ndarray
        Updated coefficients for the PDE.
    dX1, ddX1 : np.ndarray
        First and second order derivatives of the value function with respect to X1.
    ii : np.ndarray
        Updated control variable (i).
    h_k : np.ndarray
        Auxiliary variable related to climate or technology.
    """
    hX1, hX2, hX3 = steps
    K_mat, Y_mat, L_mat = states
    delta, alpha, theta, vartheta_bar, lambda_bar, mu_k, kappa, sigma_k, theta_ell, pi_c_o, pi_c, sigma_y, zeta, psi_0, psi_1, sigma_g, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g, rho, varrho = args

    i_star = controls
    # First order derivative
    dX1  = finiteDiff_3D(v0,0,1,hX1)
    dX1[dX1 <= 1e-16] = 1e-16
    dK = dX1
    dX2  = finiteDiff_3D(v0,1,1,hX2)
    dY = dX2
    dX3  = finiteDiff_3D(v0,2,1,hX3)
    dX3[dX3 <= 1e-16] = 1e-16
    dL = dX3
    ######## second order
    ddX1 = finiteDiff_3D(v0,0,2,hX1)
    ddX2 = finiteDiff_3D(v0,1,2,hX2)
    ddY = ddX2
    ddX3 = finiteDiff_3D(v0,2,2,hX3)


    temp = delta * ( (alpha - i_star)   )**(-rho) 

    i_new = (1- temp/dK)/kappa



    ii = i_new * fraction + i_star * (1 - fraction)


    
    
    h_k = -1/xi_k *sigma_k * dK 
    
    h_k[h_k>=-1e-16]=-1e-16
    

    consumption = alpha - ii 
    consumption[consumption <= 1e-16] = 1e-16
    
    
    A   = -delta * np.ones_like(K_mat)
    
    B_1 = mu_k + ii - 0.5 * kappa * ii**2 - 0.5 * sigma_k**2
    B_1 += sigma_k*h_k
    B_2 = np.zeros_like(K_mat)
    B_3 = np.zeros_like(K_mat)
    

    C_1 = 0.5 * sigma_k**2 * np.ones(K_mat.shape)
    C_2 = np.zeros_like(K_mat)
    C_3 = np.zeros_like(K_mat)
    
    D = delta * np.log(consumption) + delta * K_mat
     
    D += 1/2 * xi_k * h_k**2

    
    
    return A, B_1, B_2, B_3, C_1, C_2, C_3, D, dX1, ddX1, ii, h_k


def hjb_post_tech(
        state_grid=(), model_args=(), V_post_damage=None, 
        tol=1e-8, epsilon=0.1, fraction=0.5, max_iter=10000,
        v0=None,
        smart_guess=None,
        ):
    """
    Solve the Hamilton-Jacobi-Bellman (HJB) equation after a technology shock.

    This function sets up the state grid, initializes variables, and iteratively
    solves the HJB equation using finite difference methods and PETSc solvers
    until convergence is achieved or the maximum number of iterations is reached.

    Parameters
    ----------
    state_grid : tuple
        Tuple containing arrays for each state variable grid (K, Y, L).
    model_args : tuple
        Tuple containing all model parameters.
    V_post_damage : np.ndarray, optional
        Value function after damage shocks. Default is None.
    tol : float, optional
        Tolerance for convergence. Default is 1e-8.
    epsilon : float, optional
        Step size for updating the value function. Default is 0.1.
    fraction : float, optional
        Fraction used for updating controls to ensure stability. Default is 0.5.
    max_iter : int, optional
        Maximum number of iterations. Default is 10000.
    v0 : np.ndarray, optional
        Initial guess for the value function. If None, a default is used.
    smart_guess : dict, optional
        Dictionary containing smart guesses for initialization. Default is None.

    Returns
    -------
    res : dict
        Dictionary containing the solved value function, optimal controls,
        and other auxiliary variables.
    """
    now = datetime.now()
    current_time = now.strftime("%d-%H:%M")
    K, Y, L = state_grid

    delta, alpha, theta, vartheta_bar, lambda_bar, mu_k, kappa, sigma_k, theta_ell, pi_c_o, pi_c, sigma_y, zeta, psi_0, psi_1, sigma_g, gamma_1, gamma_2, gamma_3, y_bar, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g, rho, varrho = model_args


    X1     = K
    nX1    = len(X1)
    hX1    = X1[1] - X1[0]
    X1_min = X1.min()
    X1_max = X1.max()
    X2     = Y
    nX2    = len(X2)
    hX2    = X2[1] - X2[0]
    X2_min = X2.min()
    X2_max = X2.max()
    X3     = L
    nX3    = len(X3)
    hX3    = X3[1] - X3[0]
    X3_min = X3.min()
    X3_max = X3.max()

    print("Grid dimension: [{}, {}, {}]\n".format(nX1, nX2, nX3))
    print("Grid step: [{}, {}, {}]\n".format(hX1, hX2, hX3))
    # Discretization of the state space for numerical PDE solution.
    ######## post jump, 3 states
    (X1_mat, X2_mat, X3_mat) = np.meshgrid(X1, X2, X3, indexing = 'ij')
    stateSpace = np.hstack([X1_mat.reshape(-1,1,order = 'F'), X2_mat.reshape(-1,1,order = 'F'), X3_mat.reshape(-1, 1, order='F')])
    K_mat = X1_mat
    Y_mat = X2_mat
    L_mat = X3_mat
    # For PETSc
    X1_mat_1d = X1_mat.ravel(order='F')
    X2_mat_1d = X2_mat.ravel(order='F')
    X3_mat_1d = X3_mat.ravel(order='F')
    lowerLims = np.array([X1_min, X2_min, X3_min], dtype=np.float64)
    upperLims = np.array([X1_max, X2_max, X3_max], dtype=np.float64)
    #### Model type


    # Initial setup of HJB
    FC_Err   = 1
    epoch    = 0
    ii_scalar = 0.089998 
    
    h_k_scalar = -1/xi_k *sigma_k 
    consumption = alpha - ii_scalar
    constant = np.log(consumption) + (mu_k + ii_scalar - 0.5 * kappa * ii_scalar**2 - 0.5 * sigma_k**2 + sigma_k*h_k_scalar)/delta



    if v0 is None:
        # v0 = K_mat + L_mat - np.average(pi_c_o, axis=0) * Y_mat
        v0 = K_mat + constant
        
        
    i_star = 0.089998*np.ones(K_mat.shape)

    if smart_guess:
        v0     = smart_guess["v0"]
        # i_star = smart_guess["i_star"]

            
    

    dVec = np.array([hX1, hX2, hX3])
    increVec = np.array([1, nX1, nX1 * nX2],dtype=np.int32)

    # FOC_args = (delta, alpha, theta, vartheta_bar, lambda_bar, mu_k, kappa, sigma_k, theta_ell, pi_c_o, pi_c, sigma_y, zeta, psi_0, psi_1, psi_2, sigma_g, V_post_tech, dG, ddG, xi_a, xi_g )

    petsc_mat = PETSc.Mat().create()
    petsc_mat.setType('aij')
    petsc_mat.setSizes([nX1 * nX2 * nX3, nX1 * nX2 * nX3])
    petsc_mat.setPreallocationNNZ(13)
    petsc_mat.setUp()
    ksp = PETSc.KSP()
    ksp.create(PETSc.COMM_WORLD)
    ksp.setType('bcgs')
    ksp.getPC().setType('ilu')
    ksp.setFromOptions()

    # Enter the optimization
    while FC_Err > tol and epoch < max_iter:
        
        FOC_args = (delta, alpha, theta, vartheta_bar, lambda_bar, mu_k, kappa, sigma_k, theta_ell, pi_c_o, pi_c, sigma_y, zeta, psi_0, psi_1, sigma_g, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g, rho, varrho)

        start_ep = time.time()
        A, B_1, B_2, B_3, C_1, C_2, C_3, D, dX1, ddX1, ii, h_k = _FOC_update(v0, steps= (hX1, hX2, hX3), states = (K_mat, Y_mat, L_mat), args=FOC_args, controls=(i_star), fraction=fraction)



        out_comp,end_ksp, bpoint1 = pde_one_interation(
                ksp,
                petsc_mat,X1_mat_1d, X2_mat_1d, X3_mat_1d, 
                lowerLims, upperLims, dVec, increVec,
                v0, A, B_1, B_2, B_3, C_1, C_2, C_3, D, 1e-13, epsilon)
        # if epoch % 1 == 0 and reporterror:
            # Calculating PDE error and False Transient error
        
        PDE_rhs = A * v0 + B_1 * dX1 + C_1 * ddX1  + D
        PDE_Err = np.max(abs(PDE_rhs))
        FC_Err = np.max(abs((out_comp - v0)/ epsilon))
        

        if FC_Err < 1.2*tol:
            
            if epoch%100==0:
                print("-----------------------------------")
                print("---------Epoch {}---------------".format(epoch))
                print("-----------------------------------")
                print("min i: {},\t max i: {}\t".format(ii.min(), ii.max()))
                # print("min e: {},\t max e: {}\t".format(ee.min(), ee.max()))
                # print("min x: {},\t max x: {}\t".format(xx.min(), xx.max()))
                # print("min h: {},\t max h: {}\t".format(h.min(), h.max()))
                print("min hk: {},\t max hk: {}\t".format(h_k.min(), h_k.max()))
                # print("min hj: {},\t max hj: {}\t".format(h_j.min(), h_j.max()))
                print("petsc total: {:.3f}s, Residual Norm is {:g}".format((end_ksp - bpoint1),ksp.getResidualNorm()))
                print("Epoch {:d} (PETSc): PDE Error: {:.10f}; False Transient Error: {:.10f}" .format(epoch, PDE_Err, FC_Err))
                print("Epoch time: {:.4f}".format(time.time() - start_ep))
        elif epoch%10000==0:
            dX1  = finiteDiff_3D(v0,0,1,hX1) 

            print("-----------------------------------")
            print("---------Epoch {}---------------".format(epoch))
            print("-----------------------------------")
            print("min i: {},\t max i: {}\t".format(ii.min(), ii.max()))
            # print("min e: {},\t max e: {}\t".format(ee.min(), ee.max()))
            # print("min x: {},\t max x: {}\t".format(xx.min(), xx.max()))
            # print("min h: {},\t max h: {}\t".format(h.min(), h.max()))
            print("min hk: {},\t max hk: {}\t".format(h_k.min(), h_k.max()))
            print("min dX1: {},\t max dX1: {}\t".format(dX1.min(), dX1.max()))
            # print("min hj: {},\t max hj: {}\t".format(h_j.min(), h_j.max()))
            print("petsc total: {:.3f}s, Residual Norm is {:g}".format((end_ksp - bpoint1),ksp.getResidualNorm()))
            print("Epoch {:d} (PETSc): PDE Error: {:.10f}; False Transient Error: {:.10f}" .format(epoch, PDE_Err, FC_Err))
            print("Epoch time: {:.4f}".format(time.time() - start_ep))
        

        v0     = out_comp
        i_star = ii
        # e_star = ee
        # x_star = xx
        epoch += 1

    dX1  = finiteDiff_3D(v0,0,1,hX1)
    dX1[dX1 <= 1e-16] = 1e-16
    dK = dX1
    ######## second order
    ddX1 = finiteDiff_3D(v0,0,2,hX1)


    
    res = {
            "v0"    : v0,
            "i_star": i_star,
            "h_k": h_k,
            "FC_Err": FC_Err,
            }

    return res



