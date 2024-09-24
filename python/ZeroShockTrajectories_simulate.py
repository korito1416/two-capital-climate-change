import numpy as np
import pandas as pd
import sys
print(sys.path)

sys.path.append('./src')

import pickle
import plotly.graph_objects as go
import plotly.offline as pyo
import matplotlib as mpl
import matplotlib.pyplot as plt
import SolveLinSys
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import CubicSpline
from matplotlib.backends.backend_pdf import PdfPages
from src.Utility import finiteDiff_3D
import os
import argparse
# import src.ResultSolver_CRS
import SolveLinSys


parser = argparse.ArgumentParser(description="xi_r values")
parser.add_argument("--dataname",type=str)
parser.add_argument("--outputname",type=str,default="ReplicateSuri")
parser.add_argument("--pdfname",type=str)

parser.add_argument("--xiaarr",nargs='+', type=float)
parser.add_argument("--xikarr",nargs='+', type=float)
parser.add_argument("--xicarr",nargs='+', type=float)
parser.add_argument("--xijarr",nargs='+', type=float)
parser.add_argument("--xidarr",nargs='+', type=float)
parser.add_argument("--xigarr",nargs='+', type=float)

parser.add_argument("--varrhoarr",nargs='+', type=float)
parser.add_argument("--phi_0", type=float)
parser.add_argument("--rhoarr", nargs='+', type=float)
parser.add_argument("--delta", type=float)

parser.add_argument("--psi0arr",nargs='+',type=float)
parser.add_argument("--psi1arr",nargs='+',type=float)
parser.add_argument("--psi2arr",nargs='+',type=float)

parser.add_argument("--hXarr",nargs='+',type=float)
parser.add_argument("--Xminarr",nargs='+',type=float)
parser.add_argument("--Xmaxarr",nargs='+',type=float)

parser.add_argument("--auto",type=int)
parser.add_argument("--IntPeriod",type=int)

parser.add_argument("--scheme",type=str)
parser.add_argument("--HJB_solution",type=str)


args = parser.parse_args()



# Update = args.Update
IntPeriod = args.IntPeriod
timespan = 1/12

# psi0arr = np.array([0.006,0.009])
# # # psi0arr = np.array([0.009])
# # # psi1arr = np.array([.5,.7,.9])
# psi1arr = np.array([.3,.4])

psi0arr = args.psi0arr
psi1arr = args.psi1arr
psi2arr = args.psi2arr
xiaarr = args.xiaarr
xikarr = args.xikarr 
xicarr = args.xicarr 
xijarr = args.xijarr 
xidarr = args.xidarr 
xigarr = args.xigarr 
varrhoarr = args.varrhoarr
rhoarr = args.rhoarr
 


Xminarr = args.Xminarr
Xmaxarr = args.Xmaxarr
hXarr = args.hXarr

auto = args.auto

scheme = args.scheme
HJB_solution = args.HJB_solution


delta = args.delta
alpha = 0.115
kappa = 6.667
mu_k  = -0.045
# sigma_k = np.sqrt(0.0087**2 + 0.0038**2)
sigma_k = 0.0100
beta_f = 1.86/1000
sigma_y = 1.2 * 1.86 / 1000
zeta = 0.0
# psi_0 = 0.00025
# psi_1 = 1/2
# sigma_g   = 0.016
sigma_g   = 0.0078
gamma_1 = 1.7675 / 1000
gamma_2 = 0.0022 * 2




y_bar = 2.
y_bar_lower = 1.5


# Tech
theta = 3
phi_1 = 3
beta = 0.1206
lambda_bar = 0.1206
# vartheta_bar = 0.0453
# vartheta_bar = 0.05
# vartheta_bar = 0.056
# vartheta_bar = 0.5

phi_0 = args.phi_0
vartheta_bar = args.phi_0


lambda_bar_first = lambda_bar / 2.
vartheta_bar_first = vartheta_bar / 2.

lambda_bar_second = 1e-3
vartheta_bar_second = 0.





K_min = Xminarr[0]
K_max = Xmaxarr[0]
hK    = hXarr[0]
K     = np.arange(K_min, K_max + hK, hK)
nK    = len(K)
Y_min = Xminarr[1]
Y_max = Xmaxarr[1]
hY    = hXarr[1] # make sure it is float instead of int
Y     = np.arange(Y_min, Y_max + hY, hY)
nY    = len(Y)
L_min = Xminarr[2]
L_max = Xmaxarr[2]
hL    = hXarr[2]
L     = np.arange(L_min, L_max+hL,  hL)
nL    = len(L)


id_2 = np.abs(Y - y_bar).argmin()
Y_min_short = Xminarr[3]
Y_max_short = Xmaxarr[3]
Y_short     = np.arange(Y_min_short, Y_max_short + hY, hY)
nY_short    = len(Y_short)

n_bar1 = len(Y_short)-1
n_bar2 = np.abs(Y_short - y_bar).argmin()


# print("bY_short={:d}".format(nY_short))
(K_mat, Y_mat, L_mat) = np.meshgrid(K, Y_short, L, indexing="ij")







mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["figure.figsize"] = (16,10)
mpl.rcParams["font.size"] = 15
mpl.rcParams["legend.frameon"] = False
mpl.style.use('classic')
mpl.rcParams["lines.linewidth"] = 5


print("After, figure default size is: ", plt.rcParams["savefig.bbox"])
print("After, figure default size is: ", plt.rcParams["figure.figsize"])
print("After, figure default dpi is: ", plt.rcParams["figure.dpi"])
print("After, figure default size is: ", plt.rcParams["font.size"])
print("After, legend.frameon is: ", plt.rcParams["legend.frameon"])
print("After, lines.linewidth is: ", plt.rcParams["lines.linewidth"])

def simulate_pre(
    grid = (), 
    model_args = (), 
    controls = (),
    ME = (),
    n_bar = (),  
    initial=(np.log(85/0.115), 1.1, np.log(11.2)),   # Initial state values for K, Y, L
    T0=0, T=40, dt=1/12, # Time parameters: start time, end time, and time step
    printing=True):
    """
    Simulate the pre-jump dynamics of the model over a specified time horizon.

    Parameters:
    - grid: Tuple of arrays representing the discretized state variables (K, Y, L).
    - model_args: Tuple of model parameters used in the simulation.
    - controls: Tuple of control variables and functions obtained from solving the HJB equations.
    - ME: Marginal emission cost function or array.
    - n_bar: Index or limit used to truncate or subset the arrays (e.g., for temperature anomaly Y).
    - initial: Tuple of initial state values (K_0, Y_0, L_0).
    - T0: Initial time of the simulation.
    - T: Total time to simulate.
    - dt: Time step size.
    - printing: Boolean indicating whether to print debug information during the simulation.

    Returns:
    - res: Dictionary containing simulation results and computed variables.
    """
    
    K, Y, L = grid

    if printing==True:
        print("K_min={},K_max={},Y_min={},Y_max={},L_min={},L_max={}" .format(K.min(),K.max(),Y.min(),Y.max(),L.min(),L.max()))

    K_min, K_max, Y_min, Y_max, L_min, L_max = min(K), max(K), min(Y), max(Y), min(L), max(L)
    hK, hY, hL = K[1] - K[0], Y[1] - Y[0], L[1]-L[0]

    '''    # Unpack control variables and functions from the HJB solution
    (
        ii,           # Investment policy function
        ee,           # Emission policy function
        xx,           # R&D investment policy function
        g_tech,       # Technology jump probability adjustment
        g_damage,     # Damage jump probability adjustment
        pi_c,         # Adjusted probabilities for climate models
        h,            # Uncertainty adjustment function for temperature anomaly
        h_k,          # Uncertainty adjustment function for capital
        h_j,          # Uncertainty adjustment function for knowledge capital
        v,            # Value function before technology jump
        v_post_tech_raw,  # Value function after technology jump
    ) = controls'''

    delta, mu_k, kappa, sigma_k, beta_f, zeta, psi_0, psi_1, sigma_g, theta, lambda_bar, vartheta_bar, varrho, xi_a,xi_k,xi_c,xi_j,xi_d,xi_g,rho = model_args
    ii, ee, xx, g_tech, g_damage, pi_c, h, h_k, h_j, v, v_post_tech_raw = controls
    ME_base = ME
    n_bar = n_bar
    K_0, Y_0, L_0 = initial

    
    # #### Temporary Checks
    # print("---------------Temporary Checks Start-----------")
    # ee_modified = ee/(alpha*lambda_bar* np.exp(K_mat))
    # print("ee_modified in [{},{}]".format(ee_modified.min(), ee_modified.max()))
    # # plt.close()
    # # plt.plot(Y, ee_modified[:,:,-1].T)
    # # plt.ylim(0,1)
    # # plt.savefig("./abatement/pdf_2tech/2jump_step_4.00,9.00_0.0,4.0_1.0,6.0_SS_0.2_LR_0.1/AA_ee_modified.png")
    # # plt.close()
    # print("---------------Temporary Checks End-----------")

    Y = Y[:n_bar+1]
    
    ii = ii[:,:n_bar+1,:]
    ee = ee[:,:n_bar+1,:]
    xx = xx[:,:n_bar+1,:]
    g_tech = g_tech[:,:n_bar+1,:]
    g_damage = g_damage[:,:,:n_bar+1,:]
    pi_c = pi_c[:,:,:n_bar+1,:]
    h = h[:,:n_bar+1,:]
    h_k = h_k[:,:n_bar+1,:]
    h_j = h_j[:,:n_bar+1,:]

    v = v[:,:n_bar+1,:]
    


    v_post_tech = v_post_tech_raw
        
    (K_mat, Y_mat, L_mat) = np.meshgrid(K, Y, L, indexing = 'ij')

    # Compute the abatement cost function jj based on emissions ee and other parameters
    jj = alpha * vartheta_bar * (1 - ee / (alpha * lambda_bar * np.exp(K_mat)))**theta
    
    jj[jj <= 1e-16] = 1e-16
    consumption = alpha - ii - jj - xx
    # Compute the total marginal emission cost ME_total
    ME_total = delta/ consumption  * alpha * vartheta_bar * theta * (1 - ee / ( alpha * lambda_bar * np.exp(K_mat)))**(theta - 1) /( alpha * lambda_bar * np.exp(K_mat) )


    years  = np.arange(T0, T0 + T + dt, dt)
    pers   = len(years)
       

    # some parameters remaiend unchanged across runs
    gamma_1  = 0.00017675
    gamma_2  = 2. * 0.0022
    beta_f   = 1.86 / 1000
    sigma_y  = 1.2 * 1.86 / 1000
    
    theta_ell = pd.read_csv("./data/model144.csv", header=None).to_numpy()[:, 0]/1000.
    pi_c_o = np.ones(len(theta_ell)) / len(theta_ell)
    pi_c_o = np.array([temp * np.ones(K_mat.shape) for temp in pi_c_o])
    # theta_ell = np.array([temp * np.ones(K_mat.shape) for temp in theta_ell])

    # Compute finite differences of the value function v with respect to K, Y, and L
    dK = finiteDiff_3D(v, 0,1,hK )
    dY = finiteDiff_3D(v, 1,1,hY )
    dL = finiteDiff_3D(v, 2,1,hL )

    print("dk={},{}".format(dK.min(),dK.max()))
    
    gridpoints = (K, Y, L)

    # Create interpolating functions for control variables and value functions

    i_func = RegularGridInterpolator(gridpoints, ii) # Investment function
    e_func = RegularGridInterpolator(gridpoints, ee) # Emission function
    x_func = RegularGridInterpolator(gridpoints, xx)  # R&D function
    tech_func = RegularGridInterpolator(gridpoints, g_tech)  # Technology function
    
    v_func = RegularGridInterpolator(gridpoints, v)  # Value function
    v_post_tech_func = RegularGridInterpolator(gridpoints, v_post_tech)  # Post-tech value function

    h_func = RegularGridInterpolator(gridpoints, h)  # Uncertainty adjustment function h
    hk_func = RegularGridInterpolator(gridpoints, h_k)  # Uncertainty adjustment function h_k
    hj_func = RegularGridInterpolator(gridpoints, h_j)  # Uncertainty adjustment function h_j
    dK_func = RegularGridInterpolator(gridpoints, dK)  # Derivative of v w.r.t K
    dL_func = RegularGridInterpolator(gridpoints, dL)  # Derivative of v w.r.t L
    dY_func = RegularGridInterpolator(gridpoints, dY)  # Derivative of v w.r.t Y
    ME_total_func = RegularGridInterpolator(gridpoints, ME_total)  # Marginal emission cost function
    ME_base_func = RegularGridInterpolator(gridpoints, ME_base)  # Base marginal emission cost function

    # Create a list of interpolating functions for damage functions
    n_damage = len(g_damage)
    damage_func_list = []
    for i in range(n_damage):
        func_i = RegularGridInterpolator(gridpoints, g_damage[i])
        damage_func_list.append(func_i)
    
    # Create a list of interpolating functions for climate functions
    n_climate = len(pi_c)
    
    climate_func_list = []
    for i in range(n_climate):
        func_i = RegularGridInterpolator(gridpoints, pi_c[i])
        climate_func_list.append(func_i)

    # Define functions to retrieve control variables at a given state
    def get_i(x):
        return i_func(x)

    def get_e(x):
        return e_func(x)
    
    def get_x(x):
        return x_func(x)

    def get_dL(x):
        return dL_func(x)

    # Define drift terms for the stochastic differential equations
    def mu_K(i_x):
        return mu_k + i_x - 0.5 * kappa * i_x ** 2  - 0.5 * sigma_k ** 2
    
    def mu_L(Xt, state):
        # Drift for technology level L
        # return -zeta + psi_0 * Xt **psi_1 * (np.exp( psi_1 * state[0]) )  * np.exp( (psi_2-1) * (state[2] - np.log(448)) ) - 0.5 * sigma_g**2
        return -zeta + psi_0 * Xt **psi_1 * (np.exp( psi_1 * (state[0]-state[2])) )  - 0.5 * sigma_g**2
    
    # Initialize arrays to store simulation results
    hist      = np.zeros([pers, 3])
    i_hist    = np.zeros([pers])
    e_hist    = np.zeros([pers])
    e_annual_hist = np.zeros([T])
    x_hist    = np.zeros([pers])
    scc_hist  = np.zeros([pers])
    gt_tech   = np.zeros([pers])

    vt   = np.zeros([pers])
    v_post_techt   = np.zeros([pers])


    ht   = np.zeros([pers])
    hkt   = np.zeros([pers])
    hjt   = np.zeros([pers])
    dK_hist    = np.zeros([pers])
    dL_hist    = np.zeros([pers])
    dY_hist    = np.zeros([pers])

    gt_dmg    = np.zeros([n_damage, pers])
    pi_c_t = np.zeros([n_climate, pers])
    Ambiguity_mean_undis = np.zeros([pers])
    Ambiguity_mean_dis = np.zeros([pers])
    Ambiguity_mean_dis_h = np.zeros([pers])

    ME_base_hist = np.zeros([pers])
    ME_total_hist = np.zeros([pers])

    mu_K_hist = np.zeros([pers])
    mu_L_hist = np.zeros([pers])
    theta_ell_hist = np.zeros([len(theta_ell),pers])

    # Main simulation loop over time periods
    for tm in range(pers):
        if tm == 0:

            # initial points
            hist[0,:] = [K_0, Y_0, L_0] # logL
            i_hist[0] = get_i(hist[0, :])
            e_hist[0] = get_e(hist[0, :])
            x_hist[0] = get_x(hist[0, :])
            mu_K_hist[0] = mu_K(i_hist[0])
            mu_L_hist[0] = mu_L(x_hist[0], hist[0,:])
            gt_tech[0] = tech_func(hist[0, :])

            vt[0] = v_func(hist[0, :])
            v_post_techt[0] = v_post_tech_func(hist[0, :])


            ht[0] = h_func(hist[0, :])
            hkt[0] = hk_func(hist[0, :])
            hjt[0] = hj_func(hist[0, :])
            dK_hist[tm] = dK_func(hist[0,:])
            dL_hist[tm] = dL_func(hist[0,:])
            dY_hist[tm] = dY_func(hist[0,:])
            # Damage function and climate probability at initial state
            for i in range(n_damage):
                damage_func = damage_func_list[i]
                gt_dmg[i, 0] = damage_func(hist[0, :])
            
            for i in range(n_climate):
                climate_func = climate_func_list[i]
                pi_c_t[i, 0] = climate_func(hist[0, :])
            Ambiguity_mean_undis[tm] = np.mean(theta_ell)
            Ambiguity_mean_dis[tm] = np.average(theta_ell,weights=pi_c_t[:,tm])
            
            ME_total_hist[0] = ME_total_func(hist[0,:])
            ME_base_hist[0] = ME_base_func(hist[0,:])
            theta_ell_hist[:,tm] = theta_ell + sigma_y*ht[tm]
            
            e_annual_hist[0] = get_e(hist[0, :])
            ii=0
        else:
            # Subsequent time steps
            if tm % 12 ==0:

                e_annual_hist[ii] = get_e(hist[tm-1,:])
                ii = ii+1 
            # other periods
            # print(hist[tm-1,:])
            i_hist[tm] = get_i(hist[tm-1,:])
            e_hist[tm] = get_e(hist[tm-1,:])
            x_hist[tm] = get_x(hist[tm-1,:])
            gt_tech[tm] = tech_func(hist[tm-1,:])
            
            vt[tm] = v_func(hist[tm-1, :])
            v_post_techt[tm] = v_post_tech_func(hist[tm-1, :])

            ht[tm] = h_func(hist[tm-1, :])
            hkt[tm] = hk_func(hist[tm-1, :])
            hjt[tm] = hj_func(hist[tm-1, :])
            dK_hist[tm] = dK_func(hist[tm-1,:])
            dL_hist[tm] = dL_func(hist[tm-1,:])
            dY_hist[tm] = dY_func(hist[tm-1,:])

            for i in range(n_damage):
                damage_func = damage_func_list[i]
                gt_dmg[i, tm] = damage_func(hist[tm-1, :])

            for i in range(n_climate):
                climate_func = climate_func_list[i]
                pi_c_t[i, tm] = climate_func(hist[tm -1, :])
                

            # zero shock simulation
            mu_K_hist[tm] = mu_K(i_hist[tm])
            mu_L_hist[tm] = mu_L(x_hist[tm], hist[tm-1, :])

            hist[tm,0] = hist[tm-1,0] + mu_K_hist[tm] * dt #logK
            hist[tm,1] = hist[tm-1,1] + beta_f * e_hist[tm] * dt
            hist[tm,2] = hist[tm-1,2] + mu_L_hist[tm] * dt # logλ
            Ambiguity_mean_undis[tm] = np.mean(theta_ell)
            Ambiguity_mean_dis[tm] = np.average(theta_ell,weights=pi_c_t[:,tm])
            # Ambiguity_mean_dis_h[tm] = np.average(theta_ell + sigma_y*gt_mean[tm],weights=pi_c_t[:,tm])
            ME_total_hist[tm] = ME_total_func(hist[tm,:])
            ME_base_hist[tm] = ME_base_func(hist[tm,:])
            theta_ell_hist[:,tm] = theta_ell + sigma_y*ht[tm]

        if printing==True:
            print("time={}, K={},Y={},L={},mu_K={},mu_Y={},mu_L={},ii={},ee={},xx={},ME_total_base={:.3}" .format(tm, hist[tm,0],hist[tm,1],hist[tm,2],mu_K_hist[tm],beta_f * e_hist[tm],mu_L_hist[tm],i_hist[tm],e_hist[tm],x_hist[tm],np.log(ME_total_hist[tm]/ME_base_hist[tm])*100), flush=True)
            # print("time={}, Vg={},V={},UAD={}" .format(tm, v_post_techt[tm], vt[tm],-xi_g  * (1 - np.exp(-1/xi_g *(v_post_techt[tm]-vt[tm]))) ),  flush=True)
            # print("time={}, E={},RD={}" .format(tm, e_hist[tm], (x_hist[tm] /alpha)*100   ),  flush=True)

    
    
        # using Kt instead of K0
    jt = 1 - e_hist/ (alpha * lambda_bar * np.exp(hist[:, 0]))
    jt[jt <= 1e-16] = 1e-16
    # Left-hand side of the first-order condition
    LHS = theta * vartheta_bar / lambda_bar * jt**(theta -1)
    #  Marginal cost of consumption
    MC = delta / (alpha  - i_hist - alpha * vartheta_bar * jt**theta - x_hist)
    # Compute consumption C
    C = (alpha  - i_hist - alpha * vartheta_bar * jt**theta - x_hist) * np.exp(hist[:, 0])
    # Compute social value of R&D (svrd_hist)
    svrd_hist = np.exp(hist[:,2]) * dL_hist /(delta  * C**(-rho) * np.exp((rho-1)*vt))


    scc_hist = LHS * 1000


    # MU_RD = dL_hist * psi_0* psi_1 * x_hist**(psi_1-1) * np.exp(psi_1*(hist[:,0]-hist[:,2]))

    # scrd_hist = MU_RD/MC*1000
    # Compute social value of R&D (svrd_hist)
    scrd_hist = np.exp(hist[:,2]) * dL_hist / MC * np.exp(hist[:, 0])

    scrd_hist2 = psi_0 * psi_1 * (x_hist * np.exp(hist[:,0] - hist[:, 2]) ) ** psi_1

    spo_hist =    xi_g * np.exp(hist[:,2])/varrho * (1 - gt_tech) 

    v_diff = -xi_g * np.log(gt_tech)
    spo_hist2 =    xi_g * np.exp(hist[:,2])/varrho * (1 - gt_tech + gt_tech * np.log(gt_tech)) +  np.exp(hist[:,2])/varrho * gt_tech * v_diff
    
    uncertainty_adjusted_diff =  -xi_g  * (1 - np.exp(-1/xi_g *(v_post_techt-vt)))

    uncertainty_adjusted_diff2 =  -xi_g  * (1 - gt_tech)

    print(abs(uncertainty_adjusted_diff-uncertainty_adjusted_diff2).max())

    scgw_hist =  - dY_hist/ MC * np.exp(hist[:,0]) * 1000 
    scgw_hist2 = - dY_hist * (Ambiguity_mean_dis + sigma_y*ht)*1000 / MC * np.exp(hist[:,0]) * 1000 


    distorted_tech_intensity = np.exp(hist[:, 2]) * gt_tech/varrho

    distorted_tech_prob = 1 - np.exp(- np.cumsum(np.insert(distorted_tech_intensity * dt, 0, 0) ))[:-1]

    true_tech_intensity = np.exp(hist[:, 2]) /varrho
    true_tech_prob = 1 - np.exp(- np.cumsum(np.insert(true_tech_intensity * dt, 0, 0) ))[:-1]
        
    damage_intensity = Damage_Intensity(hist[:, 1])
    distorted_damage_intensity = np.mean(gt_dmg, axis=0) * damage_intensity
    distorted_damage_prob = 1 - np.exp(- np.cumsum(np.insert(distorted_damage_intensity * dt, 0, 0) ))[:-1]
    
    true_damage_intensity =  damage_intensity
    true_damage_prob = 1 - np.exp(- np.cumsum(np.insert(true_damage_intensity * dt, 0, 0) ))[:-1]
    TA = hist[:,1]

    RD_Plot = ((x_hist * np.exp(hist[:, 0])/(alpha*np.exp(hist[:,0])))*100)
    LogSVRD_Plot = np.log(scrd_hist)
    LogSCGW_Plot = np.log(scgw_hist)

    RelativeEntropy_hk = 1/2*xi_k*hkt**2
    RelativeEntropy_h = 1/2*xi_c*ht**2
    RelativeEntropy_hj = 1/2*xi_j*hjt**2
    
    RelativeEntropy_TechJump = xi_j * np.exp(hist[:, 2])/varrho * (1-gt_tech + gt_tech*np.log(gt_tech))
    RelativeEntropy_DamageJump = xi_d * Damage_Intensity(hist[:, 1]) * np.sum(1 - gt_dmg + gt_dmg*np.log(gt_dmg) , axis=0)/n_damage


    iota = (1-e_hist/(beta*alpha*np.exp(hist[:, 0])))*(e_hist<beta*alpha*np.exp(hist[:, 0]))
    
    c_hist = alpha*np.exp(hist[:, 0])*(1-phi_0* iota**phi_1) - i_hist * np.exp(hist[:, 0]) - x_hist * np.exp(hist[:, 0])
    
    damage_hist = np.exp(gamma_1*hist[:,1]+gamma_2*(hist[:,1])**2/2)

    print(c_hist/(alpha*np.exp(hist[:,0])))
    
    res = dict(
        states= hist, 
        i = i_hist * np.exp(hist[:, 0]), 
        e = e_hist,
        c = c_hist,
        dmg = damage_hist,
        x = x_hist * np.exp(hist[:, 0]),
        scc = scc_hist,
        scrd = scrd_hist,
        svrd_orig = svrd_hist,
        scgw = scgw_hist,
        scgw2 = scgw_hist2,
        scrd_2 = scrd_hist2,
        spo = spo_hist,
        spo2 = spo_hist2,
        uncertainty_adjusted_diff=uncertainty_adjusted_diff,
        uncertainty_adjusted_diff2=uncertainty_adjusted_diff2,
        vt = vt,
        v_post_techt = v_post_techt,
        gt_tech = gt_tech,
        ht = ht,
        hkt = hkt,
        hjt = hjt,
        gt_dmg = gt_dmg,
        distorted_damage_prob=distorted_damage_prob,
        distorted_tech_prob=distorted_tech_prob,
        pic_t = pi_c_t,
        ME_total = ME_total_hist,
        ME_base = ME_base_hist,
        ME_total_base = np.log(ME_total_hist / ME_base_hist ) * 100,
        jt = jt,
        LHS = LHS,
        years=years,
        theta_ell_new = theta_ell_hist,
        true_tech_prob = true_tech_prob,
        true_damage_prob = true_damage_prob,
        Ambiguity_mean_undis = Ambiguity_mean_undis,
        Ambiguity_mean_dis = Ambiguity_mean_dis,
        TA = TA,
        RD_Plot = RD_Plot,
        LogSVRD_Plot=LogSVRD_Plot,
        LogSCGW_Plot = LogSCGW_Plot,
        RelativeEntropy_hk = RelativeEntropy_hk,
        RelativeEntropy_h = RelativeEntropy_h,
        RelativeEntropy_hj = RelativeEntropy_hj,
        RelativeEntropy_TechJump = RelativeEntropy_TechJump,
        RelativeEntropy_DamageJump = RelativeEntropy_DamageJump,
        dL = dL_hist,
        dK = dK_hist
        )
    

    return res

def Damage_Intensity(Yt, y_bar_lower=1.5):
    '''
    Damage Jump Intensity 
    Input: 
    ---------------------------
    Yt: Temperature Anomaly
    Return:
    ---------------------
    Intensity: J^n    
    '''
    r_1 = 1.5
    r_2 = 2.5
    Intensity = r_1 * (np.exp(r_2 / 2 * (Yt - y_bar_lower)**2) -1) * (Yt > y_bar_lower)
    return Intensity






def model_simulation_generate(xi_a,xi_k,xi_c,xi_j,xi_d,xi_g,rho,psi_0,psi_1,varrho):
    """
    Generate model simulation results for given ambiguity and preference parameters.

    Parameters:
    - xi_a, xi_k, xi_c, xi_j, xi_d, xi_g: Ambiguity aversion parameters for different state variables.
    - rho: Relative risk aversion coefficient.
    - psi_0, psi_1: Parameters related to the technology jump intensity.
    - varrho: Parameter for the technology adoption rate.
    
    Returns:
    - res: A dictionary containing the simulation results.
    """
    # Output_Dir = "/scratch/bincheng/"
    Output_Dir = args.outputname
    Data_Dir = Output_Dir+"abatement/data_2tech/"+args.dataname+"/"
    # File_Dir = "xi_a_{}_xi_k_{}_xi_c_{}_xi_j_{}_xi_d_{}_xi_g_{}_psi_0_{}_psi_1_{}_varrho_{}_rho_{}_" .format(xi_a,xi_k,xi_c,xi_j,xi_d,xi_g,psi_0,psi_1,varrho,rho)
    File_Dir = "xi_a_{}_xi_k_{}_xi_c_{}_xi_j_{}_xi_d_{}_xi_g_{}_psi_0_{}_psi_1_{}_varrho_{}_rho_{}_delta_{}_" .format(xi_a, xi_k, xi_c, xi_j, xi_d, xi_g, psi_0,psi_1, varrho, rho, delta)
    


    if scheme == "macroannual":
        xi_a_pre = 100000.
        xi_g_pre = 100000.
        xi_p_pre = 100000.
        File_Name_Suffix_pre = "_xiapre_{}_xig_pre_{}_xippre_{}".format(xi_a_pre, xi_g_pre, xi_p_pre) + "_full_" + scheme + "_" +HJB_solution
        n_bar = n_bar2
        with open(Data_Dir+ File_Dir + "model_tech1_pre_damage"+File_Name_Suffix_pre, "rb") as f:
            model_tech1_pre_damage_ME_base = pickle.load(f)

    elif scheme == "newway":
        xi_a_pre = 100000.
        xi_g_pre = 100000.
        xi_p_pre = 100000.
        File_Name_Suffix_pre = "_xiapre_{}_xig_pre_{}_xippre_{}".format(xi_a_pre, xi_g_pre, xi_p_pre) + "_full_" + scheme + "_" +HJB_solution
        n_bar = n_bar2
        with open(Data_Dir+ File_Dir + "model_tech1_pre_damage"+File_Name_Suffix_pre, "rb") as f:
            model_tech1_pre_damage_ME_base = pickle.load(f)

    elif scheme == "check":
        xi_a_pre = xi_a
        xi_g_pre = xi_g
        xi_p_pre = xi_g
        File_Name_Suffix_pre = "_xiapre_{}_xig_pre_{}_xippre_{}".format(xi_a_pre, xi_g_pre, xi_p_pre) + "_full_" + scheme + "_" +HJB_solution
        n_bar = n_bar1
    
        with open(Data_Dir+ File_Dir + "model_tech1_pre_damage"+File_Name_Suffix_pre, "rb") as f:
            model_tech1_pre_damage_ME_base = pickle.load(f)
    elif scheme== 'direct':
        with open(Data_Dir+ File_Dir + "model_tech1_pre_damage", "rb") as f:
            model_tech1_pre_damage_ME_base = pickle.load(f)
        with open(Data_Dir+ File_Dir + "model_tech2_pre_damage", "rb") as f:
            model_tech2_pre_damage_ME_base = pickle.load(f)
        
        n_bar = n_bar1
        
        
    ME_base = model_tech1_pre_damage_ME_base["ME"]

    # Extract control variables and value functions from the pre-damage model solution
    v = model_tech1_pre_damage_ME_base["v0"]          # Value function before technology jump
    i = model_tech1_pre_damage_ME_base["i_star"]      # Optimal investment in capital
    e = model_tech1_pre_damage_ME_base["e_star"]      # Optimal emission control rate
    x = model_tech1_pre_damage_ME_base["x_star"]      # Optimal R&D investment
    pi_c = model_tech1_pre_damage_ME_base["pi_c"]     # Adjusted probabilities over climate models
    g_tech = model_tech1_pre_damage_ME_base["g_tech"] # Technology jump intensity
    g_damage = model_tech1_pre_damage_ME_base["g_damage"]  # Damage jump probabilities
    h = model_tech1_pre_damage_ME_base["h"]           # Uncertainty adjustment for temperature anomaly
    h_k = model_tech1_pre_damage_ME_base["h_k"]       # Uncertainty adjustment for capital
    h_j = model_tech1_pre_damage_ME_base["h_j"]       # Uncertainty adjustment for knowledge capital

    # Extract the value function after technology adoption
    v_post_tech = model_tech2_pre_damage_ME_base["v0"]



    with open(Data_Dir + File_Dir+"model_tech1_pre_damage", "rb") as f:
        tech1 = pickle.load(f)
    
    # Extract original control variables and value functions up to n_bar+1 in the Y dimension
    v_orig = tech1["v0"][:,:n_bar+1,:]
    i_orig = tech1["i_star"][:,:n_bar+1,:]
    e_orig = tech1["e_star"][:,:n_bar+1,:]
    x_orig = tech1["x_star"][:,:n_bar+1,:]
    pi_c_orig = tech1["pi_c"][:,:,:n_bar+1,:]
    g_tech_orig = tech1["g_tech"][:,:n_bar+1,:]
    g_damage_orig =  tech1["g_damage"][:,:,:n_bar+1,:]



    print("--------------Control Check Start--------------")
    print("Diff_i={}".format(np.max(abs(i-i_orig))))
    print("Diff_e={}".format(np.max(abs(e-e_orig))))
    print("Diff_x={}".format(np.max(abs(x-x_orig))))
    print("--------------Control Check End--------------")
    
    ME_family = ME_base
    
    model_args = (delta, mu_k, kappa,sigma_k, beta_f, zeta, psi_0, psi_1, sigma_g, theta, lambda_bar, vartheta_bar, varrho, xi_a,xi_k,xi_c,xi_j,xi_d,xi_g, rho)



    res = simulate_pre(grid = (K, Y_short, L), 
                       model_args = model_args, 
                       controls = (i,e,x, g_tech, g_damage, pi_c, h, h_k, h_j, v, v_post_tech),
                       ME = ME_family,
                       n_bar = n_bar,  
                       T0=0, 
                       T=IntPeriod, 
                       dt=timespan,printing=True)

    with open(Data_Dir + File_Dir+"model_tech1_pre_damage"+"_UD_simul_{}".format(IntPeriod)+ scheme + "_" +HJB_solution, "wb") as f:
        pickle.dump(res,f)


    
    return res


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):
                for id_rho in range(len(rhoarr)):

                    res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],rhoarr[id_rho],psi0arr[id_psi0],psi1arr[id_psi1],varrhoarr[id_varrho])


