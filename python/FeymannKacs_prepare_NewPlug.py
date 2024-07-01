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
from src.Utility import finiteDiff_4D
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
parser.add_argument("--xia2arr",nargs='+', type=float)
parser.add_argument("--xik2arr",nargs='+', type=float)
parser.add_argument("--xic2arr",nargs='+', type=float)
parser.add_argument("--xij2arr",nargs='+', type=float)
parser.add_argument("--xid2arr",nargs='+', type=float)
parser.add_argument("--xig2arr",nargs='+', type=float)

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
xia2arr = args.xia2arr
xik2arr = args.xik2arr 
xic2arr = args.xic2arr 
xij2arr = args.xij2arr 
xid2arr = args.xid2arr 
xig2arr = args.xig2arr 
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


y_underline = 1.5
y_overline = 2 


gamma_3_max = 0


y_limit = y_overline

LHS_ylimitlower = gamma_1 * Y_short + gamma_2/2 * Y_short**2 # y<y_limit
LHS_ylimitupper = gamma_1 * Y_short + gamma_2*y_overline * \
    (Y_short-y_limit) + (gamma_2+gamma_3_max)/2 * \
    (Y_short-y_limit)**2 + gamma_2/2 * y_limit**2

LHS_ylimit =LHS_ylimitlower*(Y_short<y_limit) + LHS_ylimitupper*(Y_short>y_limit)

# logN = np.linspace(LHS_ylimit.min(), LHS_ylimit.max(), len(K))
logN = np.linspace(LHS_ylimit.min(), LHS_ylimit.max(), 50)
nlogN = len(logN)
hlogN = logN[1]-logN[0]

print(logN)
# print("bY_short={:d}".format(nY_short))
# (K_mat, Y_mat, L_mat) = np.meshgrid(K, Y_short, L, indexing="ij")
(K_mat, Y_mat, L_mat, N_mat) = np.meshgrid(K, Y_short, L, logN, indexing="ij")


r_1 = 1.5
r_2 = 2.5
Intensity = r_1 * (np.exp(r_2 / 2 * (Y_mat - y_bar_lower)**2) -1) * (Y_mat > y_bar_lower)


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


def Damage_Intensity(Yt, y_bar_lower=1.5):
    r_1 = 1.5
    r_2 = 2.5
    Intensity = r_1 * (np.exp(r_2 / 2 * (Yt - y_bar_lower)**2) -1) * (Yt > y_bar_lower)
    return Intensity




def model_simulation_generate(xi_a,xi_k,xi_c,xi_j,xi_d,xi_g,xi_a2,xi_k2,xi_c2,xi_j2,xi_d2,xi_g2,rho,psi_0,psi_1,varrho):

    # Output_Dir = "/scratch/bincheng/"
    Output_Dir = args.outputname
    Data_Dir = Output_Dir+"abatement/data_2tech/"+args.dataname+"/"
    # File_Dir = "xi_a_{}_xi_k_{}_xi_c_{}_xi_j_{}_xi_d_{}_xi_g_{}_psi_0_{}_psi_1_{}_varrho_{}_rho_{}_" .format(xi_a,xi_k,xi_c,xi_j,xi_d,xi_g,psi_0,psi_1,varrho,rho)
    # File_Dir = "xi_a_{}_xi_k_{}_xi_c_{}_xi_j_{}_xi_d_{}_xi_g_{}_psi_0_{}_psi_1_{}_varrho_{}_rho_{}_delta_{}_" .format(xi_a, xi_k, xi_c, xi_j, xi_d, xi_g, psi_0,psi_1, varrho, rho, delta)
    File_Dir2 = "xi_k_{}_xi_c_{}_xi_j_{}_xi_d_{}_xi_g_{}_xi_k2_{}_xi_c2_{}_xi_j2_{}_xi_d2_{}_xi_g2_{}_psi_0_{}_psi_1_{}_varrho_{}_rho_{}_delta_{}_" .format(xi_k, xi_c, xi_j, xi_d, xi_g, xi_k2, xi_c2, xi_j2, xi_d2, xi_g2, psi_0,psi_1, varrho, rho, delta)


    with open(Data_Dir+ File_Dir2 + "model_tech1_pre_damage", "rb") as f:
        model_tech1_pre_damage = pickle.load(f)
    with open(Data_Dir+ File_Dir2 + "model_tech1_post_damage", "rb") as f:
        model_tech1_post_damage = pickle.load(f)
    with open(Data_Dir+ File_Dir2 + "model_tech2_pre_damage", "rb") as f:
        model_tech2_pre_damage = pickle.load(f)
        
    print("start adding damage dimension")

    # h_y = model_tech1_pre_damage["h"]
    # h_k = model_tech1_pre_damage["h_k"]
    # h_Ig = model_tech1_pre_damage["h_j"]

    # f_tech = model_tech1_pre_damage["g_tech"]
    # f_damage = model_tech1_pre_damage["g_damage"]

    # v_pre_damage = model_tech1_pre_damage["v0"]

    h_y = np.repeat(model_tech1_pre_damage["h"][..., np.newaxis], nlogN, axis=-1)
    h_k = np.repeat(model_tech1_pre_damage["h_k"][..., np.newaxis], nlogN, axis=-1)
    h_Ig = np.repeat(model_tech1_pre_damage["h_j"][..., np.newaxis], nlogN, axis=-1)

    f_tech = np.repeat(model_tech1_pre_damage["g_tech"][..., np.newaxis], nlogN, axis=-1)
    f_damage = np.repeat(model_tech1_pre_damage["g_damage"][..., np.newaxis], nlogN, axis=-1)


    gamma_3_list = np.linspace(0, 1/3, len(f_damage))
    pi_d_o = np.ones(len(gamma_3_list)) / len(gamma_3_list)
    pi_d_o = np.array([temp * np.ones(K_mat.shape) for temp in pi_d_o ])


    theta_ell = pd.read_csv('./data/model144.csv', header=None).to_numpy()[:, 0]/1000.
    pi_c_o = np.ones(len(theta_ell)) / len(theta_ell)
    pi_c_o = np.array([temp * np.ones(K_mat.shape) for temp in pi_c_o])
    pi_c = pi_c_o.copy()
    theta_ell = np.array([temp * np.ones(K_mat.shape) for temp in theta_ell])

    v_pre_damage = np.repeat(model_tech1_pre_damage["v0"][..., np.newaxis], nlogN, axis=-1) - N_mat

    v_post_tech = np.repeat(model_tech2_pre_damage["v0"][..., np.newaxis], nlogN, axis=-1) - N_mat
    
    i = np.repeat(model_tech1_pre_damage["i_star"][..., np.newaxis], nlogN, axis=-1)
    e = np.repeat(model_tech1_pre_damage["e_star"][..., np.newaxis], nlogN, axis=-1)
    x = np.repeat(model_tech1_pre_damage["x_star"][..., np.newaxis], nlogN, axis=-1)
        

    v_i = []
    v_gamma_3_i_diff = []
    
    for model in model_tech1_post_damage:
        # v_post_damage_i = model["v0"]
        v_post_damage_i = np.repeat(model["v0"][..., np.newaxis], nlogN, axis=-1)
        v_post_damage_temp = np.zeros((nK, nY_short, nL, nlogN))
        # v_gamma_3_i_diff_temp = np.zeros((nK, nY_short, nL))
        for j in range(nY_short):
            v_post_damage_temp[:, j, :, :] = v_post_damage_i[:, id_2, :, :] - N_mat[:, j, :, :]
        # v_gamma_3_i_diff_temp = v_post_damage_temp - v_pre_damage
        v_i.append(v_post_damage_temp)
        # v_gamma_3_i_diff.append(v_gamma_3_i_diff_temp)
    v_i = np.array(v_i)
    # v_gamma_3_i_diff = np.array(v_gamma_3_i_diff)


    print("start drift and diffusion")


    mu_logK_baseline = mu_k + i - 0.5 * kappa * i**2 - 0.5 * sigma_k**2
    mu_logK_distorted = mu_k + i - 0.5 * kappa * i**2 - 0.5 * sigma_k**2 + sigma_k * h_k
    sigma_logK = sigma_k*np.ones_like(K_mat)

    dmu_logK_dx1 = finiteDiff_4D(mu_logK_baseline,0,1,hK) + finiteDiff_4D(sigma_logK,0,1,hK) * h_k
    dmu_logK_dx2 = finiteDiff_4D(mu_logK_baseline,1,1,hY) + finiteDiff_4D(sigma_logK,1,1,hY) * h_k
    dmu_logK_dx3 = finiteDiff_4D(mu_logK_baseline,2,1,hL) + finiteDiff_4D(sigma_logK,2,1,hL) * h_k
    dmu_logK_dx4 = finiteDiff_4D(mu_logK_baseline,3,1,hlogN) + finiteDiff_4D(sigma_logK,3,1,hlogN) * h_k

    dsigma_logK_dx1 = finiteDiff_4D(sigma_logK,0,1,hK)
    dsigma_logK_dx2 = finiteDiff_4D(sigma_logK,1,1,hY)
    dsigma_logK_dx3 = finiteDiff_4D(sigma_logK,2,1,hL)
    dsigma_logK_dx4 = finiteDiff_4D(sigma_logK,3,1,hlogN)

    # mu_logK_baseline = mu_k + model_tech1_pre_damage["i_star"] - 0.5 * kappa * model_tech1_pre_damage["i_star"]**2 - 0.5 * sigma_k**2
    # mu_logK_distorted = mu_k + model_tech1_pre_damage["i_star"] - 0.5 * kappa * model_tech1_pre_damage["i_star"]**2 - 0.5 * sigma_k**2 + sigma_k * model_tech1_pre_damage["h_k"]


    mu_Y_baseline = np.sum(theta_ell * pi_c, axis=0) * e
    mu_Y_distorted = np.sum(theta_ell * pi_c, axis=0) * e + sigma_y * h_y * e
    sigma_Y = sigma_y* e


    dmu_Y_dx1 = finiteDiff_4D(mu_Y_baseline,0,1,hK) + finiteDiff_4D(sigma_Y,0,1,hK) * h_y
    dmu_Y_dx2 = finiteDiff_4D(mu_Y_baseline,1,1,hY) + finiteDiff_4D(sigma_Y,1,1,hY) * h_y
    dmu_Y_dx3 = finiteDiff_4D(mu_Y_baseline,2,1,hL) + finiteDiff_4D(sigma_Y,2,1,hL) * h_y
    dmu_Y_dx4 = finiteDiff_4D(mu_Y_baseline,3,1,hlogN) + finiteDiff_4D(sigma_Y,3,1,hlogN) * h_y


    dsigma_Y_dx1 = finiteDiff_4D(sigma_Y,0,1,hK)
    dsigma_Y_dx2 = finiteDiff_4D(sigma_Y,1,1,hY)
    dsigma_Y_dx3 = finiteDiff_4D(sigma_Y,2,1,hL)
    dsigma_Y_dx4 = finiteDiff_4D(sigma_Y,3,1,hlogN)

    # mu_Y_baseline = np.sum(theta_ell * pi_c, axis=0) * model_tech1_pre_damage["e_star"]
    # mu_Y_distorted = np.sum(theta_ell * pi_c, axis=0) * model_tech1_pre_damage["e_star"] + sigma_y * model_tech1_pre_damage["h"] * model_tech1_pre_damage["e_star"]
    # sigma_Y = sigma_y* model_tech1_pre_damage["e_star"]


    mu_logIg_baseline = - zeta + psi_0 * (x * np.exp(K_mat - L_mat))**psi_1 - 0.5 * sigma_g**2
    mu_logIg_distorted = - zeta + psi_0 * (x * np.exp(K_mat - L_mat))**psi_1 - 0.5 * sigma_g**2 + sigma_g*h_Ig
    sigma_logIg = sigma_g*np.ones_like(K_mat)


    dmu_logIg_dx1 = finiteDiff_4D(mu_logIg_baseline,0,1,hK) + finiteDiff_4D(sigma_logIg,0,1,hK) * h_Ig
    dmu_logIg_dx2 = finiteDiff_4D(mu_logIg_baseline,1,1,hY) + finiteDiff_4D(sigma_logIg,1,1,hY) * h_Ig
    dmu_logIg_dx3 = finiteDiff_4D(mu_logIg_baseline,2,1,hL) + finiteDiff_4D(sigma_logIg,2,1,hL) * h_Ig
    dmu_logIg_dx4 = finiteDiff_4D(mu_logIg_baseline,3,1,hlogN) + finiteDiff_4D(sigma_logIg,3,1,hlogN) * h_Ig

    dsigma_logIg_dx1 =  finiteDiff_4D(sigma_logIg,0,1,hK)
    dsigma_logIg_dx2 =  finiteDiff_4D(sigma_logIg,1,1,hY)
    dsigma_logIg_dx3 =  finiteDiff_4D(sigma_logIg,2,1,hL)
    dsigma_logIg_dx4 =  finiteDiff_4D(sigma_logIg,3,1,hlogN)

    dG  = gamma_1 + gamma_2 * Y_mat
    ddG = gamma_2 

    mu_logN_baseline = dG * (np.sum(theta_ell * pi_c, axis=0)  ) *  e   + 0.5 * ddG * sigma_y**2 * e **2  
    mu_logN_distorted = dG * (np.sum(theta_ell * pi_c, axis=0)  ) *  e   + 0.5 * ddG * sigma_y**2 * e **2  + dG * e * sigma_y *  h_y

    sigma_logN = dG * e * sigma_y


    dmu_logN_dx1 = finiteDiff_4D(mu_logN_baseline,0,1,hK) + finiteDiff_4D(sigma_logN,0,1,hK) * h_y
    dmu_logN_dx2 = finiteDiff_4D(mu_logN_baseline,1,1,hY) + finiteDiff_4D(sigma_logN,1,1,hY) * h_y
    dmu_logN_dx3 = finiteDiff_4D(mu_logN_baseline,2,1,hL) + finiteDiff_4D(sigma_logN,2,1,hL) * h_y
    dmu_logN_dx4 = finiteDiff_4D(mu_logN_baseline,3,1,hlogN) + finiteDiff_4D(sigma_logN,3,1,hlogN) * h_y

    dsigma_logN_dx1 = finiteDiff_4D(sigma_logN,0,1,hK)
    dsigma_logN_dx2 = finiteDiff_4D(sigma_logN,1,1,hY)
    dsigma_logN_dx3 = finiteDiff_4D(sigma_logN,2,1,hL)
    dsigma_logN_dx4 = finiteDiff_4D(sigma_logN,3,1,hlogN)

    # mu_logIg_baseline = - zeta + psi_0 * (model_tech1_pre_damage["x_star"] * np.exp(K_mat - L_mat))**psi_1 - 0.5 * sigma_g**2
    # mu_logIg_distorted = - zeta + psi_0 * (model_tech1_pre_damage["x_star"] * np.exp(K_mat - L_mat))**psi_1 - 0.5 * sigma_g**2 + sigma_g*model_tech1_pre_damage["h_j"]
    # sigma_logIg = sigma_g*np.ones_like(K_mat)



    mitigation =  alpha * vartheta_bar * (1 - e / (alpha * lambda_bar * np.exp(K_mat)))**theta
    # mitigation =  alpha * vartheta_bar * (1 - model_tech1_pre_damage["e_star"] / (alpha * lambda_bar * np.exp(K_mat)))**theta
    mitigation[mitigation <= 1e-16] = 1e-16

    c = alpha - mitigation - i - x
    # c = alpha - mitigation - model_tech1_pre_damage["i_star"] - model_tech1_pre_damage["x_star"]
    c[c <= 1e-16] = 1e-16

    C = c * np.exp(K_mat)

    U1 = np.log(C)

    J_1 = np.zeros_like(K_mat)
    J_2 = r_1 * (np.exp(r_2 / 2 * (Y_mat - y_bar_lower)**2) -1) * (Y_mat > y_bar_lower)
    J_3 = np.exp( (L_mat - np.log(varrho)) )
    J_4 = np.zeros_like(K_mat)

    V_1 = np.zeros_like(K_mat)
    # V_2 = np.mean(v_i, axis=0)
    V_2 = v_i
    # V_3 = model_tech2_pre_damage["v0"]
    V_3 = v_post_tech
    V_4 = np.zeros_like(K_mat)

    V   = v_pre_damage

    # dG  = gamma_1 + gamma_2 * Y_mat
    # ddG = gamma_2 
    # # U2 = - dG * (np.sum(theta_ell * pi_c, axis=0) ) *  model_tech1_pre_damage["e_star"]   - 0.5 * ddG * sigma_y**2 *  model_tech1_pre_damage["e_star"] **2  
    # U2 = - dG * (np.sum(theta_ell * pi_c, axis=0) + sigma_y *  model_tech1_pre_damage["h"] ) *  model_tech1_pre_damage["e_star"]   - 0.5 * ddG * sigma_y**2 *  model_tech1_pre_damage["e_star"] **2  


    # U3 =  np.exp( (L_mat - np.log(varrho)) )  * V_3 
    # U3_discount = -np.exp(  L_mat - np.log(varrho) )  * V
    # U4 =  J_2*V_2 
    # U4_discount = -J_2*V 

    # U5 =  xi_g * np.exp((L_mat - np.log(varrho))) * (1 - f_tech + f_tech * np.log(f_tech)) + xi_d * Intensity * np.sum( pi_d_o*(1-f_damage+f_damage*np.log(f_damage)),axis=0) 

    # U6 =  1/2 * xi_c * h_y**2 + 1/2 * xi_k * h_k**2 + 1/2 * xi_j * h_Ig**2
    # test_damage = np.abs(J_2*(V_2-V) - U4 - U4_discount)
    # test_tech = np.abs(J_3*(V_3-V) - U3 - U3_discount)
    
    Entropy_x1 = np.zeros_like(K_mat)
    Entropy_x2 = xi_d2 * r_1* r_2 *(Y_mat - y_bar_lower) * np.exp(r_2 / 2 * (Y_mat - y_bar_lower)**2) * (Y_mat>y_bar_lower) * np.sum( pi_d_o*(1-f_damage+f_damage*np.log(f_damage)),axis=0) 
    Entropy_x3 = xi_g2 * np.exp((L_mat - np.log(varrho))) * (1 - f_tech + f_tech * np.log(f_tech)) 
    Entropy_x4 = np.zeros_like(K_mat)


    J_1_x1 = np.zeros_like(K_mat)
    J_2_x2 = r_1* r_2 *(Y_mat - y_bar_lower) * np.exp(r_2 / 2 * (Y_mat - y_bar_lower)**2) * (Y_mat>y_bar_lower)
    J_3_x3 = np.exp((L_mat - np.log(varrho)))
    J_4_x4 = np.zeros_like(K_mat)



    U = U1 - N_mat
    # U = U1 + U2/ delta 

    # U_Full = U1 + U2/ delta + U5/delta+U6/delta
    # U_Full = U1 + U2/ delta + U3/ delta + U4/ delta + U3_discount/ delta + U4_discount/ delta + U5/ delta + U6/ delta
    # U_Full = U1 + U2/ delta + U3/ delta + U4/ delta + U3_discount/ delta + U4_discount/ delta 

    # U_target = U1 

    # U2 -> logN_evolution

    # U3 J2 (tech) -> partialJ + partial V +(discounting) 
    # U4 J3 (damage) -> partialJ + partial V +(discounting) 

    # U5 entropy of jump
    # U6 entropy of diffusion

    res = dict(
        logK = K,
        Y = Y_short,
        # logY = np.log(Y_short),
        logIg = L,
        logN = logN,
        U = U,
        # U_Full = U_Full,
        U1 = U1,
        # U2 = U2,
        # U3_discount = U3_discount,
        # U3 = U3,
        # U4_discount = U4_discount,
        # U4 = U4,
        # U5 = U5,
        # U6 = U6,
        J_1 = J_1,
        J_2 = J_2,
        J_3 = J_3,
        J_4 = J_4,
        J_1_x1 = J_1_x1,
        J_2_x2 = J_2_x2,
        J_3_x3 = J_3_x3,
        J_4_x4 = J_4_x4,
        Entropy_x1 = Entropy_x1,
        Entropy_x2 = Entropy_x2,
        Entropy_x3 = Entropy_x3,        
        Entropy_x4 = Entropy_x4,        
        V_1 = V_1,
        V_2 = V_2,
        V_3 = V_3,
        V_4 = V_4,
        V   = V,
        h_y = h_y,
        h_k = h_k,
        h_Ig = h_Ig,
        f_tech = f_tech,
        f_damage = f_damage,
        mu_logK_baseline = mu_logK_baseline,
        mu_logK_distorted = mu_logK_distorted,
        sigma_logK = sigma_logK,
        mu_Y_baseline = mu_Y_baseline,
        mu_Y_distorted = mu_Y_distorted,
        sigma_Y = sigma_Y,
        # mu_logY_baseline = mu_logY_baseline,
        # mu_logY_distorted = mu_logY_distorted,
        mu_logIg_baseline = mu_logIg_baseline,
        mu_logIg_distorted = mu_logIg_distorted,
        sigma_logIg = sigma_logIg,
        mu_logN_baseline = mu_logN_baseline,
        mu_logN_distorted = mu_logN_distorted,
        sigma_logN = sigma_logN,
        dmu_logK_dx1 = dmu_logK_dx1,
        dmu_logK_dx2 = dmu_logK_dx2,
        dmu_logK_dx3 = dmu_logK_dx3,
        dmu_logK_dx4 = dmu_logK_dx4,
        dmu_Y_dx1 = dmu_Y_dx1,
        dmu_Y_dx2 = dmu_Y_dx2,
        dmu_Y_dx3 = dmu_Y_dx3,
        dmu_Y_dx4 = dmu_Y_dx4,
        dmu_logIg_dx1 = dmu_logIg_dx1,
        dmu_logIg_dx2 = dmu_logIg_dx2,
        dmu_logIg_dx3 = dmu_logIg_dx3,
        dmu_logIg_dx4 = dmu_logIg_dx4,
        dmu_logN_dx1 = dmu_logN_dx1,
        dmu_logN_dx2 = dmu_logN_dx2,
        dmu_logN_dx3 = dmu_logN_dx3,
        dmu_logN_dx4 = dmu_logN_dx4,
        dsigma_logK_dx1 = dsigma_logK_dx1,
        dsigma_logK_dx2 = dsigma_logK_dx2,
        dsigma_logK_dx3 = dsigma_logK_dx3,
        dsigma_logK_dx4 = dsigma_logK_dx4,
        dsigma_Y_dx1 = dsigma_Y_dx1,
        dsigma_Y_dx2 = dsigma_Y_dx2,
        dsigma_Y_dx3 = dsigma_Y_dx3,
        dsigma_Y_dx4 = dsigma_Y_dx4,
        dsigma_logIg_dx1 = dsigma_logIg_dx1,
        dsigma_logIg_dx2 = dsigma_logIg_dx2,
        dsigma_logIg_dx3 = dsigma_logIg_dx3,
        dsigma_logIg_dx4 = dsigma_logIg_dx4,
        dsigma_logN_dx1 = dsigma_logN_dx1,
        dsigma_logN_dx2 = dsigma_logN_dx2,
        dsigma_logN_dx3 = dsigma_logN_dx3,
        dsigma_logN_dx4 = dsigma_logN_dx4,
    )
    print(res.keys())
    with open(Data_Dir + File_Dir2+"model_tech1_pre_damage"+"_prepareID", "wb") as f:
        pickle.dump(res,f)

    
    return res


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):
                for id_rho in range(len(rhoarr)):

                    res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],xia2arr[id_xiag],xik2arr[id_xiag],xic2arr[id_xiag],xij2arr[id_xiag],xid2arr[id_xiag],xig2arr[id_xiag],rhoarr[id_rho],psi0arr[id_psi0],psi1arr[id_psi1],varrhoarr[id_varrho])


