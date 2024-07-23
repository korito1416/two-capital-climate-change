import os
import sys
import numpy as np
np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=200)
from scipy.interpolate import RegularGridInterpolator as RGI
import time
import pickle
from datetime import datetime
from scipy.stats import norm
from tqdm import tqdm
from src.support import *
import argparse



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
parser.add_argument("--seed",type=int,default=1)

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

parser.add_argument("--m0",type=str)

args = parser.parse_args()

seed = args.seed


# Update = args.Update
IntPeriod = args.IntPeriod
dt = 1/12
# dt = 1

print(str(dt))

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
 
m0 = args.m0


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


gamma_3_max = 1/3


y_limit = y_overline

LHS_ylimitlower = gamma_1 * Y_short + gamma_2/2 * Y_short**2 # y<y_limit
LHS_ylimitupper = gamma_1 * Y_short + gamma_2*y_overline * \
    (Y_short-y_limit) + (gamma_2+gamma_3_max)/2 * \
    (Y_short-y_limit)**2 + gamma_2/2 * y_limit**2

LHS_ylimit =LHS_ylimitlower*(Y_short<y_limit) + LHS_ylimitupper*(Y_short>y_limit)

# logN = np.linspace(LHS_ylimit.min(), LHS_ylimit.max(), len(K))
logN = np.linspace(LHS_ylimit.min(), LHS_ylimit.max(), 30)
nlogN = len(logN)
hlogN = logN[1]-logN[0]


delta = 0.01 

(K_mat, Y_mat, L_mat, N_mat) = np.meshgrid(K, Y_short, L, logN, indexing="ij")
# (K_mat, Y_mat, L_mat) = np.meshgrid(K, Y_short, L, indexing="ij")


def simulate_pre(data, Data_Dir, File_Dir):
    print(data.keys())
    W_unique = data['logK']
    Z_unique = data['Y']
    V_unique = data['logIg']
    X_unique = data['logN']

    J1 = data['J_1']
    J2 = data['J_2']
    J3 = data['J_3']
    J4 = data['J_4']

    h_y = data['h_y']
    h_k = data['h_k']
    h_Ig = data['h_Ig']

    f_damage = data['f_damage']
    f_tech = data['f_tech']

    n_damage = len(f_damage)

    print(n_damage)

    V1 = data['V_1']
    V2 = data['V_2']
    V3 = data['V_3']
    V4 = data['V_4']

    Va = data['V']

    wscale = np.diff(W_unique).mean()
    zscale = np.diff(Z_unique).mean()
    vscale = np.diff(V_unique).mean()
    xscale = np.diff(X_unique).mean()

    # df_muW_reshaped = data['mu_logK_baseline']
    # df_muZ_reshaped = data['mu_Y_baseline']
    # df_muV_reshaped = data['mu_logIg_baseline']
    
    df_muW_reshaped = data['mu_logK_distorted']
    df_muZ_reshaped = data['mu_Y_distorted']
    df_muV_reshaped = data['mu_logIg_distorted']
    df_muX_reshaped = data['mu_logN_distorted']

    df_sigmaW0_reshaped = data['sigma_logK']
    df_sigmaZ0_reshaped = data['sigma_Y']
    df_sigmaV0_reshaped = data['sigma_logIg']
    df_sigmaX0_reshaped = data['sigma_logN']

    
    df_U = data['U']
    # df_UF = data['U_Full']
    # df_UF2 = data['U'] + data['U6']/ delta
    # df_U = data['U'] + data['U5']/ delta + data['U6']/ delta
    df_U1 = data['U1']
    # df_U2 = data['U2']
    # df_U3 = data['U3']
    # df_U3_discount = data['U3_discount']
    # df_U4 = data['U4']
    # df_U4_discount = data['U4_discount']

    df_Entropy_dx1 = data['Entropy_x1']
    df_Entropy_dx2 = data['Entropy_x2']
    df_Entropy_dx3 = data['Entropy_x3']
    df_Entropy_dx4 = data['Entropy_x4']


    # df_muW_reshaped_m1 = finiteDiff_4D(df_muW_reshaped,0,1,wscale)
    # df_muW_reshaped_m2 = finiteDiff_4D(df_muW_reshaped,1,1,zscale)
    # df_muW_reshaped_m3 = finiteDiff_4D(df_muW_reshaped,2,1,vscale)
    # df_muZ_reshaped_m1 = finiteDiff_4D(df_muZ_reshaped,0,1,wscale)
    # df_muZ_reshaped_m2 = finiteDiff_4D(df_muZ_reshaped,1,1,zscale)
    # df_muZ_reshaped_m3 = finiteDiff_4D(df_muZ_reshaped,2,1,vscale)
    # df_muV_reshaped_m1 = finiteDiff_4D(df_muV_reshaped,0,1,wscale)
    # df_muV_reshaped_m2 = finiteDiff_4D(df_muV_reshaped,1,1,zscale)
    # df_muV_reshaped_m3 = finiteDiff_4D(df_muV_reshaped,2,1,vscale)
    # df_sigmaW0_reshaped_m1 = finiteDiff_4D(df_sigmaW0_reshaped,0,1,wscale)
    # df_sigmaW0_reshaped_m2 = finiteDiff_4D(df_sigmaW0_reshaped,1,1,zscale)
    # df_sigmaW0_reshaped_m3 = finiteDiff_4D(df_sigmaW0_reshaped,2,1,vscale)
    # df_sigmaZ0_reshaped_m1 = finiteDiff_4D(df_sigmaZ0_reshaped,0,1,wscale)
    # df_sigmaZ0_reshaped_m2 = finiteDiff_4D(df_sigmaZ0_reshaped,1,1,zscale)
    # df_sigmaZ0_reshaped_m3 = finiteDiff_4D(df_sigmaZ0_reshaped,2,1,vscale)
    # df_sigmaV0_reshaped_m1 = finiteDiff_4D(df_sigmaV0_reshaped,0,1,wscale)
    # df_sigmaV0_reshaped_m2 = finiteDiff_4D(df_sigmaV0_reshaped,1,1,zscale)
    # df_sigmaV0_reshaped_m3 = finiteDiff_4D(df_sigmaV0_reshaped,2,1,vscale)


    df_muW_reshaped_m1 = data['dmu_logK_dx1']
    df_muW_reshaped_m2 = data['dmu_logK_dx2']
    df_muW_reshaped_m3 = data['dmu_logK_dx3']
    df_muW_reshaped_m4 = data['dmu_logK_dx4']
    df_muZ_reshaped_m1 = data['dmu_Y_dx1']
    df_muZ_reshaped_m2 = data['dmu_Y_dx2']
    df_muZ_reshaped_m3 = data['dmu_Y_dx3']
    df_muZ_reshaped_m4 = data['dmu_Y_dx4']
    df_muV_reshaped_m1 = data['dmu_logIg_dx1']
    df_muV_reshaped_m2 = data['dmu_logIg_dx2']
    df_muV_reshaped_m3 = data['dmu_logIg_dx3']
    df_muV_reshaped_m4 = data['dmu_logIg_dx4']
    df_muX_reshaped_m1 = data['dmu_logN_dx1']
    df_muX_reshaped_m2 = data['dmu_logN_dx2']
    df_muX_reshaped_m3 = data['dmu_logN_dx3']
    df_muX_reshaped_m4 = data['dmu_logN_dx4']
    df_sigmaW0_reshaped_m1 = data['dsigma_logK_dx1']
    df_sigmaW0_reshaped_m2 = data['dsigma_logK_dx2']
    df_sigmaW0_reshaped_m3 = data['dsigma_logK_dx3']
    df_sigmaW0_reshaped_m4 = data['dsigma_logK_dx4']
    df_sigmaZ0_reshaped_m1 = data['dsigma_Y_dx1']
    df_sigmaZ0_reshaped_m2 = data['dsigma_Y_dx2']
    df_sigmaZ0_reshaped_m3 = data['dsigma_Y_dx3']
    df_sigmaZ0_reshaped_m4 = data['dsigma_Y_dx4']
    df_sigmaV0_reshaped_m1 = data['dsigma_logIg_dx1']
    df_sigmaV0_reshaped_m2 = data['dsigma_logIg_dx2']
    df_sigmaV0_reshaped_m3 = data['dsigma_logIg_dx3']
    df_sigmaV0_reshaped_m4 = data['dsigma_logIg_dx4']
    df_sigmaX0_reshaped_m1 = data['dsigma_logN_dx1']
    df_sigmaX0_reshaped_m2 = data['dsigma_logN_dx2']
    df_sigmaX0_reshaped_m3 = data['dsigma_logN_dx3']
    df_sigmaX0_reshaped_m4 = data['dsigma_logN_dx4']

    df_U_m1 = finiteDiff_4D(df_U,0,1,wscale)
    df_U_m2 = finiteDiff_4D(df_U,1,1,zscale)
    df_U_m3 = finiteDiff_4D(df_U,2,1,vscale)
    df_U_m4 = finiteDiff_4D(df_U,3,1,xscale)

    df_J1_m1 = finiteDiff_4D(J1,0,1,wscale)
    df_J2_m2 = finiteDiff_4D(J2,1,1,zscale)
    df_J3_m3 = finiteDiff_4D(J3,2,1,vscale)
    df_J4_m4 = finiteDiff_4D(J4,3,1,vscale)


    # df_J1_m1 = data['J_1_x1']
    # df_J2_m2 = data['J_2_x2']
    # df_J3_m3 = data['J_3_x3']
    # df_J4_m4 = data['J_4_x4']


    df_V1_m1 = finiteDiff_4D(V1,0,1,wscale)
    df_V1_m2 = finiteDiff_4D(V1,1,1,zscale)
    df_V1_m3 = finiteDiff_4D(V1,2,1,vscale)
    df_V1_m4 = finiteDiff_4D(V1,3,1,xscale)



    df_V2_m1 = np.zeros_like(V2)
    df_V2_m2 = np.zeros_like(V2)
    df_V2_m3 = np.zeros_like(V2)
    df_V2_m4 = np.zeros_like(V2)

    # df_V2_m1 = finiteDiff_4D(V2,0,1,wscale)
    # df_V2_m2 = finiteDiff_4D(V2,1,1,zscale)
    # df_V2_m3 = finiteDiff_4D(V2,2,1,vscale)

    df_f_damage_m1 = np.zeros_like(f_damage)
    df_f_damage_m2 = np.zeros_like(f_damage)
    df_f_damage_m3 = np.zeros_like(f_damage)
    df_f_damage_m4 = np.zeros_like(f_damage)

    for i in range(n_damage):

        
        df_V2_m1_slice = finiteDiff_4D(V2[i,:,:,:,:],0,1,wscale)
        df_V2_m2_slice = finiteDiff_4D(V2[i,:,:,:,:],1,1,zscale)
        df_V2_m3_slice = finiteDiff_4D(V2[i,:,:,:,:],2,1,vscale)
        df_V2_m4_slice = finiteDiff_4D(V2[i,:,:,:,:],3,1,xscale)

        df_V2_m1[i,:,:,:,:] = df_V2_m1_slice
        df_V2_m2[i,:,:,:,:] = df_V2_m2_slice
        df_V2_m3[i,:,:,:,:] = df_V2_m3_slice
        df_V2_m4[i,:,:,:,:] = df_V2_m4_slice


        
        df_f_damage_m1_slice = finiteDiff_4D(f_damage[i,:,:,:,:],0,1,wscale)
        df_f_damage_m2_slice = finiteDiff_4D(f_damage[i,:,:,:,:],1,1,zscale)
        df_f_damage_m3_slice = finiteDiff_4D(f_damage[i,:,:,:,:],2,1,vscale)
        df_f_damage_m4_slice = finiteDiff_4D(f_damage[i,:,:,:,:],3,1,xscale)

        df_f_damage_m1[i,:,:,:,:] = df_f_damage_m1_slice
        df_f_damage_m2[i,:,:,:,:] = df_f_damage_m2_slice
        df_f_damage_m3[i,:,:,:,:] = df_f_damage_m3_slice
        df_f_damage_m4[i,:,:,:,:] = df_f_damage_m4_slice

    df_f_tech_m1 = finiteDiff_4D(f_tech,0,1,wscale)
    df_f_tech_m2 = finiteDiff_4D(f_tech,1,1,zscale)
    df_f_tech_m3 = finiteDiff_4D(f_tech,2,1,vscale)
    df_f_tech_m4 = finiteDiff_4D(f_tech,3,1,xscale)


    df_V3_m1 = finiteDiff_4D(V3,0,1,wscale)
    df_V3_m2 = finiteDiff_4D(V3,1,1,zscale)
    df_V3_m3 = finiteDiff_4D(V3,2,1,vscale)
    df_V3_m4 = finiteDiff_4D(V3,3,1,xscale)

    df_V4_m1 = finiteDiff_4D(V4,0,1,wscale)
    df_V4_m2 = finiteDiff_4D(V4,1,1,zscale)
    df_V4_m3 = finiteDiff_4D(V4,2,1,vscale)
    df_V4_m4 = finiteDiff_4D(V4,3,1,xscale)

    df_V_m1 = finiteDiff_4D(Va,0,1,wscale)
    df_V_m2 = finiteDiff_4D(Va,1,1,zscale)
    df_V_m3 = finiteDiff_4D(Va,2,1,vscale)
    df_V_m4 = finiteDiff_4D(Va,3,1,xscale)

    df_muW_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], df_muW_reshaped, fill_value=None, bounds_error=True)
    df_muZ_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], df_muZ_reshaped, fill_value=None, bounds_error=True)
    df_muV_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], df_muV_reshaped, fill_value=None, bounds_error=True)
    df_muX_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], df_muX_reshaped, fill_value=None, bounds_error=True)

    df_muW_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muW_reshaped_m1, fill_value=None, bounds_error=True)
    df_muW_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muW_reshaped_m2, fill_value=None, bounds_error=True)
    df_muW_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muW_reshaped_m3, fill_value=None, bounds_error=True)
    df_muW_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muW_reshaped_m4, fill_value=None, bounds_error=True)

    df_muZ_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muZ_reshaped_m1, fill_value=None, bounds_error=True)
    df_muZ_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muZ_reshaped_m2, fill_value=None, bounds_error=True)
    df_muZ_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muZ_reshaped_m3, fill_value=None, bounds_error=True)
    df_muZ_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muZ_reshaped_m4, fill_value=None, bounds_error=True)

    df_muV_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muV_reshaped_m1, fill_value=None, bounds_error=True)
    df_muV_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muV_reshaped_m2, fill_value=None, bounds_error=True)
    df_muV_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muV_reshaped_m3, fill_value=None, bounds_error=True)
    df_muV_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muV_reshaped_m4, fill_value=None, bounds_error=True)

    df_muX_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muX_reshaped_m1, fill_value=None, bounds_error=True)
    df_muX_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muX_reshaped_m2, fill_value=None, bounds_error=True)
    df_muX_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muX_reshaped_m3, fill_value=None, bounds_error=True)
    df_muX_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_muX_reshaped_m4, fill_value=None, bounds_error=True)

    df_sigmaW0_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaW0_reshaped, fill_value=None, bounds_error=True)
    df_sigmaZ0_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaZ0_reshaped, fill_value=None, bounds_error=True)
    df_sigmaV0_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaV0_reshaped, fill_value=None, bounds_error=True)
    df_sigmaX0_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaX0_reshaped, fill_value=None, bounds_error=True)

    df_sigmaW0_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaW0_reshaped_m1, fill_value=None, bounds_error=True)
    df_sigmaW0_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaW0_reshaped_m2, fill_value=None, bounds_error=True)
    df_sigmaW0_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaW0_reshaped_m3, fill_value=None, bounds_error=True)
    df_sigmaW0_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaW0_reshaped_m4, fill_value=None, bounds_error=True)

    df_sigmaZ0_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaZ0_reshaped_m1, fill_value=None, bounds_error=True)
    df_sigmaZ0_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaZ0_reshaped_m2, fill_value=None, bounds_error=True)
    df_sigmaZ0_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaZ0_reshaped_m3, fill_value=None, bounds_error=True)
    df_sigmaZ0_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaZ0_reshaped_m4, fill_value=None, bounds_error=True)

    df_sigmaV0_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaV0_reshaped_m1, fill_value=None, bounds_error=True)
    df_sigmaV0_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaV0_reshaped_m2, fill_value=None, bounds_error=True)
    df_sigmaV0_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaV0_reshaped_m3, fill_value=None, bounds_error=True)
    df_sigmaV0_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaV0_reshaped_m4, fill_value=None, bounds_error=True)

    df_sigmaX0_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaX0_reshaped_m1, fill_value=None, bounds_error=True)
    df_sigmaX0_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaX0_reshaped_m2, fill_value=None, bounds_error=True)
    df_sigmaX0_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaX0_reshaped_m3, fill_value=None, bounds_error=True)
    df_sigmaX0_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_sigmaX0_reshaped_m4, fill_value=None, bounds_error=True)

    df_U_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U_m1, fill_value=None, bounds_error=True)
    df_U_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U_m2, fill_value=None, bounds_error=True)
    df_U_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U_m3, fill_value=None, bounds_error=True)
    df_U_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U_m4, fill_value=None, bounds_error=True)

    # df_U1_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U1_m1, fill_value=None, bounds_error=True)
    # df_U1_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U1_m2, fill_value=None, bounds_error=True)
    # df_U1_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U1_m3, fill_value=None, bounds_error=True)

    # df_U2_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U2_m1, fill_value=None, bounds_error=True)
    # df_U2_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U2_m2, fill_value=None, bounds_error=True)
    # df_U2_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U2_m3, fill_value=None, bounds_error=True)

    df_Entropy_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_Entropy_dx1, fill_value=None, bounds_error=True)
    df_Entropy_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_Entropy_dx2, fill_value=None, bounds_error=True)
    df_Entropy_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_Entropy_dx3, fill_value=None, bounds_error=True)
    df_Entropy_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_Entropy_dx4, fill_value=None, bounds_error=True)


    # df_U3_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U3_m1, fill_value=None, bounds_error=True)
    # df_U3_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U3_m2, fill_value=None, bounds_error=True)
    # df_U3_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U3_m3, fill_value=None, bounds_error=True)

    # df_U4_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U4_m1, fill_value=None, bounds_error=True)
    # df_U4_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U4_m2, fill_value=None, bounds_error=True)
    # df_U4_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_U4_m3, fill_value=None, bounds_error=True)

    # df_UF_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_UF_m1, fill_value=None, bounds_error=True)
    # df_UF_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_UF_m2, fill_value=None, bounds_error=True)
    # df_UF_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_UF_m3, fill_value=None, bounds_error=True)

    # df_UF2_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_UF2_m1, fill_value=None, bounds_error=True)
    # df_UF2_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_UF2_m2, fill_value=None, bounds_error=True)
    # df_UF2_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_UF2_m3, fill_value=None, bounds_error=True)

    df_J1_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_J1_m1, fill_value=None, bounds_error=True)
    df_J2_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_J2_m2, fill_value=None, bounds_error=True)
    df_J3_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_J3_m3, fill_value=None, bounds_error=True)
    df_J4_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_J4_m4, fill_value=None, bounds_error=True)
    
    df_V1_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V1_m1, fill_value=None, bounds_error=True)
    df_V1_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V1_m2, fill_value=None, bounds_error=True)
    df_V1_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V1_m3, fill_value=None, bounds_error=True)
    df_V1_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V1_m4, fill_value=None, bounds_error=True)


    df_V2_interpolated_m1 = []
    df_V2_interpolated_m2 = []
    df_V2_interpolated_m3 = []
    df_V2_interpolated_m4 = []

    df_f_damage_interpolated_m1 = []
    df_f_damage_interpolated_m2 = []
    df_f_damage_interpolated_m3 = []
    df_f_damage_interpolated_m4 = []

    f_damage_interpolated=[]

    for i in range(n_damage):

        df_V2_interpolated_m1_slice = RGI([W_unique,Z_unique, V_unique, X_unique], df_V2_m1[i,:,:,:,:], fill_value=None, bounds_error=True)
        df_V2_interpolated_m2_slice = RGI([W_unique,Z_unique, V_unique, X_unique], df_V2_m2[i,:,:,:,:], fill_value=None, bounds_error=True)
        df_V2_interpolated_m3_slice = RGI([W_unique,Z_unique, V_unique, X_unique], df_V2_m3[i,:,:,:,:], fill_value=None, bounds_error=True)
        df_V2_interpolated_m4_slice = RGI([W_unique,Z_unique, V_unique, X_unique], df_V2_m4[i,:,:,:,:], fill_value=None, bounds_error=True)

        df_f_damage_interpolated_m1_slice = RGI([W_unique,Z_unique, V_unique, X_unique], df_f_damage_m1[i,:,:,:,:], fill_value=None, bounds_error=True)
        df_f_damage_interpolated_m2_slice = RGI([W_unique,Z_unique, V_unique, X_unique], df_f_damage_m2[i,:,:,:,:], fill_value=None, bounds_error=True)
        df_f_damage_interpolated_m3_slice = RGI([W_unique,Z_unique, V_unique, X_unique], df_f_damage_m3[i,:,:,:,:], fill_value=None, bounds_error=True)
        df_f_damage_interpolated_m4_slice = RGI([W_unique,Z_unique, V_unique, X_unique], df_f_damage_m4[i,:,:,:,:], fill_value=None, bounds_error=True)

        f_damage_interpolated_slice = RGI([W_unique,Z_unique, V_unique, X_unique], f_damage[i,:,:,:,:], fill_value=None, bounds_error=True)

        df_V2_interpolated_m1.append(df_V2_interpolated_m1_slice)
        df_V2_interpolated_m2.append(df_V2_interpolated_m2_slice)
        df_V2_interpolated_m3.append(df_V2_interpolated_m3_slice)
        df_V2_interpolated_m4.append(df_V2_interpolated_m4_slice)

        df_f_damage_interpolated_m1.append(df_f_damage_interpolated_m1_slice)
        df_f_damage_interpolated_m2.append(df_f_damage_interpolated_m2_slice)
        df_f_damage_interpolated_m3.append(df_f_damage_interpolated_m3_slice)
        df_f_damage_interpolated_m4.append(df_f_damage_interpolated_m4_slice)

        f_damage_interpolated.append(f_damage_interpolated_slice)

    df_f_tech_m1_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], df_f_tech_m1, fill_value=None, bounds_error=True)
    df_f_tech_m2_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], df_f_tech_m2, fill_value=None, bounds_error=True)
    df_f_tech_m3_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], df_f_tech_m3, fill_value=None, bounds_error=True)
    df_f_tech_m4_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], df_f_tech_m4, fill_value=None, bounds_error=True)

    f_tech_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], f_tech, fill_value=None, bounds_error=True)

    df_V3_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V3_m1, fill_value=None, bounds_error=True)
    df_V3_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V3_m2, fill_value=None, bounds_error=True)
    df_V3_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V3_m3, fill_value=None, bounds_error=True)
    df_V3_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V3_m4, fill_value=None, bounds_error=True)

    df_V4_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V4_m1, fill_value=None, bounds_error=True)
    df_V4_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V4_m2, fill_value=None, bounds_error=True)
    df_V4_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V4_m3, fill_value=None, bounds_error=True)
    df_V4_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V4_m4, fill_value=None, bounds_error=True)

    df_V_interpolated_m1 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V_m1, fill_value=None, bounds_error=True)
    df_V_interpolated_m2 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V_m2, fill_value=None, bounds_error=True)
    df_V_interpolated_m3 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V_m3, fill_value=None, bounds_error=True)
    df_V_interpolated_m4 = RGI([W_unique,Z_unique, V_unique, X_unique], df_V_m4, fill_value=None, bounds_error=True)

    df_V_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], Va, fill_value=None, bounds_error=True)

    df_J1_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], J1, fill_value=None, bounds_error=True)
    df_J2_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], J2, fill_value=None, bounds_error=True)
    df_J3_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], J3, fill_value=None, bounds_error=True)
    df_J4_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], J4, fill_value=None, bounds_error=True)

    df_V1_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], V1, fill_value=None, bounds_error=True)


    df_V2_interpolated = []

    for i in range(n_damage):

        df_V2_interpolated_slice = RGI([W_unique,Z_unique, V_unique, X_unique], V2[i,:,:,:,:], fill_value=None, bounds_error=True)

        df_V2_interpolated.append(df_V2_interpolated_slice)


    # df_V2_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], V2, fill_value=None, bounds_error=True)


    df_V3_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], V3, fill_value=None, bounds_error=True)
    df_V4_interpolated = RGI([W_unique,Z_unique, V_unique, X_unique], V4, fill_value=None, bounds_error=True)


    w_processes = []
    z_processes = []
    v_processes = []
    x_processes = []
    m1_processes = []
    m2_processes = []
    m3_processes = []
    m4_processes = []
    u1_processes = []
    u2_processes = []
    u3_processes = []
    u4_processes = []
    
    # u11_processes = []
    # u12_processes = []
    # u13_processes = []
    
    # u21_processes = []
    # u22_processes = []
    # u23_processes = []


    # u31_processes = []
    # u32_processes = []
    # u33_processes = []
    
    # u41_processes = []
    # u42_processes = []
    # u43_processes = []

    # uF1_processes = []
    # uF2_processes = []
    # uF3_processes = []

    # uF21_processes = []
    # uF22_processes = []
    # uF23_processes = []


    j1_processes = []
    j2_processes = []
    j3_processes = []
    j4_processes = []

    v1_processes = []
    v2_processes = []
    v3_processes = []
    v4_processes = []

    va_processes = []

    dj1_dx1_processes = []
    dj2_dx2_processes = []
    dj3_dx3_processes = []
    dj4_dx4_processes = []

    dv1_dx1_processes = []
    dv1_dx2_processes = []
    dv1_dx3_processes = []
    dv1_dx4_processes = []
    
    dv2_dx1_processes = []
    dv2_dx2_processes = []
    dv2_dx3_processes = []
    dv2_dx4_processes = []

    dv3_dx1_processes = []
    dv3_dx2_processes = []
    dv3_dx3_processes = []
    dv3_dx4_processes = []
    
    dv4_dx1_processes = []
    dv4_dx2_processes = []
    dv4_dx3_processes = []
    dv4_dx4_processes = []
    
    dv_dx1_processes = []
    dv_dx2_processes = []
    dv_dx3_processes = []
    dv_dx4_processes = []

    dfdamage_dx1_processes = []
    dfdamage_dx2_processes = []
    dfdamage_dx3_processes = []
    dfdamage_dx4_processes = []

    dftech_dx1_processes = []
    dftech_dx2_processes = []
    dftech_dx3_processes = []
    dftech_dx4_processes = []

    damage_processes = []
    tech_processes = []

    Entropy_x1_processes = []
    Entropy_x2_processes = []
    Entropy_x3_processes = []
    Entropy_x4_processes = []

    first_terms = []
    first_term_fs = []
    second_terms = []
    second_term_1s = []
    second_term_2s = []
    second_term_3s = []
    second_term_4s = []
    third_term_1s = []
    third_term_2s = []
    third_term_3s = []
    third_term_4s = []
    fourth_term_entropys = []
    fourth_term_2s = []
    fourth_term_3s = []
    undiscount_processes = []
    undiscount_o_processes = []
    undiscount_pg_processes = []
    undiscount_pg_f2_processes = []
    undiscount_f_processes = []
    undiscount_f_pg_processes = []
    discount_factor_temp_processes = []
    discount_factor_processes = []
    discount_processes = []
    discount_o_processes = []
    discount_pg_processes = []
    discount_pg_f2_processes = []
    discount_f_processes = []
    discount_f_pg_processes = []

    discount_factor_nodelta_DisSep_Damage_processes = []
    discount_factor_nodelta_DisSep_Tech_processes = []
    discount_factor_nodelta_processes = []


    discount_factor_nodeltadt_DisSep_Damage_processes = []
    discount_factor_nodeltadt_DisSep_Tech_processes= []
    discount_factor_nodeltadt_processes = []

    first_terms2 = []
    second_terms2 = []
    third_term_1s2 = []
    third_term_2s2 = []
    third_term_3s2 = []
    undiscount_processes2 = []
    discount_factor_temp_processes2 = []
    discount_factor_processes2 = []
    discount_processes2 = []
    
    time_processes = []
    N = int(IntPeriod/dt)
    # for i in tqdm(range(50)):
    for i in tqdm(range(100)):
        
        # i = seed + i*100
        # np.random.seed(i)
        print('Seed:', i)

        try:

        # N = IntPeriod

            W1 = np.random.normal(0,1,N)*np.sqrt(dt)
            W2 = np.random.normal(0,1,N)*np.sqrt(dt)
            W3 = np.random.normal(0,1,N)*np.sqrt(dt)
            W4 = np.random.normal(0,1,N)*np.sqrt(dt)

            w_process = np.zeros(N)
            z_process = np.zeros(N)
            v_process = np.zeros(N)
            x_process = np.zeros(N)


            m1_process = np.zeros(N)
            m2_process = np.zeros(N)
            m3_process = np.zeros(N)
            m4_process = np.zeros(N)

            u1_process = np.zeros(N)
            u2_process = np.zeros(N)
            u3_process = np.zeros(N)
            u4_process = np.zeros(N)
            
            # u11_process = np.zeros(N)
            # u12_process = np.zeros(N)
            # u13_process = np.zeros(N)
            
            # u21_process = np.zeros(N)
            # u22_process = np.zeros(N)
            # u23_process = np.zeros(N)


            # u31_process = np.zeros(N)
            # u32_process = np.zeros(N)
            # u33_process = np.zeros(N)
            
            # u41_process = np.zeros(N)
            # u42_process = np.zeros(N)
            # u43_process = np.zeros(N)


            # uF1_process = np.zeros(N)
            # uF2_process = np.zeros(N)
            # uF3_process = np.zeros(N)

            # uF21_process = np.zeros(N)
            # uF22_process = np.zeros(N)
            # uF23_process = np.zeros(N)

            j1_process = np.zeros(N)
            j2_process = np.zeros(N)
            j3_process = np.zeros(N)
            j4_process = np.zeros(N)
            
            v1_process = np.zeros(N)
            # v2_process = np.zeros(N)
            v2_process = np.zeros([n_damage,N])
            v3_process = np.zeros(N)
            v4_process = np.zeros(N)

            va_process = np.zeros(N)  

            dj1_dx1_process = np.zeros(N)
            dj2_dx2_process = np.zeros(N)
            dj3_dx3_process = np.zeros(N)
            dj4_dx4_process = np.zeros(N)
            
            dv1_dx1_process = np.zeros(N)
            dv1_dx2_process = np.zeros(N)
            dv1_dx3_process = np.zeros(N)
            dv1_dx4_process = np.zeros(N)

            # dv2_dx1_process = np.zeros(N)
            # dv2_dx2_process = np.zeros(N)
            # dv2_dx3_process = np.zeros(N)

            dv2_dx1_process = np.zeros([n_damage,N])
            dv2_dx2_process = np.zeros([n_damage,N])
            dv2_dx3_process = np.zeros([n_damage,N])
            dv2_dx4_process = np.zeros([n_damage,N])

            dfdamage_dx1_process = np.zeros([n_damage,N])
            dfdamage_dx2_process = np.zeros([n_damage,N])
            dfdamage_dx3_process = np.zeros([n_damage,N])
            dfdamage_dx4_process = np.zeros([n_damage,N])

            dftech_dx1_process = np.zeros(N)
            dftech_dx2_process = np.zeros(N)
            dftech_dx3_process = np.zeros(N)
            dftech_dx4_process = np.zeros(N)

            damage_process = np.zeros([n_damage,N])
            tech_process = np.zeros(N)

            Entropy_x1_process = np.zeros(N)
            Entropy_x2_process = np.zeros(N)
            Entropy_x3_process = np.zeros(N)
            Entropy_x4_process = np.zeros(N)

            dv3_dx1_process = np.zeros(N)
            dv3_dx2_process = np.zeros(N)
            dv3_dx3_process = np.zeros(N)
            dv3_dx4_process = np.zeros(N)

            dv4_dx1_process = np.zeros(N)
            dv4_dx2_process = np.zeros(N)
            dv4_dx3_process = np.zeros(N)
            dv4_dx4_process = np.zeros(N)

            dv_dx1_process = np.zeros(N)
            dv_dx2_process = np.zeros(N)
            dv_dx3_process = np.zeros(N)
            dv_dx4_process = np.zeros(N)

            first_term = np.zeros(N)
            first_term_f = np.zeros(N)
            first_term_f2 = np.zeros(N)
            second_term = np.zeros(N)
            second_term_1 = np.zeros(N)
            second_term_2 = np.zeros(N)
            second_term_3 = np.zeros(N)
            second_term_4 = np.zeros(N)            
            third_term_1 = np.zeros(N)
            third_term_2 = np.zeros(N)
            third_term_3 = np.zeros(N)
            third_term_4 = np.zeros(N)

            fourth_term_entropy = np.zeros(N)
            fourth_term_2 = np.zeros(N)
            fourth_term_3 = np.zeros(N)

            undiscount_process = np.zeros(N)
            undiscount_o_process = np.zeros(N)
            undiscount_pg_process = np.zeros(N)
            undiscount_pg_f2_process = np.zeros(N)
            undiscount_f_process = np.zeros(N)
            undiscount_f_pg_process = np.zeros(N)
            discount_factors = np.zeros(N)
            discount_factor_temps = np.zeros(N)
            discount_factor_temps_nodelta = np.zeros(N)
            discount_factor_temps_nodeltadt = np.zeros(N)

            discount_factor_nodelta_DisSep_Damage = np.zeros(N)
            discount_factor_nodelta_DisSep_Tech = np.zeros(N)
            discount_factor_nodelta = np.zeros(N)

            

            discount_factor_nodeltadt_DisSep_Damage = np.zeros(N)
            discount_factor_nodeltadt_DisSep_Tech = np.zeros(N)
            discount_factor_nodeltadt = np.zeros(N)

            discount_process = np.zeros(N)
            discount_o_process = np.zeros(N)
            discount_pg_process = np.zeros(N)
            discount_pg_f2_process = np.zeros(N)
            discount_f_process = np.zeros(N)
            discount_f_pg_process = np.zeros(N)

            first_term2 = np.zeros(N)
            second_term2 = np.zeros(N)

            third_term_12 = np.zeros(N)
            third_term_22 = np.zeros(N)
            third_term_32 = np.zeros(N)
            
            undiscount_process2 = np.zeros(N)
            discount_factors2 = np.zeros(N)
            discount_factor_temps2 = np.zeros(N)
            discount_process2 = np.zeros(N)
            time_process = np.zeros(N)

            w_process[0] = np.log(85/0.115)
            z_process[0] = 1.1
            v_process[0] = np.log(11.2)
            x_process[0] = gamma_1 * z_process[0] + gamma_2/2 * z_process[0]**2 # y<y_limit

            if m0=="Temperature":
                m1_process[0] = 0
                m2_process[0] = 1
                m3_process[0] = 0
                m4_process[0] = 0


            if m0=="Technology":
                m1_process[0] = 0
                m2_process[0] = 0
                m3_process[0] = 1
                m4_process[0] = 0

            u1_process[0] = df_U_interpolated_m1([w_process[0],z_process[0],v_process[0],x_process[0]])
            u2_process[0] = df_U_interpolated_m2([w_process[0],z_process[0],v_process[0],x_process[0]])
            u3_process[0] = df_U_interpolated_m3([w_process[0],z_process[0],v_process[0],x_process[0]])
            u4_process[0] = df_U_interpolated_m4([w_process[0],z_process[0],v_process[0],x_process[0]])
            

            # uF1_process[0] = df_UF_interpolated_m1([w_process[0],z_process[0],v_process[0],x_process[0]])
            # uF2_process[0] = df_UF_interpolated_m2([w_process[0],z_process[0],v_process[0],x_process[0]])
            # uF3_process[0] = df_UF_interpolated_m3([w_process[0],z_process[0],v_process[0],x_process[0]])
            
            # uF21_process[0] = df_UF2_interpolated_m1([w_process[0],z_process[0],v_process[0],x_process[0]])
            # uF22_process[0] = df_UF2_interpolated_m2([w_process[0],z_process[0],v_process[0],x_process[0]])
            # uF23_process[0] = df_UF2_interpolated_m3([w_process[0],z_process[0],v_process[0],x_process[0]])
            
            j1_process[0] = df_J1_interpolated([w_process[0],z_process[0],v_process[0],x_process[0]])
            j2_process[0] = df_J2_interpolated([w_process[0],z_process[0],v_process[0],x_process[0]])
            j3_process[0] = df_J3_interpolated([w_process[0],z_process[0],v_process[0],x_process[0]])
            j4_process[0] = df_J4_interpolated([w_process[0],z_process[0],v_process[0],x_process[0]])

            v1_process[0] = df_V1_interpolated([w_process[0],z_process[0],v_process[0],x_process[0]])
            # v2_process[0] = df_V2_interpolated([w_process[0],z_process[0],v_process[0],x_process[0]])


            for i in range(n_damage):

                v2_process[i,0] = df_V2_interpolated[i]([w_process[0],z_process[0],v_process[0],x_process[0]])
                damage_process[i,0] = f_damage_interpolated[i]([w_process[0],z_process[0],v_process[0],x_process[0]])
                dfdamage_dx1_process[i,0] = df_f_damage_interpolated_m1[i]([w_process[0],z_process[0],v_process[0],x_process[0]])
                dfdamage_dx2_process[i,0] = df_f_damage_interpolated_m2[i]([w_process[0],z_process[0],v_process[0],x_process[0]])
                dfdamage_dx3_process[i,0] = df_f_damage_interpolated_m3[i]([w_process[0],z_process[0],v_process[0],x_process[0]])
                dfdamage_dx4_process[i,0] = df_f_damage_interpolated_m4[i]([w_process[0],z_process[0],v_process[0],x_process[0]])

            dftech_dx1_process[0] = df_f_tech_m1_interpolated([w_process[0],z_process[0],v_process[0],x_process[0]])
            dftech_dx2_process[0] = df_f_tech_m2_interpolated([w_process[0],z_process[0],v_process[0],x_process[0]])
            dftech_dx3_process[0] = df_f_tech_m3_interpolated([w_process[0],z_process[0],v_process[0],x_process[0]])
            dftech_dx4_process[0] = df_f_tech_m4_interpolated([w_process[0],z_process[0],v_process[0],x_process[0]])

            tech_process[0] = f_tech_interpolated([w_process[0],z_process[0],v_process[0],x_process[0]])

            Entropy_x1_process[0] = df_Entropy_interpolated_m1([w_process[0],z_process[0],v_process[0],x_process[0]])
            Entropy_x2_process[0] = df_Entropy_interpolated_m2([w_process[0],z_process[0],v_process[0],x_process[0]])
            Entropy_x3_process[0] = df_Entropy_interpolated_m3([w_process[0],z_process[0],v_process[0],x_process[0]])
            Entropy_x4_process[0] = df_Entropy_interpolated_m4([w_process[0],z_process[0],v_process[0],x_process[0]])

            v3_process[0] = df_V3_interpolated([w_process[0],z_process[0],v_process[0],x_process[0]])
            v4_process[0] = df_V4_interpolated([w_process[0],z_process[0],v_process[0],x_process[0]])


            va_process[0] = df_V_interpolated([w_process[0],z_process[0],v_process[0],x_process[0]])

            dj1_dx1_process[0] = df_J1_interpolated_m1([w_process[0],z_process[0],v_process[0],x_process[0]])
            dj2_dx2_process[0] = df_J2_interpolated_m2([w_process[0],z_process[0],v_process[0],x_process[0]])
            dj3_dx3_process[0] = df_J3_interpolated_m3([w_process[0],z_process[0],v_process[0],x_process[0]])
            dj4_dx4_process[0] = df_J4_interpolated_m4([w_process[0],z_process[0],v_process[0],x_process[0]])
            
            dv1_dx1_process[0] = df_V1_interpolated_m1([w_process[0],z_process[0],v_process[0],x_process[0]])
            dv1_dx2_process[0] = df_V1_interpolated_m2([w_process[0],z_process[0],v_process[0],x_process[0]])
            dv1_dx3_process[0] = df_V1_interpolated_m3([w_process[0],z_process[0],v_process[0],x_process[0]])
            dv1_dx4_process[0] = df_V1_interpolated_m4([w_process[0],z_process[0],v_process[0],x_process[0]])

            # dv2_dx1_process[0] = df_V2_interpolated_m1([w_process[0],z_process[0],v_process[0],x_process[0]])
            # dv2_dx2_process[0] = df_V2_interpolated_m2([w_process[0],z_process[0],v_process[0],x_process[0]])
            # dv2_dx3_process[0] = df_V2_interpolated_m3([w_process[0],z_process[0],v_process[0],x_process[0]])

            for i in range(n_damage):
                dv2_dx1_process[i,0] = df_V2_interpolated_m1[i]([w_process[0],z_process[0],v_process[0],x_process[0]])
                dv2_dx2_process[i,0] = df_V2_interpolated_m2[i]([w_process[0],z_process[0],v_process[0],x_process[0]])
                dv2_dx3_process[i,0] = df_V2_interpolated_m3[i]([w_process[0],z_process[0],v_process[0],x_process[0]])
                dv2_dx4_process[i,0] = df_V2_interpolated_m4[i]([w_process[0],z_process[0],v_process[0],x_process[0]])



            dv3_dx1_process[0] = df_V3_interpolated_m1([w_process[0],z_process[0],v_process[0],x_process[0]])
            dv3_dx2_process[0] = df_V3_interpolated_m2([w_process[0],z_process[0],v_process[0],x_process[0]])
            dv3_dx3_process[0] = df_V3_interpolated_m3([w_process[0],z_process[0],v_process[0],x_process[0]])
            dv3_dx4_process[0] = df_V3_interpolated_m4([w_process[0],z_process[0],v_process[0],x_process[0]])


            dv4_dx1_process[0] = df_V4_interpolated_m1([w_process[0],z_process[0],v_process[0],x_process[0]])
            dv4_dx2_process[0] = df_V4_interpolated_m2([w_process[0],z_process[0],v_process[0],x_process[0]])
            dv4_dx3_process[0] = df_V4_interpolated_m3([w_process[0],z_process[0],v_process[0],x_process[0]])
            dv4_dx4_process[0] = df_V4_interpolated_m4([w_process[0],z_process[0],v_process[0],x_process[0]])

            dv_dx1_process[0] = df_V_interpolated_m1([w_process[0],z_process[0],v_process[0],x_process[0]])
            dv_dx2_process[0] = df_V_interpolated_m2([w_process[0],z_process[0],v_process[0],x_process[0]])
            dv_dx3_process[0] = df_V_interpolated_m3([w_process[0],z_process[0],v_process[0],x_process[0]])
            dv_dx4_process[0] = df_V_interpolated_m4([w_process[0],z_process[0],v_process[0],x_process[0]])
            
            first_term[0] = u1_process[0]*m1_process[0]+u2_process[0]*m2_process[0]+u3_process[0]*m3_process[0]+u4_process[0]*m4_process[0]
            # first_term_f[0] = uF1_process[0]*m1_process[0]+uF2_process[0]*m2_process[0]+uF3_process[0]*m3_process[0]
            # first_term_f2[0] = uF21_process[0]*m1_process[0]+uF22_process[0]*m2_process[0]+uF23_process[0]*m3_process[0]
            # second_term[0] = dj1_dx1_process[0]*(v1_process[0] - va_process[0])*m1_process[0]+\
            #                 dj2_dx2_process[0]*(v2_process[0] - va_process[0])*m2_process[0]+\
            #                 dj3_dx3_process[0]*(v3_process[0] - va_process[0])*m3_process[0]
            second_term[0] = dj1_dx1_process[0]*(v1_process[0] - va_process[0])*m1_process[0]+\
                            dj2_dx2_process[0]*np.mean(damage_process[:,0]*(v2_process[:,0] - va_process[0]),axis=0)*m2_process[0]+\
                            dj3_dx3_process[0]*tech_process[0]*(v3_process[0] - va_process[0])*m3_process[0]+\
                            dj4_dx4_process[0]*(v4_process[0] - va_process[0])*m4_process[0]
            second_term_1[0] = dj1_dx1_process[0]*(v1_process[0] - va_process[0])*m1_process[0]
            second_term_2[0] = dj2_dx2_process[0]*np.mean(damage_process[:,0]*(v2_process[:,0] - va_process[0]),axis=0)*m2_process[0]
            second_term_3[0] = dj3_dx3_process[0]*tech_process[0]*(v3_process[0] - va_process[0])*m3_process[0]
            second_term_4[0] = dj4_dx4_process[0]*(v4_process[0] - va_process[0])*m4_process[0]

            third_term_1[0] = j1_process[0]*(dv1_dx1_process[0])*m1_process[0]+\
                            j1_process[0]*(dv1_dx2_process[0])*m2_process[0]+\
                            j1_process[0]*(dv1_dx3_process[0])*m3_process[0]+\
                            j1_process[0]*(dv1_dx4_process[0])*m4_process[0]
            # third_term_2[0] = j2_process[0]*(dv2_dx1_process[0])*m1_process[0]+\
            #                 j2_process[0]*(dv2_dx2_process[0])*m2_process[0]+\
            #                 j2_process[0]*(dv2_dx3_process[0])*m3_process[0]
            third_term_2[0] = j2_process[0]*np.mean(damage_process[:,0]*dv2_dx1_process[:,0],axis=0)*m1_process[0]+\
                            j2_process[0]*np.mean(damage_process[:,0]*dv2_dx2_process[:,0],axis=0)*m2_process[0]+\
                            j2_process[0]*np.mean(damage_process[:,0]*dv2_dx3_process[:,0],axis=0)*m3_process[0]+\
                            j2_process[0]*np.mean(damage_process[:,0]*dv2_dx4_process[:,0],axis=0)*m4_process[0]
            # third_term_3[0] = j3_process[0]*(dv3_dx1_process[0])*m1_process[0]+\
            #                 j3_process[0]*(dv3_dx2_process[0])*m2_process[0]+\
            #                 j3_process[0]*(dv3_dx3_process[0])*m3_process[0]
            third_term_3[0] = j3_process[0]*tech_process[0]*(dv3_dx1_process[0])*m1_process[0]+\
                            j3_process[0]*tech_process[0]*(dv3_dx2_process[0])*m2_process[0]+\
                            j3_process[0]*tech_process[0]*(dv3_dx3_process[0])*m3_process[0]+\
                            j3_process[0]*tech_process[0]*(dv3_dx4_process[0])*m4_process[0]
            third_term_4[0] = j4_process[0]*(dv4_dx1_process[0])*m1_process[0]+\
                            j4_process[0]*(dv4_dx2_process[0])*m2_process[0]+\
                            j4_process[0]*(dv4_dx3_process[0])*m3_process[0]+\
                            j4_process[0]*(dv4_dx4_process[0])*m4_process[0]

            fourth_term_entropy[0] = Entropy_x1_process[0]*m1_process[0]+\
                                    Entropy_x2_process[0]*m2_process[0]+\
                                    Entropy_x3_process[0]*m3_process[0]+\
                                    Entropy_x4_process[0]*m4_process[0]
            time_process[0] = 0
            # fourth_term_2[0] = j2_process[0]*np.mean(dfdamage_dx1_process[:,0]*(v2_process[:,0] - va_process[0]),axis=0)*m1_process[0]+\
            #                 j2_process[0]*np.mean(dfdamage_dx2_process[:,0]*(v2_process[:,0] - va_process[0]),axis=0)*m2_process[0]+\
            #                 j2_process[0]*np.mean(dfdamage_dx3_process[:,0]*(v2_process[:,0] - va_process[0]),axis=0)*m3_process[0]
            # fourth_term_3[0] = j3_process[0]*dftech_dx1_process[0]*(v3_process[0] - va_process[0])*m1_process[0]+\
            #                 j3_process[0]*dftech_dx2_process[0]*(v3_process[0] - va_process[0])*m2_process[0]+\
            #                 j3_process[0]*dftech_dx3_process[0]*(v3_process[0] - va_process[0])*m3_process[0]

            # fourth_term_3[0] = j3_process[0]*dftech_dx1_process[0]*(v3_process[0] - va_process[0])*m1_process[0]+\
            #                 j3_process[0]*dftech_dx2_process[0]*(v3_process[0] - va_process[0])*m2_process[0]+\
            #                 j3_process[0]*dftech_dx3_process[0]*(v3_process[0] - va_process[0])*m3_process[0]

            discount_factor_temps[0] = -j1_process[0]-j2_process[0]*np.mean(damage_process[:,0],axis=0)-j3_process[0]*tech_process[0]-j4_process[0]
            discount_factor_temps_nodelta[0] = -j1_process[0]-j2_process[0]*np.mean(damage_process[:,0],axis=0)-j3_process[0]*tech_process[0]-j4_process[0]


            discount_factors[0] = -j1_process[0]-j2_process[0]*np.mean(damage_process[:,0],axis=0)-j3_process[0]*tech_process[0]-j4_process[0]


            discount_factor_nodelta_DisSep_Damage[0] = 0
            discount_factor_nodelta_DisSep_Tech[0] = 0
            discount_factor_nodelta[0] = 0


            discount_factor_nodeltadt_DisSep_Damage[0] = -np.exp(discount_factor_nodelta[0]) * (-j2_process[0]*np.mean(damage_process[:,0],axis=0)) 
            discount_factor_nodeltadt_DisSep_Tech[0] = -np.exp(discount_factor_nodelta[0]) * (-j3_process[0]*tech_process[0]) 
            discount_factor_nodeltadt[0] = - np.exp(discount_factor_nodelta[0]) * (-j2_process[0]*np.mean(damage_process[:,0],axis=0)-j3_process[0]*tech_process[0]) 



            undiscount_process[0] = delta*first_term[0]+second_term[0]+third_term_1[0]+third_term_2[0]+third_term_3[0]+third_term_4[0]+fourth_term_entropy[0]
            # undiscount_o_process[0] = delta*first_term[0]+second_term[0]+third_term_1[0]+third_term_2[0]+third_term_3[0]
            # undiscount_pg_process[0] = delta*first_term[0]+second_term[0]+third_term_1[0]+third_term_2[0]+third_term_3[0]+fourth_term_entropy[0]+fourth_term_2[0] +fourth_term_3[0]
            # undiscount_pg_f2_process[0] = delta*first_term_f2[0]+second_term[0]+third_term_1[0]+third_term_2[0]+third_term_3[0]+fourth_term_entropy[0]+fourth_term_2[0] +fourth_term_3[0]
            # undiscount_f_process[0] = delta*first_term_f[0]+second_term[0]+third_term_1[0]+third_term_2[0]+third_term_3[0]
            # undiscount_f_pg_process[0] = delta*first_term_f[0]+second_term[0]+third_term_1[0]+third_term_2[0]+third_term_3[0]+fourth_term_2[0] +fourth_term_3[0]
            discount_process[0] = undiscount_process[0] * np.exp(discount_factors[0])
            # discount_o_process[0] = undiscount_o_process[0] * np.exp(discount_factors[0])
            # discount_pg_process[0] = undiscount_pg_process[0] * np.exp(discount_factors[0])
            # discount_pg_f2_process[0] = undiscount_pg_f2_process[0] * np.exp(discount_factors[0])
            # discount_f_process[0] = undiscount_f_process[0] * np.exp(discount_factors[0])
            # discount_f_pg_process[0] = undiscount_f_pg_process[0] * np.exp(discount_factors[0])


            # first_term2[0] = uF1_process[0]*m1_process[0]+uF2_process[0]*m2_process[0]+uF3_process[0]*m3_process[0]
            # second_term2[0] = 0
            # third_term_12[0] = 0
            # third_term_22[0] =0
            # third_term_32[0] = 0
            # discount_factor_temps2[0] = 0
            # discount_factors2[0] = 0
            # undiscount_process2[0] = delta*first_term2[0]+second_term2[0]+third_term_12[0]+third_term_22[0]+third_term_32[0]
            # discount_process2[0] = undiscount_process2[0] * np.exp(discount_factors2[0])





            for t in range(N-1):
                # print(w_process[t], z_process[t], v_process[t])
                # print(w_process[t],z_process[t],v_process[t],x_process[t])
                muW_t = df_muW_interpolated([w_process[t],z_process[t],v_process[t],x_process[t]])
                sigmaW0_t = df_sigmaW0_interpolated([w_process[t],z_process[t],v_process[t],x_process[t]])
                
                muZ_t = df_muZ_interpolated([w_process[t],z_process[t],v_process[t],x_process[t]])
                sigmaZ0_t = df_sigmaZ0_interpolated([w_process[t],z_process[t],v_process[t],x_process[t]])
                
                muV_t = df_muV_interpolated([w_process[t],z_process[t],v_process[t],x_process[t]])
                sigmaV0_t = df_sigmaV0_interpolated([w_process[t],z_process[t],v_process[t],x_process[t]])

                muX_t = df_muX_interpolated([w_process[t],z_process[t],v_process[t],x_process[t]])
                sigmaX0_t = df_sigmaX0_interpolated([w_process[t],z_process[t],v_process[t],x_process[t]])

                mu_m1_t = m1_process[t]*df_muW_interpolated_m1([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m2_process[t]*df_muW_interpolated_m2([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m3_process[t]*df_muW_interpolated_m3([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m4_process[t]*df_muW_interpolated_m4([w_process[t],z_process[t],v_process[t],x_process[t]])

                sigma_m1_t = m1_process[t]*df_sigmaW0_interpolated_m1([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m2_process[t]*df_sigmaW0_interpolated_m2([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m3_process[t]*df_sigmaW0_interpolated_m3([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m4_process[t]*df_sigmaW0_interpolated_m4([w_process[t],z_process[t],v_process[t],x_process[t]])
                            
                mu_m2_t = m1_process[t]*df_muZ_interpolated_m1([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m2_process[t]*df_muZ_interpolated_m2([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m3_process[t]*df_muZ_interpolated_m3([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m4_process[t]*df_muZ_interpolated_m4([w_process[t],z_process[t],v_process[t],x_process[t]])
                            
                sigma_m2_t = m1_process[t]*df_sigmaZ0_interpolated_m1([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m2_process[t]*df_sigmaZ0_interpolated_m2([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m3_process[t]*df_sigmaZ0_interpolated_m3([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m4_process[t]*df_sigmaZ0_interpolated_m4([w_process[t],z_process[t],v_process[t],x_process[t]])
                            
                mu_m3_t = m1_process[t]*df_muV_interpolated_m1([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m2_process[t]*df_muV_interpolated_m2([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m3_process[t]*df_muV_interpolated_m3([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m4_process[t]*df_muV_interpolated_m4([w_process[t],z_process[t],v_process[t],x_process[t]])

                sigma_m3_t = m1_process[t]*df_sigmaV0_interpolated_m1([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m2_process[t]*df_sigmaV0_interpolated_m2([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m3_process[t]*df_sigmaV0_interpolated_m3([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m4_process[t]*df_sigmaV0_interpolated_m4([w_process[t],z_process[t],v_process[t],x_process[t]])
                

                mu_m4_t = m1_process[t]*df_muX_interpolated_m1([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m2_process[t]*df_muX_interpolated_m2([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m3_process[t]*df_muX_interpolated_m3([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m4_process[t]*df_muX_interpolated_m4([w_process[t],z_process[t],v_process[t],x_process[t]])

                sigma_m4_t = m1_process[t]*df_sigmaX0_interpolated_m1([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m2_process[t]*df_sigmaX0_interpolated_m2([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m3_process[t]*df_sigmaX0_interpolated_m3([w_process[t],z_process[t],v_process[t],x_process[t]])+\
                            m4_process[t]*df_sigmaX0_interpolated_m4([w_process[t],z_process[t],v_process[t],x_process[t]])
                

                m1_process[t+1] = m1_process[t] + mu_m1_t*dt + W1[t]*sigma_m1_t
                m2_process[t+1] = m2_process[t] + mu_m2_t*dt + W2[t]*sigma_m2_t
                m3_process[t+1] = m3_process[t] + mu_m3_t*dt + W3[t]*sigma_m3_t
                m4_process[t+1] = m4_process[t] + mu_m4_t*dt + W4[t]*sigma_m4_t
                        
                w_process[t+1] = w_process[t] + muW_t*dt + W1[t]*sigmaW0_t 
                z_process[t+1] = z_process[t] + muZ_t*dt + W2[t]*sigmaZ0_t
                v_process[t+1] = v_process[t] + muV_t*dt + W3[t]*sigmaV0_t
                x_process[t+1] = x_process[t] + muX_t*dt + W4[t]*sigmaX0_t

                u1_process[t+1] = df_U_interpolated_m1([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                u2_process[t+1] = df_U_interpolated_m2([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                u3_process[t+1] = df_U_interpolated_m3([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                u4_process[t+1] = df_U_interpolated_m4([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])


                # u11_process[t+1] = df_U1_interpolated_m1([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                # u12_process[t+1] = df_U1_interpolated_m2([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                # u13_process[t+1] = df_U1_interpolated_m3([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])


                # u21_process[t+1] = df_U2_interpolated_m1([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                # u22_process[t+1] = df_U2_interpolated_m2([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                # u23_process[t+1] = df_U2_interpolated_m3([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])

                # u31_process[t+1] = df_U3_interpolated_m1([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                # u32_process[t+1] = df_U3_interpolated_m2([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                # u33_process[t+1] = df_U3_interpolated_m3([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])


                # u41_process[t+1] = df_U4_interpolated_m1([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                # u42_process[t+1] = df_U4_interpolated_m2([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                # u43_process[t+1] = df_U4_interpolated_m3([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])



                # uF1_process[t+1] = df_UF_interpolated_m1([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                # uF2_process[t+1] = df_UF_interpolated_m2([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                # uF3_process[t+1] = df_UF_interpolated_m3([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])

                # uF21_process[t+1] = df_UF2_interpolated_m1([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                # uF22_process[t+1] = df_UF2_interpolated_m2([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                # uF23_process[t+1] = df_UF2_interpolated_m3([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])

                j1_process[t+1] = df_J1_interpolated([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                j2_process[t+1] = df_J2_interpolated([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                j3_process[t+1] = df_J3_interpolated([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                j4_process[t+1] = df_J4_interpolated([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])


                v1_process[t+1] = df_V1_interpolated([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])



                for i in range(n_damage):

                    v2_process[i,t+1] = df_V2_interpolated[i]([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                    damage_process[i,t+1] = f_damage_interpolated[i]([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])


                tech_process[t+1] = f_tech_interpolated([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])

                Entropy_x1_process[t+1] = df_Entropy_interpolated_m1([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                Entropy_x2_process[t+1] = df_Entropy_interpolated_m2([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                Entropy_x3_process[t+1] = df_Entropy_interpolated_m3([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                Entropy_x4_process[t+1] = df_Entropy_interpolated_m4([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])

                # v2_process[t+1] = df_V2_interpolated([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])

                v3_process[t+1] = df_V3_interpolated([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                v4_process[t+1] = df_V4_interpolated([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])

                va_process[t+1] = df_V_interpolated([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])

                dj1_dx1_process[t+1] = df_J1_interpolated_m1([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dj2_dx2_process[t+1] = df_J2_interpolated_m2([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dj3_dx3_process[t+1] = df_J3_interpolated_m3([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dj4_dx4_process[t+1] = df_J4_interpolated_m4([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                
                dv1_dx1_process[t+1] = df_V1_interpolated_m1([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dv1_dx2_process[t+1] = df_V1_interpolated_m2([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dv1_dx3_process[t+1] = df_V1_interpolated_m3([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dv1_dx4_process[t+1] = df_V1_interpolated_m4([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])


                for i in range(n_damage):
                    dv2_dx1_process[i,t+1] = df_V2_interpolated_m1[i]([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                    dv2_dx2_process[i,t+1] = df_V2_interpolated_m2[i]([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                    dv2_dx3_process[i,t+1] = df_V2_interpolated_m3[i]([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                    dv2_dx4_process[i,t+1] = df_V2_interpolated_m4[i]([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])


                    dfdamage_dx1_process[i,t+1] = df_f_damage_interpolated_m1[i]([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                    dfdamage_dx2_process[i,t+1] = df_f_damage_interpolated_m2[i]([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                    dfdamage_dx3_process[i,t+1] = df_f_damage_interpolated_m3[i]([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                    dfdamage_dx4_process[i,t+1] = df_f_damage_interpolated_m4[i]([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])

                dftech_dx1_process[t+1] = df_f_tech_m1_interpolated([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dftech_dx2_process[t+1] = df_f_tech_m2_interpolated([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dftech_dx3_process[t+1] = df_f_tech_m3_interpolated([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dftech_dx4_process[t+1] = df_f_tech_m4_interpolated([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])

                # dv2_dx1_process[t+1] = df_V2_interpolated_m1([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                # dv2_dx2_process[t+1] = df_V2_interpolated_m2([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                # dv2_dx3_process[t+1] = df_V2_interpolated_m3([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])




                dv3_dx1_process[t+1] = df_V3_interpolated_m1([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dv3_dx2_process[t+1] = df_V3_interpolated_m2([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dv3_dx3_process[t+1] = df_V3_interpolated_m3([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dv3_dx4_process[t+1] = df_V3_interpolated_m4([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                
                dv4_dx1_process[t+1] = df_V4_interpolated_m1([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dv4_dx2_process[t+1] = df_V4_interpolated_m2([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dv4_dx3_process[t+1] = df_V4_interpolated_m3([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dv4_dx4_process[t+1] = df_V4_interpolated_m4([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                
                dv_dx1_process[t+1] = df_V_interpolated_m1([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dv_dx2_process[t+1] = df_V_interpolated_m2([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dv_dx3_process[t+1] = df_V_interpolated_m3([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])
                dv_dx4_process[t+1] = df_V_interpolated_m4([w_process[t+1],z_process[t+1],v_process[t+1],x_process[t+1]])



                first_term[t+1] = u1_process[t+1]*m1_process[t+1]+\
                                u2_process[t+1]*m2_process[t+1]+\
                                u3_process[t+1]*m3_process[t+1]+\
                                u4_process[t+1]*m4_process[t+1]
                # first_term_f[t+1] = uF1_process[t+1]*m1_process[t+1]+uF2_process[t+1]*m2_process[t+1]+uF3_process[t+1]*m3_process[t+1]
                # first_term_f2[t+1] = uF21_process[t+1]*m1_process[t+1]+uF22_process[t+1]*m2_process[t+1]+uF23_process[t+1]*m3_process[t+1]

                # second_term[t+1] = dj1_dx1_process[t+1]*(v1_process[t+1] - va_process[t+1])*m1_process[t+1]+\
                #                 dj2_dx2_process[t+1]*(v2_process[t+1] - va_process[t+1])*m2_process[t+1]+\
                #                 dj3_dx3_process[t+1]*(v3_process[t+1] - va_process[t+1])*m3_process[t+1]
                second_term[t+1] = dj1_dx1_process[t+1]*(v1_process[t+1] - va_process[t+1])*m1_process[t+1]+\
                            dj2_dx2_process[t+1]*np.mean(damage_process[:,t+1]*(v2_process[:,t+1] - va_process[t+1]),axis=0)*m2_process[t+1]+\
                            dj3_dx3_process[t+1]*tech_process[t+1]*(v3_process[t+1] - va_process[t+1])*m3_process[t+1]+\
                            dj4_dx4_process[t+1]*(v4_process[t+1] - va_process[t+1])*m4_process[t+1]
                
                second_term_1[t+1] = dj1_dx1_process[t+1]*(v1_process[t+1] - va_process[t+1])*m1_process[t+1]
                second_term_2[t+1] = dj2_dx2_process[t+1]*np.mean(damage_process[:,t+1]*(v2_process[:,t+1] - va_process[t+1]),axis=0)*m2_process[t+1]
                second_term_3[t+1] = dj3_dx3_process[t+1]*tech_process[t+1]*(v3_process[t+1] - va_process[t+1])*m3_process[t+1]
                second_term_4[t+1] = dj4_dx4_process[t+1]*(v4_process[t+1] - va_process[t+1])*m4_process[t+1]


                third_term_1[t+1] = j1_process[t+1]*(dv1_dx1_process[t+1])*m1_process[t+1]+\
                                j1_process[t+1]*(dv1_dx2_process[t+1])*m2_process[t+1]+\
                                j1_process[t+1]*(dv1_dx3_process[t+1])*m3_process[t+1]+\
                                j1_process[t+1]*(dv1_dx4_process[t+1])*m4_process[t+1]
                # third_term_2[t+1] = j2_process[t+1]*(dv2_dx1_process[t+1])*m1_process[t+1]+\
                #                 j2_process[t+1]*(dv2_dx2_process[t+1])*m2_process[t+1]+\
                #                 j2_process[t+1]*(dv2_dx3_process[t+1])*m3_process[t+1]
                third_term_2[t+1] = j2_process[t+1]*np.mean(damage_process[:,t+1]*dv2_dx1_process[:,t+1],axis=0)*m1_process[t+1]+\
                                j2_process[t+1]*np.mean(damage_process[:,t+1]*dv2_dx2_process[:,t+1],axis=0)*m2_process[t+1]+\
                                j2_process[t+1]*np.mean(damage_process[:,t+1]*dv2_dx3_process[:,t+1],axis=0)*m3_process[t+1]+\
                                j2_process[t+1]*np.mean(damage_process[:,t+1]*dv2_dx4_process[:,t+1],axis=0)*m4_process[t+1]
                # third_term_3[t+1] = j3_process[t+1]*(dv3_dx1_process[t+1])*m1_process[t+1]+\
                #                 j3_process[t+1]*(dv3_dx2_process[t+1])*m2_process[t+1]+\
                #                 j3_process[t+1]*(dv3_dx3_process[t+1])*m3_process[t+1]
                third_term_3[t+1] = j3_process[t+1]*tech_process[t+1]*(dv3_dx1_process[t+1])*m1_process[t+1]+\
                                j3_process[t+1]*tech_process[t+1]*(dv3_dx2_process[t+1])*m2_process[t+1]+\
                                j3_process[t+1]*tech_process[t+1]*(dv3_dx3_process[t+1])*m3_process[t+1]+\
                                j3_process[t+1]*tech_process[t+1]*(dv3_dx4_process[t+1])*m4_process[t+1]

                third_term_4[t+1] = j4_process[t+1]*(dv4_dx1_process[t+1])*m1_process[t+1]+\
                                j4_process[t+1]*(dv4_dx2_process[t+1])*m2_process[t+1]+\
                                j4_process[t+1]*(dv4_dx3_process[t+1])*m3_process[t+1]+\
                                j4_process[t+1]*(dv4_dx4_process[t+1])*m4_process[t+1]
                

                # fourth_term_2[t+1] = j2_process[t+1]*np.mean(dfdamage_dx1_process[:,t+1]*(v2_process[:,t+1] - va_process[t+1]),axis=0)*m1_process[t+1]+\
                #                 j2_process[t+1]*np.mean(dfdamage_dx2_process[:,t+1]*(v2_process[:,t+1] - va_process[t+1]),axis=0)*m2_process[t+1]+\
                #                 j2_process[t+1]*np.mean(dfdamage_dx3_process[:,t+1]*(v2_process[:,t+1] - va_process[t+1]),axis=0)*m3_process[t+1]
                # fourth_term_3[t+1] = j3_process[t+1]*dftech_dx1_process[t+1]*(v3_process[t+1] - va_process[t+1])*m1_process[t+1]+\
                #                 j3_process[t+1]*dftech_dx2_process[t+1]*(v3_process[t+1] - va_process[t+1])*m2_process[t+1]+\
                #                 j3_process[t+1]*dftech_dx3_process[t+1]*(v3_process[t+1] - va_process[t+1])*m3_process[t+1]

                fourth_term_entropy[t+1] = Entropy_x1_process[t+1]*m1_process[t+1]+\
                                            Entropy_x2_process[t+1]*m2_process[t+1]+\
                                            Entropy_x3_process[t+1]*m3_process[t+1]+\
                                            Entropy_x4_process[t+1]*m4_process[t+1]

                # fourth_term_2[t+1] = j2_process[t+1]*np.mean(dfdamage_dx1_process[:,t+1]*(v2_process[:,t+1] - va_process[t+1]),axis=0)*m1_process[t+1]+\
                #                 j2_process[t+1]*np.mean(dfdamage_dx2_process[:,t+1]*(v2_process[:,t+1] - va_process[t+1]),axis=0)*m2_process[t+1]+\
                #                 j2_process[t+1]*np.mean(dfdamage_dx3_process[:,t+1]*(v2_process[:,t+1] - va_process[t+1]),axis=0)*m3_process[t+1]
                # fourth_term_3[0] = j3_process[0]*dftech_dx1_process[0]*(v3_process[0] - va_process[0])*m1_process[0]+\
                #                 j3_process[0]*dftech_dx2_process[0]*(v3_process[0] - va_process[0])*m2_process[0]+\
                #                 j3_process[0]*dftech_dx3_process[0]*(v3_process[0] - va_process[0])*m3_process[0]

                # fourth_term_3[t+1] = j3_process[t+1]*dftech_dx1_process[t+1]*(v3_process[t+1] - va_process[t+1])*m1_process[t+1]+\
                #                 j3_process[t+1]*dftech_dx2_process[t+1]*(v3_process[t+1] - va_process[t+1])*m2_process[t+1]+\
                #                 j3_process[t+1]*dftech_dx3_process[t+1]*(v3_process[t+1] - va_process[t+1])*m3_process[t+1]

                discount_factor_temps[t+1] =  -delta-j1_process[t+1]-j2_process[t+1]*np.mean(damage_process[:,t+1],axis=0)-j3_process[t+1]*tech_process[t+1]-j4_process[t+1]
                discount_factor_temps_nodelta[t+1] =  -j1_process[t+1]-j2_process[t+1]*np.mean(damage_process[:,t+1],axis=0)-j3_process[t+1]*tech_process[t+1]-j4_process[t+1]

                discount_factors[t+1] = discount_factors[t]+discount_factor_temps[t+1]*dt


                discount_factor_nodelta_DisSep_Damage[t+1] = discount_factor_nodelta_DisSep_Damage[t]+(-j2_process[t+1]*np.mean(damage_process[:,t+1],axis=0))*dt
                discount_factor_nodelta_DisSep_Tech[t+1] = discount_factor_nodelta_DisSep_Tech[t]+(-j3_process[t+1]*tech_process[t+1])*dt
                discount_factor_nodelta[t+1] = discount_factor_nodelta[t]+discount_factor_temps_nodelta[t+1]*dt



                
                discount_factor_nodeltadt_DisSep_Damage[t+1] = -np.exp(discount_factor_nodelta[t+1]) * ( - j2_process[t+1]*np.mean(damage_process[:,t+1],axis=0)) 
                discount_factor_nodeltadt_DisSep_Tech[t+1] = -np.exp(discount_factor_nodelta[t+1]) * (  - j3_process[t+1]*tech_process[t+1]) 
                discount_factor_nodeltadt[t+1] = -np.exp(discount_factor_nodelta[t+1]) * discount_factor_temps_nodelta[t+1]


                undiscount_process[t+1] = delta*first_term[t+1]+second_term[t+1]+third_term_1[t+1]+third_term_2[t+1]+third_term_3[t+1]+third_term_4[t+1]+fourth_term_entropy[t+1]
                # undiscount_pg_process[t+1] = delta*first_term[t+1]+second_term[t+1]+third_term_1[t+1]+third_term_2[t+1]+third_term_3[t+1]+fourth_term_entropy[t+1]+fourth_term_2[t+1] +fourth_term_3[t+1]
                # undiscount_pg_f2_process[t+1] = delta*first_term_f2[t+1]+second_term[t+1]+third_term_1[t+1]+third_term_2[t+1]+third_term_3[t+1]+fourth_term_entropy[t+1]+fourth_term_2[t+1] +fourth_term_3[t+1]
                # undiscount_f_process[t+1] = delta*first_term_f[t+1]+second_term[t+1]+third_term_1[t+1]+third_term_2[t+1]+third_term_3[t+1]
                # undiscount_f_pg_process[t+1] = delta*first_term_f[t+1]+second_term[t+1]+third_term_1[t+1]+third_term_2[t+1]+third_term_3[t+1]+fourth_term_2[t+1] +fourth_term_3[t+1]
                discount_process[t+1] = undiscount_process[t+1] * np.exp(discount_factors[t+1])
                # discount_pg_process[t+1] = undiscount_pg_process[t+1] * np.exp(discount_factors[t+1])
                # discount_pg_f2_process[t+1] = undiscount_pg_f2_process[t+1] * np.exp(discount_factors[t+1])
                # discount_f_process[t+1] = undiscount_f_process[t+1] * np.exp(discount_factors[t+1])
                # discount_f_pg_process[t+1] = undiscount_f_pg_process[t+1] * np.exp(discount_factors[t+1])

                # first_term2[t+1] = uF1_process[t+1]*m1_process[t+1]+uF2_process[t+1]*m2_process[t+1]+uF3_process[t+1]*m3_process[t+1]
                # second_term2[t+1] = 0
                # third_term_12[t+1] = 0
                # third_term_22[t+1] =0
                # third_term_32[t+1] = 0
                # discount_factor_temps2[t+1] = -delta
                # discount_factors2[t+1] = discount_factors2[t]+discount_factor_temps2[t+1]
                # undiscount_process2[t+1] = delta*first_term2[t+1]+second_term2[t+1]+third_term_12[t+1]+third_term_22[t+1]+third_term_32[t+1]
                # discount_process2[t+1] = undiscount_process2[t+1] * np.exp(discount_factors2[t+1])
                time_process[t+1] = time_process[t]+dt




            w_processes.append(w_process)
            z_processes.append(z_process)
            v_processes.append(v_process)
            x_processes.append(x_process)
            m1_processes.append(m1_process)
            m2_processes.append(m2_process)
            m3_processes.append(m3_process)
            m4_processes.append(m4_process)
            u1_processes.append(u1_process)
            u2_processes.append(u2_process)
            u3_processes.append(u3_process)
            u4_processes.append(u4_process)
            
            # u11_processes.append(u11_process)
            # u12_processes.append(u12_process)
            # u13_processes.append(u13_process)
            
            # u21_processes.append(u21_process)
            # u22_processes.append(u22_process)
            # u23_processes.append(u23_process)


            # u31_processes.append(u31_process)
            # u32_processes.append(u32_process)
            # u33_processes.append(u33_process)
            
            # u41_processes.append(u41_process)
            # u42_processes.append(u42_process)
            # u43_processes.append(u43_process)


            # uF1_processes.append(uF1_process)
            # uF2_processes.append(uF2_process)
            # uF3_processes.append(uF3_process)

            damage_processes.append(damage_process)
            tech_processes.append(tech_process)

            j1_processes.append(j1_process)
            j2_processes.append(j2_process)
            j3_processes.append(j3_process)
            j4_processes.append(j4_process)
            v1_processes.append(v1_process)
            v2_processes.append(v2_process)
            v3_processes.append(v3_process)
            v4_processes.append(v4_process)
            va_processes.append(va_process)
            dj1_dx1_processes.append(dj1_dx1_process)
            dj2_dx2_processes.append(dj2_dx2_process)
            dj3_dx3_processes.append(dj3_dx3_process)
            dj4_dx4_processes.append(dj4_dx4_process)
            dv1_dx1_processes.append(dv1_dx1_process)
            dv1_dx2_processes.append(dv1_dx2_process)
            dv1_dx3_processes.append(dv1_dx3_process)
            dv1_dx4_processes.append(dv1_dx4_process)
            dv2_dx1_processes.append(dv2_dx1_process)
            dv2_dx2_processes.append(dv2_dx2_process)
            dv2_dx3_processes.append(dv2_dx3_process)
            dv2_dx4_processes.append(dv2_dx4_process)
            dv3_dx1_processes.append(dv3_dx1_process)
            dv3_dx2_processes.append(dv3_dx2_process)
            dv3_dx3_processes.append(dv3_dx3_process)
            dv3_dx4_processes.append(dv3_dx4_process)
            dv4_dx1_processes.append(dv4_dx1_process)
            dv4_dx2_processes.append(dv4_dx2_process)
            dv4_dx3_processes.append(dv4_dx3_process)
            dv4_dx4_processes.append(dv4_dx4_process)
            dv_dx1_processes.append(dv_dx1_process)
            dv_dx2_processes.append(dv_dx2_process)
            dv_dx3_processes.append(dv_dx3_process)
            dv_dx4_processes.append(dv_dx4_process)
            first_terms.append(first_term)
            # first_term_fs.append(first_term_f)
            second_terms.append(second_term)
            second_term_1s.append(second_term_1)
            second_term_2s.append(second_term_2)
            second_term_3s.append(second_term_3)
            second_term_4s.append(second_term_4)
            third_term_1s.append(third_term_1)
            third_term_2s.append(third_term_2)
            third_term_3s.append(third_term_3)
            third_term_4s.append(third_term_4)
            fourth_term_entropys.append(fourth_term_entropy)
            # fourth_term_2s.append(fourth_term_2)
            # fourth_term_3s.append(fourth_term_3)
            # fourth_term_4s.append(fourth_term_4)
            undiscount_processes.append(undiscount_process)
            # undiscount_o_processes.append(undiscount_o_process)
            # undiscount_pg_processes.append(undiscount_pg_process)
            # undiscount_pg_f2_processes.append(undiscount_pg_f2_process)
            # undiscount_f_processes.append(undiscount_f_process)
            # undiscount_f_pg_processes.append(undiscount_f_pg_process)
            discount_factor_temp_processes.append(discount_factor_temps)
            discount_factor_processes.append(discount_factors)

            discount_factor_nodelta_DisSep_Damage_processes.append(discount_factor_nodelta_DisSep_Damage)
            discount_factor_nodelta_DisSep_Tech_processes.append(discount_factor_nodelta_DisSep_Tech)
            discount_factor_nodelta_processes.append(discount_factor_nodelta)

            discount_factor_nodeltadt_DisSep_Damage_processes.append(discount_factor_nodeltadt_DisSep_Damage)
            discount_factor_nodeltadt_DisSep_Tech_processes.append(discount_factor_nodeltadt_DisSep_Tech)
            discount_factor_nodeltadt_processes.append(discount_factor_nodeltadt)

            discount_processes.append(discount_process)
            # discount_o_processes.append(discount_o_process)
            # discount_pg_processes.append(discount_pg_process)
            # discount_pg_f2_processes.append(discount_pg_f2_process)
            # discount_f_processes.append(discount_f_process)
            # discount_f_pg_processes.append(discount_f_pg_process)

            # first_terms2.append(first_term2)
            # second_terms2.append(second_term2)
            # third_term_1s2.append(third_term_12)
            # third_term_2s2.append(third_term_22)
            # third_term_3s2.append(third_term_32)
            # undiscount_processes2.append(undiscount_process2)
            # discount_factor_temp_processes2.append(discount_factor_temps2)
            # discount_factor_processes2.append(discount_factors2)
            # discount_processes2.append(discount_process2)
            time_processes.append(time_process)

            # print(time_process)
            # print(f"integral = {np.trapz(discount_factor_nodeltadt, time_process)}")
            # Integral = np.trapz(np.mean(discount_factor_nodeltadt_processes, axis=0), np.mean(time_process,axis=0))
        except:
            print('Error at seed', i)
            continue
        np.savez(Data_Dir + File_Dir+'Sim_'+str(seed)+f'_{IntPeriod}_{str(dt)}_{m0}_ID.npz', W_process=w_processes, Z_process=z_processes, V_process=v_processes, X_process=x_processes, m1_process=m1_processes, m2_process=m2_processes, m3_process=m3_processes, m4_process=m4_processes, u1_process=u1_processes, u2_process=u2_processes, u3_process=u3_processes, u4_process=u4_processes,\
                # u11_process=u11_processes, u12_process=u12_processes, u13_process=u13_processes, u21_process=u21_processes, u22_process=u22_processes, u23_process=u23_processes, \
                # u31_process=u31_processes, u32_process=u32_processes, u33_process=u33_processes, u41_process=u41_processes, u42_process=u42_processes, u43_process=u43_processes, uF1_process=uF1_processes, uF2_process=uF2_processes, uF3_process=uF3_processes, \
                time_process= time_process, j1_process=j1_processes, j2_process=j2_processes, j3_process=j3_processes, j4_process=j4_processes, v1_process=v1_processes, v2_process=v2_processes, v3_process=v3_processes, v4_process=v4_processes, va_process=va_processes,\
                dj1_dx1_process=dj1_dx1_processes, dj2_dx2_process=dj2_dx2_processes, dj3_dx3_process=dj3_dx3_processes,dj4_dx4_process=dj4_dx4_processes,\
                dv1_dx1_process=dv1_dx1_processes, dv1_dx2_process=dv1_dx2_processes, dv1_dx3_process=dv1_dx3_processes,dv1_dx4_process=dv1_dx4_processes,\
                dv2_dx1_process=dv2_dx1_processes, dv2_dx2_process=dv2_dx2_processes, dv2_dx3_process=dv2_dx3_processes,dv2_dx4_process=dv2_dx4_processes,\
                dv3_dx1_process=dv3_dx1_processes, dv3_dx2_process=dv3_dx2_processes, dv3_dx3_process=dv3_dx3_processes,dv3_dx4_process=dv3_dx4_processes,\
                dv4_dx1_process=dv4_dx1_processes, dv4_dx2_process=dv4_dx2_processes, dv4_dx3_process=dv4_dx3_processes,dv4_dx4_process=dv4_dx4_processes,\
                dv_dx1_process=dv_dx1_processes, dv_dx2_process=dv_dx2_processes, dv_dx3_process=dv_dx3_processes, dv_dx4_process=dv_dx4_processes,\
                damage_processes = damage_processes, tech_processes=tech_processes,\
                # first_term=first_terms, first_term_f=first_term_fs, second_term=second_terms, third_term_1=third_term_1s, third_term_2=third_term_2s, third_term_3=third_term_3s, undiscount_process=undiscount_processes, undiscount_pg_process=undiscount_pg_processes, undiscount_pg_f2_process=undiscount_pg_f2_processes, undiscount_f_process=undiscount_f_processes, undiscount_f_pg_process=undiscount_f_pg_processes, discount_factor_temp=discount_factor_temp_processes, discount_factor=discount_factor_processes, discount_process=discount_processes, discount_pg_process=discount_pg_processes, discount_pg_f2_process=discount_pg_f2_processes, discount_f_process=discount_f_processes, discount_f_pg_process=discount_f_pg_processes)
                first_term=first_terms, second_term=second_terms, third_term_1=third_term_1s, third_term_2=third_term_2s, third_term_3=third_term_3s, third_term_4=third_term_4s, fourth_term_entropy=fourth_term_entropys, undiscount_process=undiscount_processes, \
                discount_factor_temp=discount_factor_temp_processes, discount_factor=discount_factor_processes, \
                discount_factor_nodelta=discount_factor_nodelta_processes, discount_factor_nodelta_DisSep_Damage=discount_factor_nodelta_DisSep_Damage_processes, discount_factor_nodelta_DisSep_Tech=discount_factor_nodelta_DisSep_Tech_processes, \
                discount_factor_nodeltadt=discount_factor_nodeltadt_processes, discount_factor_nodeltadt_DisSep_Damage = discount_factor_nodeltadt_DisSep_Damage_processes, discount_factor_nodeltadt_DisSep_Tech= discount_factor_nodeltadt_DisSep_Tech_processes, \
                discount_process=discount_processes, \
                second_term_1=second_term_1s,second_term_2=second_term_2s,second_term_3=second_term_3s,second_term_4=second_term_4s)
                # first_term2=first_terms2, second_term2=second_terms2, third_term_12=third_term_1s2, third_term_22=third_term_2s2, third_term_32=third_term_3s2, undiscount_process2=undiscount_processes2, discount_factor_temp2=discount_factor_temp_processes2, discount_factor2=discount_factor_processes2, discount_process2=discount_processes2)
        print('Data saved')

    return 0


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
        with open(Data_Dir+ File_Dir2 + "model_tech1_pre_damage"+"_prepareID", "rb") as f:
            data = pickle.load(f)
        # with open(Data_Dir+ File_Dir + "model_tech2_pre_damage"+"_UD_prepare", "rb") as f:
        #     model_tech2_pre_damage_ME_base = pickle.load(f)
        # 
        # n_bar = n_bar1
        

    simulate_pre(data, Data_Dir, File_Dir2)

    # with open(Data_Dir + File_Dir+"model_tech1_pre_damage"+"_UD_simul_{}".format(IntPeriod)+ scheme + "_" +HJB_solution, "wb") as f:
    #     pickle.dump(res,f)


    
    return 0


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):
                for id_rho in range(len(rhoarr)):

                    res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],xia2arr[id_xiag],xik2arr[id_xiag],xic2arr[id_xiag],xij2arr[id_xiag],xid2arr[id_xiag],xig2arr[id_xiag],rhoarr[id_rho],psi0arr[id_psi0],psi1arr[id_psi1],varrhoarr[id_varrho])


