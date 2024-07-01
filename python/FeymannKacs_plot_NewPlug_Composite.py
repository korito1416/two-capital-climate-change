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
import plotly.graph_objects as go
import plotly.offline as pyo
import matplotlib.pyplot as plt
from src.Utility import finite_1D

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

parser.add_argument("--m0",type=str)

args = parser.parse_args()



# Update = args.Update
IntPeriod = args.IntPeriod
dt = 1/12
# dt = 1

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


# print("bY_short={:d}".format(nY_short))
(K_mat, Y_mat, L_mat) = np.meshgrid(K, Y_short, L, indexing="ij")




# # if len(xicarr)==3 and min(xikarr)==0.075:
# labellist = ['More Aversion', 'Less Aversion', 'Neutrality']
# Filename = 'Aversion Intensity'
# # Filename = 'Aversion Intensity_old'
# # Filename = 'Aversion Intensity_onlyj'
# # Filename = 'Aversion Intensity_onlyk'
# colors = ['blue','red', 'green', 'cyan', 'purple']
# colors = ['blue','red', 'green', 'cyan', 'purple']
# print("define success")


# if min(xikarr)==0.150:
#     labellist = ['Less Aversion']
#     Filename = 'Aversion Intensity'
#     print("define success")


# if min(xikarr)>1:
#     labellist = ['Neutrality']
#     Filename = 'Aversion Intensity'
#     print("define success")

# if min(xij2arr)==0.150 and min(xijarr)>1:
#     labellist = ['Pre Less Aversion Post Neutrality']
#     Filename = 'Aversion Intensity'
#     print("define success")

# if min(xij2arr)==0.150 and min(xijarr)==0.150:
#     labellist = ['Pre Less Aversion Post Less Aversion']
#     Filename = 'Aversion Intensity'
#     print("define success")


# if min(xij2arr)>1 and min(xijarr)>1:
#     labellist = ['Pre Neutrality Post Neutrality']
#     Filename = 'Aversion Intensity'
#     print("define success")


# if min(xij2arr)>1 and min(xijarr)==0.150:
#     labellist = ['Pre Neutrality Post Less Aversion']
#     Filename = 'Aversion Intensity'
#     print("define success")

# labellist = ['pre neutrality post aversion', 'pre aversion post neutrality', 'pre aversion post aversion', 'pre neutrality post neutrality']
# Filename = 'Aversion Intensity'
# print("define success")


labellist = ['pre neutrality post aversion', 'pre aversion post neutrality', 'pre aversion post aversion', 'pre neutrality post neutrality']

if min(xij2arr)==0.075 and min(xik2arr)>1:

    Filename = 'Aversion Intensity_More'
elif min(xij2arr)==0.150 and min(xik2arr)>1:

    Filename = 'Aversion Intensity_Less'
elif min(xij2arr)==0.005 and min(xik2arr)>1:

    Filename = 'Aversion Intensity_Extreme'
    
elif min(xij2arr)==0.150 and min(xik2arr)==0.150 :

    Filename = 'All Channel On Less'

print("define success")

# Filename = 'Aversion Intensity_old'
# Filename = 'Aversion Intensity_onlyj'
# Filename = 'Aversion Intensity_onlyk'
colors = ['blue','red', 'green', 'cyan', 'purple']
colors = ['blue','red', 'green', 'cyan', 'purple']


plt.style.use('classic')
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["figure.figsize"] = (16,9)
plt.rcParams["figure.dpi"] = 1000
plt.rcParams["font.size"] = 18
plt.rcParams["legend.frameon"] = False
plt.rcParams["lines.linewidth"] = 10


print("After, figure default size is: ", plt.rcParams["savefig.bbox"])
print("After, figure default size is: ", plt.rcParams["figure.figsize"])
print("After, figure default dpi is: ", plt.rcParams["figure.dpi"])
print("After, figure default size is: ", plt.rcParams["font.size"])
print("After, legend.frameon is: ", plt.rcParams["legend.frameon"])
print("After, lines.linewidth is: ", plt.rcParams["lines.linewidth"])


os.makedirs("./figure/"+args.dataname+"/", exist_ok=True)

Plot_Dir = "./figure/"+args.dataname+"/"


def Damage_Intensity(Yt, y_bar_lower=1.5):
    r_1 = 1.5
    r_2 = 2.5
    Intensity = r_1 * (np.exp(r_2 / 2 * (Yt - y_bar_lower)**2) -1) * (Yt > y_bar_lower)
    return Intensity






def model_simulation_generate(id, xi_a,xi_k,xi_c,xi_j,xi_d,xi_g,xi_a2,xi_k2,xi_c2,xi_j2,xi_d2,xi_g2,rho,psi_0,psi_1,varrho):

    # Output_Dir = "/scratch/bincheng/"
    Output_Dir = args.outputname
    Data_Dir = Output_Dir+"abatement/data_2tech/"+args.dataname+"/"
    # File_Dir = "xi_a_{}_xi_k_{}_xi_c_{}_xi_j_{}_xi_d_{}_xi_g_{}_psi_0_{}_psi_1_{}_varrho_{}_rho_{}_" .format(xi_a,xi_k,xi_c,xi_j,xi_d,xi_g,psi_0,psi_1,varrho,rho)
    # File_Dir = "xi_a_{}_xi_k_{}_xi_c_{}_xi_j_{}_xi_d_{}_xi_g_{}_psi_0_{}_psi_1_{}_varrho_{}_rho_{}_delta_{}_" .format(xi_a, xi_k, xi_c, xi_j, xi_d, xi_g, psi_0,psi_1, varrho, rho, delta)
    File_Dir2 = "xi_k_{}_xi_c_{}_xi_j_{}_xi_d_{}_xi_g_{}_xi_k2_{}_xi_c2_{}_xi_j2_{}_xi_d2_{}_xi_g2_{}_psi_0_{}_psi_1_{}_varrho_{}_rho_{}_delta_{}_" .format(xi_k, xi_c, xi_j, xi_d, xi_g, xi_k2, xi_c2, xi_j2, xi_d2, xi_g2, psi_0,psi_1, varrho, rho, delta)

    # print(Data_Dir)
    N = int(IntPeriod/dt)
    W_processes = np.empty((0,N))
    Z_processes = np.empty((0,N))
    V_processes = np.empty((0,N))
    X_processes = np.empty((0,N))
    m1_processes = np.empty((0,N))
    m2_processes = np.empty((0,N))
    m3_processes = np.empty((0,N))
    m4_processes = np.empty((0,N))
    # u1_processes = np.empty((0,N))
    # u2_processes = np.empty((0,N))
    # u3_processes = np.empty((0,N))
    
    # u11_processes = np.empty((0,N))
    # u12_processes = np.empty((0,N))
    # u13_processes = np.empty((0,N))
    
    # u21_processes = np.empty((0,N))
    # u22_processes = np.empty((0,N))
    # u23_processes = np.empty((0,N))


    # u31_processes = np.empty((0,N))
    # u32_processes = np.empty((0,N))
    # u33_processes = np.empty((0,N))
    
    # u41_processes = np.empty((0,N))
    # u42_processes = np.empty((0,N))
    # u43_processes = np.empty((0,N))


    # uF1_processes = np.empty((0,N))
    # uF2_processes = np.empty((0,N))
    # uF3_processes = np.empty((0,N))

    # j1_processes = np.empty((0,N))
    # j2_processes = np.empty((0,N))
    # j3_processes = np.empty((0,N))
    # v1_processes = np.empty((0,N))
    # v2_processes = np.empty((0,N))
    # v3_processes = np.empty((0,N))
    # v_processes = np.empty((0,N))
    # dj1_dx1_processes = np.empty((0,N))
    # dj2_dx2_processes = np.empty((0,N))
    # dj3_dx3_processes = np.empty((0,N))
    # dv1_dx1_processes = np.empty((0,N))
    # dv1_dx2_processes = np.empty((0,N))
    # dv1_dx3_processes = np.empty((0,N))
    # dv2_dx1_processes = np.empty((0,N))
    # dv2_dx2_processes = np.empty((0,N))
    # dv2_dx3_processes = np.empty((0,N))
    # dv3_dx1_processes = np.empty((0,N))
    # dv3_dx2_processes = np.empty((0,N))
    # dv3_dx3_processes = np.empty((0,N))
    dv_dx1_processes = np.empty((0,N))
    dv_dx2_processes = np.empty((0,N))
    dv_dx3_processes = np.empty((0,N))
    dv_dx4_processes = np.empty((0,N))
    first_terms = np.empty((0,N))
    second_terms = np.empty((0,N))
    third_term_1s = np.empty((0,N))
    third_term_2s = np.empty((0,N))
    third_term_3s = np.empty((0,N))
    third_term_4s = np.empty((0,N))
    fourth_term_entropys = np.empty((0,N))
    undiscount_processes = np.empty((0,N))
    undiscount_o_processes = np.empty((0,N))
    undiscount_pg_processes = np.empty((0,N))
    undiscount_pg_f2_processes = np.empty((0,N))
    undiscount_f_processes = np.empty((0,N))
    undiscount_f_pg_processes = np.empty((0,N))
    discount_factor_temp_processes = np.empty((0,N))
    discount_factor_processes = np.empty((0,N))

    discount_factor_nodelta_DisSep_Damage_processes = np.empty((0,N))
    discount_factor_nodelta_DisSep_Tech_processes = np.empty((0,N))
    discount_factor_nodelta_processes = np.empty((0,N))
    
    discount_factor_nodeltadt_DisSep_Damage_processes = np.empty((0,N))
    discount_factor_nodeltadt_DisSep_Tech_processes = np.empty((0,N))
    discount_factor_nodeltadt_processes = np.empty((0,N))

    discount_processes = np.empty((0,N))
    discount_o_processes = np.empty((0,N))
    discount_pg_processes = np.empty((0,N))
    discount_pg_f2_processes = np.empty((0,N))
    discount_f_processes = np.empty((0,N))
    discount_f_pg_processes = np.empty((0,N))
    time_processes = np.empty((0,N))

    # first_terms2 = np.empty((0,N))
    # second_terms2 = np.empty((0,N))
    # third_term_1s2 = np.empty((0,N))
    # third_term_2s2 = np.empty((0,N))
    # third_term_3s2 = np.empty((0,N))
    # undiscount_processes2 = np.empty((0,N))
    # discount_factor_temp_processes2 = np.empty((0,N))
    # discount_factor_processes2 = np.empty((0,N))
    # discount_processes2 = np.empty((0,N))




    for i in range(1,101):
    # for i in range(1,51):
        
        # print(Data_Dir + File_Dir+'simulation_'+str(i)+f'_{IntPeriod}_2d_coarse_{m0}.npz')
        # simbatch = np.load(Data_Dir + File_Dir+'simulation_'+str(i)+f'_{IntPeriod}_2d_coarse_{m0}.npz')

        # try:



        # print(i)
        simbatch = np.load(Data_Dir + File_Dir2+'Sim_'+str(i)+f'_{IntPeriod}_{str(dt)}_{m0}_ID.npz')
        W_process = simbatch['W_process']
        Z_process = simbatch['Z_process']
        V_process = simbatch['V_process']
        X_process = simbatch['X_process']
        m1_process = simbatch['m1_process']
        m2_process = simbatch['m2_process']
        m3_process = simbatch['m3_process']
        m4_process = simbatch['m4_process']
        # u1_process = simbatch['u1_process']
        # u2_process = simbatch['u2_process']
        # u3_process = simbatch['u3_process']
        
        # u11_process = simbatch['u11_process']
        # u12_process = simbatch['u12_process']
        # u13_process = simbatch['u13_process']
        
        # u21_process = simbatch['u21_process']
        # u22_process = simbatch['u22_process']
        # u23_process = simbatch['u23_process']


        # u31_process = simbatch['u31_process']
        # u32_process = simbatch['u32_process']
        # u33_process = simbatch['u33_process']
        
        # u41_process = simbatch['u41_process']
        # u42_process = simbatch['u42_process']
        # u43_process = simbatch['u43_process']


        # uF1_process = simbatch['uF1_process']
        # uF2_process = simbatch['uF2_process']
        # uF3_process = simbatch['uF3_process']

        # j1_process = simbatch['j1_process']
        # j2_process = simbatch['j2_process']
        # j3_process = simbatch['j3_process']
        # v1_process = simbatch['v1_process']
        # v2_process = simbatch['v2_process']
        # v3_process = simbatch['v3_process']
        # v_process = simbatch['va_process']
        # dj1_dx1_process = simbatch['dj1_dx1_process']
        # dj2_dx2_process = simbatch['dj2_dx2_process']
        # dj3_dx3_process = simbatch['dj3_dx3_process']
        # dv1_dx1_process = simbatch['dv1_dx1_process']
        # dv1_dx2_process = simbatch['dv1_dx2_process']
        # dv1_dx3_process = simbatch['dv1_dx3_process']
        # dv2_dx1_process = simbatch['dv2_dx1_process']
        # dv2_dx2_process = simbatch['dv2_dx2_process']
        # dv2_dx3_process = simbatch['dv2_dx3_process']
        # dv3_dx1_process = simbatch['dv3_dx1_process']
        # dv3_dx2_process = simbatch['dv3_dx2_process']
        # dv3_dx3_process = simbatch['dv3_dx3_process']
        dv_dx1_process = simbatch['dv_dx1_process']
        dv_dx2_process = simbatch['dv_dx2_process']
        dv_dx3_process = simbatch['dv_dx3_process']
        dv_dx4_process = simbatch['dv_dx4_process']
        first_term = simbatch['first_term']
        second_term = simbatch['second_term']
        third_term_1 = simbatch['third_term_1']
        third_term_2 = simbatch['third_term_2']
        third_term_3 = simbatch['third_term_3']
        third_term_4 = simbatch['third_term_4']
        fourth_term_entropy = simbatch['fourth_term_entropy']
        undiscount_process = simbatch['undiscount_process']
        # undiscount_pg_process = simbatch['undiscount_pg_process']
        # undiscount_pg_f2_process = simbatch['undiscount_pg_f2_process']
        # undiscount_f_process = simbatch['undiscount_f_process']
        # undiscount_f_pg_process = simbatch['undiscount_f_pg_process']
        discount_factor_temp = simbatch['discount_factor_temp']
        discount_factor = simbatch['discount_factor']
        discount_process = simbatch['discount_process']

        discount_factor_nodelta_DisSep_Damage = simbatch['discount_factor_nodelta_DisSep_Damage']
        discount_factor_nodelta_DisSep_Tech = simbatch['discount_factor_nodelta_DisSep_Tech']
        discount_factor_nodelta = simbatch['discount_factor_nodelta']

        discount_factor_nodeltadt_DisSep_Damage = simbatch['discount_factor_nodeltadt_DisSep_Damage']
        discount_factor_nodeltadt_DisSep_Tech = simbatch['discount_factor_nodeltadt_DisSep_Tech']
        discount_factor_nodeltadt = simbatch['discount_factor_nodeltadt']

        time_process = simbatch['time_process']
        # discount_pg_process = simbatch['discount_pg_process']
        # discount_pg_f2_process = simbatch['discount_pg_f2_process']
        # discount_f_process = simbatch['discount_f_process']
        # discount_f_pg_process = simbatch['discount_f_pg_process']
        # first_term2 = simbatch['first_term2']
        # second_term2 = simbatch['second_term2']
        # third_term_12 = simbatch['third_term_12']
        # third_term_22 = simbatch['third_term_22']
        # third_term_32 = simbatch['third_term_32']
        # undiscount_process2 = simbatch['undiscount_process2']
        # discount_factor_temp2 = simbatch['discount_factor_temp2']
        # discount_factor2 = simbatch['discount_factor2']
        # discount_process2 = simbatch['discount_process2']
        W_processes = np.concatenate((W_processes,W_process), axis=0)
        Z_processes = np.concatenate((Z_processes,Z_process), axis=0)
        V_processes = np.concatenate((V_processes,V_process), axis=0)
        X_processes = np.concatenate((X_processes,X_process), axis=0)
        m1_processes = np.concatenate((m1_processes,m1_process), axis=0)
        m2_processes = np.concatenate((m2_processes,m2_process), axis=0)
        m3_processes = np.concatenate((m3_processes,m3_process), axis=0)
        m4_processes = np.concatenate((m4_processes,m4_process), axis=0)
        # u1_processes = np.concatenate((u1_processes,u1_process), axis=0)
        # u2_processes = np.concatenate((u2_processes,u2_process), axis=0)
        # u3_processes = np.concatenate((u3_processes,u3_process), axis=0)
        
        # u11_processes = np.concatenate((u11_processes,u11_process), axis=0)
        # u12_processes = np.concatenate((u12_processes,u12_process), axis=0)
        # u13_processes = np.concatenate((u13_processes,u13_process), axis=0)
        
        # u21_processes = np.concatenate((u21_processes,u21_process), axis=0)
        # u22_processes = np.concatenate((u22_processes,u22_process), axis=0)
        # u23_processes = np.concatenate((u23_processes,u23_process), axis=0)


        # u31_processes = np.concatenate((u31_processes,u31_process), axis=0)
        # u32_processes = np.concatenate((u32_processes,u32_process), axis=0)
        # u33_processes = np.concatenate((u33_processes,u33_process), axis=0)
        
        # u41_processes = np.concatenate((u41_processes,u41_process), axis=0)
        # u42_processes = np.concatenate((u42_processes,u42_process), axis=0)
        # u43_processes = np.concatenate((u43_processes,u43_process), axis=0)


        # uF1_processes = np.concatenate((uF1_processes,uF1_process), axis=0)
        # uF2_processes = np.concatenate((uF2_processes,uF2_process), axis=0)
        # uF3_processes = np.concatenate((uF3_processes,uF3_process), axis=0)

        # j1_processes = np.concatenate((j1_processes,j1_process), axis=0)
        # j2_processes = np.concatenate((j2_processes,j2_process), axis=0)
        # j3_processes = np.concatenate((j3_processes,j3_process), axis=0)
        # v1_processes = np.concatenate((v1_processes,v1_process), axis=0)
        # v2_processes = np.concatenate((v2_processes,v2_process), axis=0)
        # v3_processes = np.concatenate((v3_processes,v3_process), axis=0)
        # v_processes = np.concatenate((v_processes,v_process), axis=0)
        # dj1_dx1_processes = np.concatenate((dj1_dx1_processes,dj1_dx1_process), axis=0)
        # dj2_dx2_processes = np.concatenate((dj2_dx2_processes,dj2_dx2_process), axis=0)
        # dj3_dx3_processes = np.concatenate((dj3_dx3_processes,dj3_dx3_process), axis=0)
        # dv1_dx1_processes = np.concatenate((dv1_dx1_processes,dv1_dx1_process), axis=0)
        # dv1_dx2_processes = np.concatenate((dv1_dx2_processes,dv1_dx2_process), axis=0)
        # dv1_dx3_processes = np.concatenate((dv1_dx3_processes,dv1_dx3_process), axis=0)
        # dv2_dx1_processes = np.concatenate((dv2_dx1_processes,dv2_dx1_process), axis=0)
        # dv2_dx2_processes = np.concatenate((dv2_dx2_processes,dv2_dx2_process), axis=0)
        # dv2_dx3_processes = np.concatenate((dv2_dx3_processes,dv2_dx3_process), axis=0)
        # dv3_dx1_processes = np.concatenate((dv3_dx1_processes,dv3_dx1_process), axis=0)
        # dv3_dx2_processes = np.concatenate((dv3_dx2_processes,dv3_dx2_process), axis=0)
        # dv3_dx3_processes = np.concatenate((dv3_dx3_processes,dv3_dx3_process), axis=0)
        dv_dx1_processes = np.concatenate((dv_dx1_processes,dv_dx1_process), axis=0)
        dv_dx2_processes = np.concatenate((dv_dx2_processes,dv_dx2_process), axis=0)
        dv_dx3_processes = np.concatenate((dv_dx3_processes,dv_dx3_process), axis=0)
        dv_dx4_processes = np.concatenate((dv_dx4_processes,dv_dx4_process), axis=0)
        first_terms = np.concatenate((first_terms,first_term), axis=0)
        second_terms = np.concatenate((second_terms,second_term), axis=0)
        third_term_1s = np.concatenate((third_term_1s,third_term_1), axis=0)
        third_term_2s = np.concatenate((third_term_2s,third_term_2), axis=0)
        third_term_3s = np.concatenate((third_term_3s,third_term_3), axis=0)
        third_term_4s = np.concatenate((third_term_4s,third_term_4), axis=0)
        fourth_term_entropys = np.concatenate((fourth_term_entropys,fourth_term_entropy), axis=0)
        undiscount_processes = np.concatenate((undiscount_processes,undiscount_process), axis=0)
        # undiscount_pg_processes = np.concatenate((undiscount_pg_processes,undiscount_pg_process), axis=0)
        # undiscount_pg_f2_processes = np.concatenate((undiscount_pg_f2_processes,undiscount_pg_f2_process), axis=0)
        # undiscount_f_processes = np.concatenate((undiscount_f_processes,undiscount_f_process), axis=0)
        # undiscount_f_pg_processes = np.concatenate((undiscount_f_pg_processes,undiscount_f_pg_process), axis=0)
        discount_factor_temp_processes = np.concatenate((discount_factor_temp_processes,discount_factor_temp), axis=0)
        discount_factor_processes = np.concatenate((discount_factor_processes,discount_factor), axis=0)

        discount_factor_nodelta_DisSep_Damage_processes = np.concatenate((discount_factor_nodelta_DisSep_Damage_processes,discount_factor_nodelta_DisSep_Damage), axis=0)
        discount_factor_nodelta_DisSep_Tech_processes = np.concatenate((discount_factor_nodelta_DisSep_Tech_processes,discount_factor_nodelta_DisSep_Tech), axis=0)
        discount_factor_nodelta_processes = np.concatenate((discount_factor_nodelta_processes,discount_factor_nodelta), axis=0)

        discount_factor_nodeltadt_DisSep_Damage_processes = np.concatenate((discount_factor_nodeltadt_DisSep_Damage_processes,discount_factor_nodeltadt_DisSep_Damage), axis=0)
        discount_factor_nodeltadt_DisSep_Tech_processes = np.concatenate((discount_factor_nodeltadt_DisSep_Tech_processes,discount_factor_nodeltadt_DisSep_Tech), axis=0)
        discount_factor_nodeltadt_processes = np.concatenate((discount_factor_nodeltadt_processes,discount_factor_nodeltadt), axis=0)

        discount_processes = np.concatenate((discount_processes,discount_process), axis=0)
        # discount_pg_processes = np.concatenate((discount_pg_processes,discount_pg_process), axis=0)
        # discount_pg_f2_processes = np.concatenate((discount_pg_f2_processes,discount_pg_f2_process), axis=0)
        # discount_f_processes = np.concatenate((discount_f_processes,discount_f_process), axis=0)
        # discount_f_pg_processes = np.concatenate((discount_f_pg_processes,discount_f_pg_process), axis=0)

        # first_terms2 = np.concatenate((first_terms2,first_term2), axis=0)
        # second_terms2 = np.concatenate((second_terms2,second_term2), axis=0)
        # third_term_1s2 = np.concatenate((third_term_1s2,third_term_12), axis=0)
        # third_term_2s2 = np.concatenate((third_term_2s2,third_term_22), axis=0)
        # third_term_3s2 = np.concatenate((third_term_3s2,third_term_32), axis=0)
        # undiscount_processes2 = np.concatenate((undiscount_processes2,undiscount_process2), axis=0)
        # discount_factor_temp_processes2 = np.concatenate((discount_factor_temp_processes2,discount_factor_temp2), axis=0)
        # discount_factor_processes2 = np.concatenate((discount_factor_processes2,discount_factor2), axis=0)
        # discount_processes2 = np.concatenate((discount_processes2,discount_process2), axis=0)
        # time_processes.append(time_process)




        # except:
        #     print(f"error loading {i}")

    # plt.figure()
    # plt.plot(np.mean(undiscount_processes, axis=0),label = r'$\delta m \cdot \frac{\partial U}{\partial x} + m \cdot \frac{\partial J}{\partial x} (V^j - V) + m \cdot J \frac{\partial V^j}{\partial x} $')
    # plt.legend(loc='lower left')
    # plt.title(str(m0)+" "+labellist[id])
    # plt.savefig(Plot_Dir+"/"+Filename+labellist[id]+str(m0)+"_UndiscountedProcess.pdf")
    # plt.close()


    # plt.figure()
    # plt.plot( np.mean(discount_processes, axis=0),label = r'$\exp( \int_0^t(-\delta - J)du )(\delta m \cdot \frac{\partial U}{\partial x} + m \cdot \frac{\partial J}{\partial x} (V^j - V) + m \cdot J \frac{\partial V^j}{\partial x}) $')
    # plt.legend(loc='lower left')
    # plt.title(str(m0)+" "+labellist[id])
    # if m0=="Technology":
    #     # plt.ylim(0,0.002)
    #     print("ylim Applied Success")
    # if m0=="Temperature":
    #     plt.ylim(-0.005,0)
    #     plt.xlim(0,50)
    #     print("ylim Applied Success")
    # plt.savefig(Plot_Dir+"/"+Filename+labellist[id]+str(m0)+"_DiscountedProcess.pdf")
    # plt.close()

    # plt.figure()
    # plt.plot( np.mean(discount_pg_processes, axis=0),label = r'$\exp( \int_0^t(-\delta - J)du )(\delta m \cdot \frac{\partial U}{\partial x} + m \cdot \frac{\partial J}{\partial x} (V^j - V) + m \cdot J \frac{\partial V^j}{\partial x}) + m \cdot J \frac{\partial g}{\partial x} (V^j - V)$')
    # plt.legend(loc='lower left')
    # plt.title(str(m0)+" "+labellist[id])
    # if m0=="Technology":
    #     plt.ylim(0,0.002)
    #     print("ylim Applied Success")
    # if m0=="Temperature":
    #     plt.ylim(-0.005,0)
    #     plt.xlim(0,50)
    #     print("ylim Applied Success")
    # plt.savefig(Plot_Dir+"/"+Filename+labellist[id]+str(m0)+"_DiscountedProcess_partialg.pdf")
    # plt.close()

    # plt.figure()
    # plt.plot( np.mean(discount_pg_f2_processes, axis=0),label = r'$\exp( \int_0^t(-\delta - J)du )(\delta m \cdot \frac{\partial U}{\partial x} + m \cdot \frac{\partial J}{\partial x} (V^j - V) + m \cdot J \frac{\partial V^j}{\partial x}) + m \cdot J \frac{\partial g}{\partial x} (V^j - V) + \frac{\partial h^2}{\partial x}$')
    # plt.legend(loc='lower left')
    # plt.title(str(m0)+" "+labellist[id])
    # if m0=="Technology":
    #     plt.ylim(0,0.002)
    #     print("ylim Applied Success")
    # if m0=="Temperature":
    #     plt.ylim(-0.005,0)
    #     plt.xlim(0,50)
    #     print("ylim Applied Success")
    # plt.savefig(Plot_Dir+"/"+Filename+labellist[id]+str(m0)+"_DiscountedProcess_partialg_f2.pdf")
    # plt.close()

    # plt.figure()
    # plt.plot( np.mean(discount_pg_f2_processes, axis=0),label = r'Draft + Jump Derivative + Entropy Derivative')
    # plt.plot( np.mean(discount_pg_processes, axis=0),label = r'Draft + Jump Derivative')
    # plt.plot( np.mean(discount_processes, axis=0),label = r'Draft')
    # plt.legend(loc='lower left')
    # plt.title(str(m0)+" "+labellist[id])
    # if m0=="Technology":
    #     plt.ylim(0,0.002)
    #     print("ylim Applied Success")
    # if m0=="Temperature":
    #     plt.ylim(-0.005,0)
    #     plt.xlim(0,50)
    #     print("ylim Applied Success")
    # plt.savefig(Plot_Dir+"/"+Filename+labellist[id]+str(m0)+"_DiscountedProcess_Compare.pdf")
    # plt.close()


    # plt.figure()
    # plt.plot(np.mean(delta*first_terms, axis=0),label = r'$\delta m \cdot \frac{\partial U}{\partial x}$')
    # plt.plot(np.mean(second_terms, axis=0),label = r'$m \cdot \frac{\partial J}{\partial x} (V^j - V)$')
    # # plt.plot(np.mean(third_term_1s, axis=0),label = 'third term 1')
    # # plt.plot(np.mean(third_term_2s, axis=0),label = 'third term 2')
    # # plt.plot(np.mean(third_term_3s, axis=0),label = 'third term 3'
    # plt.plot(np.mean(third_term_1s, axis=0)+np.mean(third_term_2s, axis=0)+np.mean(third_term_3s, axis=0)+np.mean(third_term_4s, axis=0),label = r'$m \cdot J \frac{\partial V^j}{\partial x} $')
    # plt.plot(np.mean(fourth_term_entropys, axis=0),label = r'$\xi m \cdot \frac{\partial J}{\partial x} (1-g^* + g^* \log g^*)$')
    # # plt.plot(,label = 'third term 3')
    # plt.legend(loc='upper left')
    # plt.title(str(m0)+" "+labellist[id])
    # plt.savefig(Plot_Dir+"/"+Filename+labellist[id]+str(m0)+"_Undiscount_Term1234_decomposition.pdf")
    # plt.close()


    # plt.figure()
    # plt.plot(np.mean(np.exp(discount_factor_processes), axis=0),label = r'$\exp(\int_0^t (-\delta-J) dt)$')
    # # plt.plot(,label = 'third term 3')
    # plt.legend(loc='upper right')
    # # plt.ylim(0,0.001)
    # plt.title(str(m0)+" "+labellist[id])
    # plt.savefig(Plot_Dir+"/"+Filename+labellist[id]+str(m0)+"_Discount_Factor.pdf")
    # plt.close()


    # plt.figure()
    # plt.plot(np.mean(np.exp(discount_factor_processes)*(delta*first_terms), axis=0),label = r'$\delta m \cdot \frac{\partial U}{\partial x}$')
    # plt.plot(np.mean(np.exp(discount_factor_processes)*(second_terms), axis=0),label = r'$m \cdot \frac{\partial J}{\partial x} (V^j - V)$')
    # # plt.plot(np.mean(third_term_1s, axis=0),label = 'third term 1')
    # # plt.plot(np.mean(third_term_2s, axis=0),label = 'third term 2')
    # # plt.plot(np.mean(third_term_3s, axis=0),label = 'third term 3'
    # # plt.plot(np.mean(third_term_1s, axis=0)+np.mean(third_term_2s, axis=0)+np.mean(third_term_3s, axis=0)+np.mean(third_term_4s, axis=0),label = r'$m \cdot J \frac{\partial V^j}{\partial x} $')
    # plt.plot(np.mean(np.exp(discount_factor_processes)*(third_term_1s+third_term_2s+third_term_3s+third_term_4s), axis=0),label = r'$m \cdot J \frac{\partial V^j}{\partial x} $')
    # plt.plot(np.mean(np.exp(discount_factor_processes)*(fourth_term_entropys), axis=0),label = r'$\xi m \cdot \frac{\partial J}{\partial x} (1-g^* + g^* \log g^*)$')
    # # plt.plot(,label = 'third term 3')
    # plt.legend(loc='upper right')
    # plt.ylim(0,0.0012)
    # plt.title(str(m0)+" "+labellist[id])
    # plt.savefig(Plot_Dir+"/"+Filename+labellist[id]+str(m0)+"_Discount_Term1234_decomposition.pdf")
    # plt.close()

    # plt.figure()
    # plt.plot(np.mean(W_processes, axis=0),label = r'$\log K$')
    # plt.plot(np.mean(Z_processes, axis=0),label = 'Y')
    # plt.plot(np.mean(V_processes, axis=0),label = r'$\log Ig$')    # plt.plot(,label = 'third term 2')
    # plt.plot(np.mean(X_processes, axis=0),label = r'$\log N$')    # plt.plot(,label = 'third term 2')
    # # plt.plot(,label = 'third term 3')
    # plt.legend()
    # plt.title(str(m0)+" "+labellist[id])
    # plt.savefig(Plot_Dir+"/"+Filename+labellist[id]+str(m0)+"_StateProcesss.pdf")
    # plt.close()


    # plt.figure()
    # plt.plot(np.mean(m1_processes, axis=0),label = r'$m_1$')
    # plt.plot(np.mean(m2_processes, axis=0),label = r'$m_2$')
    # plt.plot(np.mean(m3_processes, axis=0),label = r'$m_3$')
    # plt.plot(np.mean(m4_processes, axis=0),label = r'$m_4$')
    # # plt.plot(,label = 'third term 3')
    # plt.legend()
    # plt.title(str(m0)+" "+labellist[id])
    # if m0=="Temperature":
    #     plt.xlim(0,50)
    # plt.savefig(Plot_Dir+"/"+Filename+labellist[id]+str(m0)+"_MProcesss.pdf")
    # plt.close()

    # plt.figure()
    # plt.plot(np.mean(u1_processes, axis=0),label = 'ux1')
    # plt.plot(np.mean(u2_processes, axis=0),label = 'ux2')
    # plt.plot(np.mean(u3_processes, axis=0),label = 'ux3')
    # # plt.plot(,label = 'third term 3')
    # plt.legend()
    # plt.title(str(m0)+r'$\xi=\infty$')
    # plt.savefig(Plot_Dir+"/"+Filename+labellist[id]+str(m0)+"_PartialUProcesss.pdf")
    # plt.close()

 
    # plt.figure()
    # plt.plot(np.mean(dv_dx1_processes, axis=0),label = 'dv_x1')
    # plt.plot(np.mean(dv_dx2_processes, axis=0),label = 'dv_x2')
    # plt.plot(np.mean(dv_dx3_processes, axis=0),label = 'dv_x3')
    # # plt.plot(,label = 'third term 3')
    # plt.legend()
    # plt.title(str(m0)+r'$\xi=\infty$')
    # plt.savefig(Plot_Dir+"/"+Filename+labellist[id]+str(m0)+"_PartialVProcesss.pdf")
    # plt.close()

    # plt.figure()
    # plt.plot(np.mean(u1_processes*m1_processes +u2_processes*m2_processes + u3_processes*m3_processes, axis=0),label = 'u')
    # plt.plot(np.mean(u11_processes*m1_processes +u12_processes*m2_processes + u13_processes*m3_processes, axis=0),label = 'u1')
    # plt.plot(np.mean(u21_processes*m1_processes +u22_processes*m2_processes + u23_processes*m3_processes, axis=0),label = 'u2')
    # plt.plot(np.mean(u31_processes*m1_processes +u32_processes*m2_processes + u33_processes*m3_processes, axis=0),label = 'u3')
    # plt.plot(np.mean(u41_processes*m1_processes +u42_processes*m2_processes + u43_processes*m3_processes, axis=0),label = 'u4')
    # # plt.plot(,label = 'third term 3')
    # plt.legend()
    # plt.title(str(m0)+r'$\xi=\infty$')
    # plt.savefig(Plot_Dir+"/"+Filename+labellist[id]+str(m0)+"_UndiscountU_U1_U2_Contribution_Processs.pdf")
    # plt.close()
    LHS = np.mean(dv_dx1_processes*m1_processes +dv_dx2_processes*m2_processes + dv_dx3_processes*m3_processes +dv_dx4_processes*m4_processes,axis=0)[0]
    # RHS = np.sum(np.mean(discount_processes, axis=0))
    # RHS_1 = np.sum(np.mean(np.exp(discount_factor_processes)*(delta*first_terms), axis=0))
    # RHS_2 = np.sum(np.mean(np.exp(discount_factor_processes)*(second_terms), axis=0))
    # RHS_3 = np.sum(np.mean(np.exp(discount_factor_processes)*(third_term_1s+third_term_2s+third_term_3s+third_term_4s), axis=0))
    # RHS_4 = np.sum(np.mean(np.exp(discount_factor_processes)*(fourth_term_entropys), axis=0))
    RHS = np.trapz(np.mean(discount_processes, axis=0), time_process)
    RHS_1 = np.trapz(np.mean(np.exp(discount_factor_processes)*(delta*first_terms), axis=0), time_process)
    RHS_2 = np.trapz(np.mean(np.exp(discount_factor_processes)*(second_terms), axis=0), time_process)
    RHS_3 = np.trapz(np.mean(np.exp(discount_factor_processes)*(third_term_1s+third_term_2s+third_term_3s+third_term_4s), axis=0), time_process)
    RHS_3a = np.trapz(np.mean(np.exp(discount_factor_processes)*(third_term_2s), axis=0), time_process)
    RHS_3b = np.trapz(np.mean(np.exp(discount_factor_processes)*(third_term_3s), axis=0), time_process)
    RHS_4 = np.trapz(np.mean(np.exp(discount_factor_processes)*(fourth_term_entropys), axis=0), time_process)
    # print(time_process.shape)
    Integral = np.trapz(np.mean(discount_factor_nodeltadt_processes, axis=0), time_process)
    Integral_d = np.trapz(np.mean(discount_factor_nodeltadt_DisSep_Damage_processes, axis=0), time_process)
    Integral_t = np.trapz(np.mean(discount_factor_nodeltadt_DisSep_Tech_processes, axis=0), time_process)

    print("###########################################")
    print("Start######################################")
    print(f"shape = {discount_factor_nodeltadt_processes.shape}")
    print(f"RHS_Undiscount = {np.sum(np.mean(undiscount_processes, axis=0))}, RHS_Discount = {RHS}, LHS = {LHS}, Difference (RHS-LHS)/LHS %= {(RHS-LHS)/LHS*100}")
    print(f"RHS_Discounted_1={RHS_1}, RHS_Discounted_2={RHS_2}, RHS_Discounted_3={RHS_3}, RHS_Discounted_3a={RHS_3a}, RHS_Discounted_3b={RHS_3b}, RHS_Discounted_4={RHS_4}")
    print(f"Integral of Time derivative={Integral}, Damage Integral={Integral_d}, Tech Integral={Integral_t}")
    # RHS = np.sum(np.mean(discount_pg_processes, axis=0))
    # print(f"RHS_Undiscount_partial g = {np.sum(np.mean(undiscount_pg_processes, axis=0))}, RHS_Discount = {RHS}, LHS = {LHS}, Difference (RHS-LHS)/LHS %= {(RHS-LHS)/LHS*100}")
    # RHS = np.sum(np.mean(discount_pg_f2_processes, axis=0))
    # print(f"RHS_Undiscount_partial g f2= {np.sum(np.mean(undiscount_pg_f2_processes, axis=0))}, RHS_Discount = {RHS}, LHS = {LHS}, Difference (RHS-LHS)/LHS %= {(RHS-LHS)/LHS*100}")
    # RHS = np.sum(np.mean(discount_f_processes, axis=0))
    # print(f"RHS_Undiscount_Only U Full = {np.sum(np.mean(undiscount_f_processes, axis=0))}, RHS_Discount = {RHS}, LHS = {LHS}, Difference (RHS-LHS)/LHS %= {(RHS-LHS)/LHS*100}")
    # RHS = np.sum(np.mean(discount_f_pg_processes, axis=0))
    # print(f"RHS_Undiscount_U full+npartialg = {np.sum(np.mean(undiscount_f_pg_processes, axis=0))}, RHS_Discount = {RHS}, LHS = {LHS}, Difference (RHS-LHS)/LHS %= {(RHS-LHS)/LHS*100}")
    # print(f"RHS_Undiscount_Method2 = {np.sum(np.mean(undiscount_processes2, axis=0))}, RHS_Discount_Method2 = {np.sum(np.mean(discount_processes2, axis=0))}, LHS = {np.mean(dv_dx1_processes*m1_processes +dv_dx2_processes*m2_processes + dv_dx3_processes*m3_processes,axis=0)[0]}")
    print(f"dvdY = {np.mean(dv_dx2_processes,axis=0)[0]}, dvdlogIg = {np.mean(dv_dx3_processes,axis=0)[0]}")
    

    Second_Term_0 = np.mean(np.exp(discount_factor_processes)*(second_terms), axis=0)[0]
    Second_Term_20 = np.mean(np.exp(discount_factor_processes)*(second_terms), axis=0)[int(20/dt)]
    Second_Term_25 = np.mean(np.exp(discount_factor_processes)*(second_terms), axis=0)[int(25/dt)]
    Second_Term_30 = np.mean(np.exp(discount_factor_processes)*(second_terms), axis=0)[int(30/dt)]
    print(f"Second Year 0={Second_Term_0}, 20={Second_Term_20}, 25={Second_Term_25}, 30={Second_Term_30}, Year 20 Ratio={Second_Term_20/Second_Term_0}, Year 25 Ratio={Second_Term_25/Second_Term_0}, Year 30 Ratio={Second_Term_30/Second_Term_0}")
    print("#########################################")
    print("End######################################")

    return np.mean(np.exp(discount_factor_processes), axis=0), np.mean(np.exp(discount_factor_nodelta_processes),axis=0), np.mean(discount_factor_nodeltadt_processes,axis=0), np.mean(discount_factor_nodelta_processes,axis=0),time_process,np.mean(np.exp(discount_factor_nodelta_DisSep_Damage_processes),axis=0),np.mean(np.exp(discount_factor_nodelta_DisSep_Tech_processes),axis=0),np.mean(discount_factor_nodeltadt_DisSep_Damage_processes,axis=0),np.mean(discount_factor_nodeltadt_DisSep_Tech_processes,axis=0)

res_list = []
res2_list = []
res3_list = []
res4_list = []

res5_list = []
res6_list = []
res7_list = []
res8_list = []

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):
                for id_rho in range(len(rhoarr)):

                    res,  res2, res3, res4, time_process, res_ProbDamage, res_ProbTech, res_probdamage, res_probtech = model_simulation_generate(id_xiag, xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],xia2arr[id_xiag],xik2arr[id_xiag],xic2arr[id_xiag],xij2arr[id_xiag],xid2arr[id_xiag],xig2arr[id_xiag],rhoarr[id_rho],psi0arr[id_psi0],psi1arr[id_psi1],varrhoarr[id_varrho])
                    
                    res_list.append(res)
                    res2_list.append(res2)
                    res3_list.append(res3)
                    res4_list.append(res4)

                    res5_list.append(res_ProbDamage)
                    res6_list.append(res_ProbTech)
                    res7_list.append(res_probdamage)
                    res8_list.append(res_probtech)

plt.figure()
for id_xiag in range(len(xiaarr)): 
    plt.plot(time_process,res_list[id_xiag],label=labellist[id_xiag])
plt.legend(loc='upper right')
plt.savefig(Plot_Dir+"/"+Filename+"Composite"+str(m0)+"{:.3f}".format(dt)+"_DiscountProcesss_dt.pdf")
plt.close()


plt.figure()
for id_xiag in range(len(xiaarr)): 
    plt.plot(time_process,res2_list[id_xiag],label=labellist[id_xiag])
plt.legend(loc='upper right')
plt.savefig(Plot_Dir+"/"+Filename+"Composite"+str(m0)+"{:.3f}".format(dt)+"_Discount_nodelta_Processs_dt.pdf")
plt.close()


plt.figure()
for id_xiag in range(len(xiaarr)): 
    plt.plot(time_process,res3_list[id_xiag],label=labellist[id_xiag])
plt.legend(loc='upper right')
plt.ylim(0, 0.06)
plt.savefig(Plot_Dir+"/"+Filename+"Composite"+str(m0)+"{:.3f}".format(dt)+"_Discount_nodeltadt_Processs_dt.pdf")
plt.close()


plt.figure()
for id_xiag in range(len(xiaarr)): 
    plt.plot(time_process,res4_list[id_xiag],label=labellist[id_xiag])
plt.legend(loc='upper right')
plt.savefig(Plot_Dir+"/"+Filename+"Composite"+str(m0)+"{:.3f}".format(dt)+"_Discount_nodelta_noexp_Processs_dt.pdf")
plt.close()


plt.figure()
for id_xiag in range(len(xiaarr)): 
    plt.plot(time_process,res5_list[id_xiag],label=labellist[id_xiag])
plt.legend(loc='upper right')
plt.ylim(0,1)
plt.savefig(Plot_Dir+"/"+Filename+"Composite"+str(m0)+"{:.3f}".format(dt)+"_Discount_nodelta_damage_Processs_dt.pdf")
plt.close()

plt.figure()
for id_xiag in range(len(xiaarr)): 
    plt.plot(time_process,res6_list[id_xiag],label=labellist[id_xiag])
plt.legend(loc='upper right')
plt.ylim(0,1)
plt.savefig(Plot_Dir+"/"+Filename+"Composite"+str(m0)+"{:.3f}".format(dt)+"_Discount_nodelta_tech_Processs_dt.pdf")
plt.close()

plt.figure()
for id_xiag in range(len(xiaarr)): 
    plt.plot(time_process,res7_list[id_xiag],label=labellist[id_xiag])
plt.legend(loc='upper right')
plt.ylim(0, 0.06)
plt.savefig(Plot_Dir+"/"+Filename+"Composite"+str(m0)+"{:.3f}".format(dt)+"_Discount_nodeltadt_damage_Processs_dt.pdf")
plt.close()

plt.figure()
for id_xiag in range(len(xiaarr)): 
    plt.plot(time_process,res8_list[id_xiag],label=labellist[id_xiag])
plt.legend(loc='upper right')
plt.ylim(0, 0.06)
plt.savefig(Plot_Dir+"/"+Filename+"Composite"+str(m0)+"{:.3f}".format(dt)+"_Discount_nodeltadt_tech_Processs_dt.pdf")
plt.close()

plt.figure()
for id_xiag in range(len(xiaarr)): 
    # print(res5_list[id_xiag].shape, -finite_1D(res5_list[id_xiag],dt).shape)
    plt.plot(time_process,-finite_1D(res2_list[id_xiag],dt),label=labellist[id_xiag])
plt.legend(loc='upper right')
plt.ylim(0, 0.06)
plt.savefig(Plot_Dir+"/"+Filename+"Composite"+str(m0)+"{:.3f}".format(dt)+"_Discount_nodeltadt_Processs_numerical_dt.pdf")
plt.close()


plt.figure()
for id_xiag in range(len(xiaarr)): 
    # print(res5_list[id_xiag].shape, -finite_1D(res5_list[id_xiag],dt).shape)
    plt.plot(time_process,-finite_1D(res5_list[id_xiag],dt),label=labellist[id_xiag])
plt.legend(loc='upper right')
plt.ylim(0, 0.06)
plt.savefig(Plot_Dir+"/"+Filename+"Composite"+str(m0)+"{:.3f}".format(dt)+"_Discount_nodeltadt_damage_Processs_numerical_dt.pdf")
plt.close()

plt.figure()
for id_xiag in range(len(xiaarr)): 
    plt.plot(time_process,-finite_1D(res6_list[id_xiag],dt),label=labellist[id_xiag])
plt.ylim(0, 0.06)
plt.legend(loc='upper right')
plt.savefig(Plot_Dir+"/"+Filename+"Composite"+str(m0)+"{:.3f}".format(dt)+"_Discount_nodeltadt_tech_Processs_numerical_dt.pdf")
plt.close()