"""
post_damage.py
======================
Solver for solving post damage HJBs, with different values of gamma_3 

python3 -u /home/bcheng4/TwoCapital_Shrink/abatement_UD/postdamage_2jump_CRS.py --num_gamma 3 --xi_a 0.0002 --xi_g 0.025  --epsilonarr 0.1 0.01  --fractionarr 0.1 0.01   --maxiterarr 80000 200000  --id 3 --psi_0 0.105830 --psi_1 0.5 --name 2jump_step_4.00,9.00_0.0,4.0_1.0,6.0_SS_0.2,0.2,0.2_LR_0.1_CRS_PETSCFK --hXarr 0.2 0.2 0.2 --Xminarr 4.00 0.0 1.0 0.0 --Xmaxarr 9.00 4.0 6.0 3.0


"""
# Optimization of post jump HJB
# Required packages
import os
import sys
sys.path.append('./src')
import csv
from src.Utility import *
from src.Utility import finiteDiff_3D
sys.stdout.flush()
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import petsclinearsystem
from scipy.sparse import spdiags
from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix
from datetime import datetime
# from solver import solver_3d
from src.PostSolver_new import hjb_post_tech
from src.PreSolver_CRS2_new import hjb_pre_tech
from src.PostSolver_new_rho1 import hjb_post_tech as hjb_post_tech_rho1
from src.PreSolver_CRS2_new_rho1 import hjb_pre_tech as hjb_pre_tech_rho1
import argparse

reporterror = True
# Linear solver choices
# Chosse among petsc, petsc4py, eigen, both
# petsc: matrix assembled in C
# petsc4py: matrix assembled in Python
# eigen: matrix assembled in C++
# both: petsc+petsc4py
#
now = datetime.now()
current_time = now.strftime("%d-%H:%M")


parser = argparse.ArgumentParser(description="xi_r values")
parser.add_argument("--xi_a", type=float, default=1000.)
parser.add_argument("--xi_k", type=float, default=1000.)
parser.add_argument("--xi_c", type=float, default=1000.)
parser.add_argument("--xi_j", type=float, default=1000.)
parser.add_argument("--xi_d", type=float, default=1000.)
parser.add_argument("--xi_g", type=float, default=1000.)
parser.add_argument("--varrho", type=float, default=1000.)
parser.add_argument("--phi_0", type=float)
parser.add_argument("--rho", type=float)
parser.add_argument("--delta", type=float)
parser.add_argument("--id", type=int, default=0)
parser.add_argument("--id_sub", type=int)
parser.add_argument("--psi_0", type=float, default=0.003)
parser.add_argument("--psi_1", type=float, default=0.5)
parser.add_argument("--num_gamma",type=int,default=6)
parser.add_argument("--name",type=str,default="ReplicateSuri")
parser.add_argument("--outputname",type=str,default="ReplicateSuri")
parser.add_argument("--hXarr",nargs='+',type=float)
parser.add_argument("--Xminarr",nargs='+',type=float)
parser.add_argument("--Xmaxarr",nargs='+',type=float)
parser.add_argument("--epsilonarr",nargs='+',type=float,default=(0.1))
parser.add_argument("--fractionarr",nargs='+',type=float, default=(0.1, 0.1))
parser.add_argument("--maxiterarr",nargs='+',type=int, default=(80000, 200000))

args = parser.parse_args()
# args, unknown = parser.parse_known_args()

epsilonarr = args.epsilonarr
fractionarr = args.fractionarr
maxiterarr = args.maxiterarr


start_time = time.time()
# Parameters as defined in the paper
xi_a = args.xi_a  # Smooth ambiguity
xi_b = 1000. # Brownian misspecification
xi_k = args.xi_k  # Technology jump
xi_c = args.xi_c  # Technology jump
xi_j = args.xi_j  # Technology jump
xi_d = args.xi_d # Hold place for arguments, no real effects 
xi_g = args.xi_g # Hold place for arguments, no real effects 
varrho = args.varrho # Hold place for arguments, no real effects 
rho = args.rho

# Model parameters
# delta   = 0.010
delta = args.delta
alpha   = 0.115
kappa   = 6.667
mu_k    = -0.043
# sigma_k = np.sqrt(0.0087**2 + 0.0038**2)
sigma_k = 0.0100

# Technology
theta        = 3
lambda_bar   = 0.1206
# vartheta_bar = 0.0453
# vartheta_bar = 0.05
# vartheta_bar = 0.056
# vartheta_bar = 0.5
vartheta_bar = args.phi_0

# Damage function
gamma_1 = 1.7675/10000
gamma_2 = 0.0022 * 2
# gamma_3 = 0.3853 * 2

# num_gamma = args.num_gamma
# gamma_3_list = np.linspace(0,1./3.,num_gamma)

num_gamma = args.num_gamma
gamma_3_list = np.linspace(0,1./3.,num_gamma)

id_damage = args.id
id_damage_sub = args.id_sub

gamma_3_i = gamma_3_list[id_damage]
gamma_3_i_sub = gamma_3_list[id_damage_sub]
# gamma_3_list = np.array([0.])
y_bar = 2.
y_bar_lower = 1.5


theta_ell = pd.read_csv('./data/model144.csv', header=None).to_numpy()[:, 0]/1000.
pi_c_o    = np.ones_like(theta_ell)/len(theta_ell)
sigma_y   = 1.2 * np.mean(theta_ell)
beta_f    = np.mean(theta_ell)
# Jump intensity
zeta      = 0.00
psi_0     = args.psi_0
psi_1     = args.psi_1
# sigma_g   = 0.016
sigma_g   = 0.0078
# Tech jump
lambda_bar_first = lambda_bar / 2
vartheta_bar_first = vartheta_bar / 2
lambda_bar_second = 1e-9
vartheta_bar_second = 0.

Xminarr = args.Xminarr
Xmaxarr = args.Xmaxarr
hXarr = args.hXarr

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

# Output_Dir = "/scratch/bincheng/"
Output_Dir = args.outputname
# Output_Dir = "./TwoCapital_Shrink/"
Data_Dir = Output_Dir+"abatement/data_2tech/"+args.name+"/"

File_Name = "xi_a_{}_xi_k_{}_xi_c_{}_xi_j_{}_xi_d_{}_xi_g_{}_psi_0_{}_psi_1_{}_varrho_{}_rho_{}_delta_{}_" .format(xi_a, xi_k, xi_c, xi_j, xi_d, xi_g, psi_0,psi_1, varrho, rho, delta)
#if not os.path.exists(DataDir):
#    os.mkdir(DataDir)

os.makedirs(Data_Dir, exist_ok=True)

# filename =  "post_damage_" + str(gamma_3)  + '_{}'.format(current_time)
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


# Post damage, tech II
print("-------------------------------------------")
print("------------Post damage, Tech II----------")
print("-------------------------------------------")

# model_tech2_post_damage = pickle.load(open(Data_Dir + File_Name + "model_tech2_post_damage_gamma_{:.4f}".format(gamma_3_i), "rb"))

# Post damage, tech II
pi_c = np.array([temp * np.ones(K_mat.shape) for temp in pi_c_o])
pi_c_o = pi_c.copy()
theta_ell = np.array([temp * np.ones(K_mat.shape) for temp in theta_ell])

# v = model_tech2_post_damage["v"]
Guess = pickle.load(open(Data_Dir + File_Name + "model_tech2_post_damage_gamma_{:.4f}".format(gamma_3_i_sub), "rb"))

# Guess = None


if rho==1:
    print("hjb_post_damage_post_tech_rho1")
    
    model_tech2_post_damage = hjb_post_tech_rho1(
            state_grid=(K, Y, L), 
            model_args=(delta, alpha, theta, vartheta_bar, lambda_bar, mu_k, kappa, sigma_k, theta_ell, pi_c_o, pi_c, sigma_y, zeta, psi_0, psi_1, sigma_g, gamma_1, gamma_2, gamma_3_i, y_bar, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g,rho, varrho),
            V_post_damage=None,
            tol=1e-7, epsilon=epsilonarr[0], fraction=fractionarr[0], 
            smart_guess=Guess, 
            max_iter=maxiterarr[0],
            )
else:
    print("hjb_post_damage_post_tech")
    # model_args = (delta, alpha, kappa, mu_k, sigma_k, theta_ell, pi_c_o, sigma_y, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g,rho, gamma_1, gamma_2, gamma_3_i, y_bar, theta, lambda_bar_second, vartheta_bar_second)
    # model_args=(delta, alpha, theta, vartheta_bar, lambda_bar, mu_k, kappa, sigma_k, theta_ell, pi_c_o, pi_c, sigma_y, zeta, psi_0, psi_1, sigma_g, V_post_tech2, gamma_1, gamma_2, gamma_3_i, y_bar, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g, rho, varrho)
    
    model_tech2_post_damage = hjb_post_tech(
            state_grid=(K, Y, L), 
            model_args=(delta, alpha, theta, vartheta_bar, lambda_bar, mu_k, kappa, sigma_k, theta_ell, pi_c_o, pi_c, sigma_y, zeta, psi_0, psi_1, sigma_g, gamma_1, gamma_2, gamma_3_i, y_bar, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g,rho, varrho),
            V_post_damage=None,
            tol=1e-7, epsilon=epsilonarr[0], fraction=fractionarr[0], 
            smart_guess=Guess, 
            max_iter=maxiterarr[0],
            )
    
with open(Data_Dir+ File_Name + "model_tech2_post_damage_gamma_{:.4f}".format(gamma_3_i), "wb") as f:
   pickle.dump(model_tech2_post_damage, f)



print("--------------------------------")
print("-----------------Checking Post Tech---------------")
print("--------------------------------")

model_tech2_post_damage = pickle.load(open(Data_Dir + File_Name + "model_tech2_post_damage_gamma_{:.4f}".format(gamma_3_i), "rb"))

Guess = model_tech2_post_damage

if rho==1:
    print("hjb_post_damage_post_tech_rho1")
    model_args = (delta, alpha, kappa, mu_k, sigma_k, theta_ell, pi_c_o, sigma_y, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g, gamma_1, gamma_2, gamma_3_i, y_bar, theta, lambda_bar_second, vartheta_bar_second)

    model_tech2_post_damage = hjb_post_tech_rho1(
            state_grid=(K, Y, L), 
            model_args=(delta, alpha, theta, vartheta_bar, lambda_bar, mu_k, kappa, sigma_k, theta_ell, pi_c_o, pi_c, sigma_y, zeta, psi_0, psi_1, sigma_g, gamma_1, gamma_2, gamma_3_i, y_bar, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g,rho, varrho),
            V_post_damage=None,
            tol=1e-7, epsilon=epsilonarr[0], fraction=fractionarr[0], 
            smart_guess=Guess, 
            max_iter=maxiterarr[0],
            )
    

else:
    print("hjb_post_damage_post_tech")
    model_args = (delta, alpha, kappa, mu_k, sigma_k, theta_ell, pi_c_o, sigma_y, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g,rho, gamma_1, gamma_2, gamma_3_i, y_bar, theta, lambda_bar_second, vartheta_bar_second)

    model_tech2_post_damage = hjb_post_tech(
            state_grid=(K, Y, L), 
            model_args=(delta, alpha, theta, vartheta_bar, lambda_bar, mu_k, kappa, sigma_k, theta_ell, pi_c_o, pi_c, sigma_y, zeta, psi_0, psi_1, sigma_g, gamma_1, gamma_2, gamma_3_i, y_bar, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g,rho, varrho),
            V_post_damage=None,
            tol=1e-7, epsilon=epsilonarr[0], fraction=fractionarr[0], 
            smart_guess=Guess, 
            max_iter=maxiterarr[0],
            )
    

# model_tech2_post_damage_check = hjb_post_damage_post_tech(
#         K, Y, model_args, v0=v,
#        epsilon=epsilonarr[0], fraction=fractionarr[0] ,tol=1e-8, max_iter=maxiterarr[0], print_iteration=False)


V_post_3D = model_tech2_post_damage["v0"]

# dvdy = finiteDiff_3D(v_post,1,1,hY)

# print("check post damage post tech function shape:dvdY=[{},{}]".format(dvdy.min(), dvdy.max()))

# V_post_3D = np.zeros_like(K_mat)
# for j in range(nL):
#     V_post_3D[:,:,j] = v_post


print(Data_Dir + File_Name + "model_tech2_post_damage_gamma_{:.4f}".format(gamma_3_i))

print(model_tech2_post_damage.keys())



print("-------------------------------------------")
print("------------Post damage, Tech I-----------")
print("-------------------------------------------")

V_post_tech2 = V_post_3D

with open(Data_Dir+ File_Name + "model_tech1_post_damage_gamma_{:.4f}".format(gamma_3_i_sub), "rb") as f:
    Guess = pickle.load(f)

# Guess = None

if rho==1:
    print("hjb_pre_tech_rho1")

    res = hjb_pre_tech_rho1(
            state_grid=(K, Y, L), 
            model_args=(delta, alpha, theta, vartheta_bar, lambda_bar, mu_k, kappa, sigma_k, theta_ell, pi_c_o, pi_c, sigma_y, zeta, psi_0, psi_1, sigma_g, V_post_tech2, gamma_1, gamma_2, gamma_3_i, y_bar, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g,rho, varrho),
            V_post_damage=None,
            tol=1e-7, epsilon=epsilonarr[1], fraction=fractionarr[1], 
            smart_guess=Guess, 
            max_iter=maxiterarr[1],
            )
else:
    res = hjb_pre_tech(
            state_grid=(K, Y, L), 
            model_args=(delta, alpha, theta, vartheta_bar, lambda_bar, mu_k, kappa, sigma_k, theta_ell, pi_c_o, pi_c, sigma_y, zeta, psi_0, psi_1, sigma_g, V_post_tech2, gamma_1, gamma_2, gamma_3_i, y_bar, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g,rho, varrho),
            V_post_damage=None,
            tol=1e-7, epsilon=epsilonarr[1], fraction=fractionarr[1], 
            smart_guess=Guess, 
            max_iter=maxiterarr[1],
            )


with open(Data_Dir+ File_Name  + "model_tech1_post_damage_gamma_{:.4f}".format(gamma_3_i), "wb") as f:
    pickle.dump(res, f)


    
print("--------------------------------")
print("-----------------Checking Pre Tech---------------")
print("--------------------------------")

with open(Data_Dir+ File_Name + "model_tech1_post_damage_gamma_{:.4f}".format(gamma_3_i), "rb") as f:
    Guess = pickle.load(f)
    
# res_check = hjb_pre_tech(
#         state_grid=(K, Y, L), 
#         model_args=(delta, alpha, theta, vartheta_bar, lambda_bar, mu_k, kappa, sigma_k, theta_ell, pi_c_o, pi_c, sigma_y, zeta, psi_0, psi_1, sigma_g, V_post_tech2, gamma_1, gamma_2, gamma_3_i, y_bar, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g,rho, varrho),
#         V_post_damage=None,
#         tol=1e-7, epsilon=epsilonarr[1], fraction=fractionarr[1], 
#         smart_guess=Guess, 
#         max_iter=maxiterarr[1],
#         )


if rho==1:
    res = hjb_pre_tech_rho1(
        state_grid=(K, Y, L), 
        model_args=(delta, alpha, theta, vartheta_bar, lambda_bar, mu_k, kappa, sigma_k, theta_ell, pi_c_o, pi_c, sigma_y, zeta, psi_0, psi_1, sigma_g, V_post_tech2, gamma_1, gamma_2, gamma_3_i, y_bar, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g,rho, varrho),
        V_post_damage=None,
        tol=1e-7, epsilon=epsilonarr[1], fraction=fractionarr[1], 
        smart_guess=Guess, 
        max_iter=100,
        )


else:
    res = hjb_pre_tech(
        state_grid=(K, Y, L), 
        model_args=(delta, alpha, theta, vartheta_bar, lambda_bar, mu_k, kappa, sigma_k, theta_ell, pi_c_o, pi_c, sigma_y, zeta, psi_0, psi_1, sigma_g, V_post_tech2, gamma_1, gamma_2, gamma_3_i, y_bar, xi_a, xi_k, xi_c, xi_j, xi_d, xi_g,rho, varrho),
        V_post_damage=None,
        tol=1e-7, epsilon=epsilonarr[1], fraction=fractionarr[1], 
        smart_guess=Guess, 
        max_iter=100,
        )


