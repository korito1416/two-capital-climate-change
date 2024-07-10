import numpy as np
import pandas as pd
import sys
print(sys.path)

sys.path.append('./src')

import pickle
import plotly.graph_objects as go
import plotly.offline as pyo
import matplotlib.pyplot as plt
import SolveLinSys
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import CubicSpline
from matplotlib.backends.backend_pdf import PdfPages
# from src.supportfunctions import finiteDiff_3D
import os
import argparse
import time
import petsc4py
from petsc4py import PETSc
import petsclinearsystem
# from Result_support import *
sys.stdout.flush()


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
parser.add_argument("--rhoarr", type=float)
parser.add_argument("--delta", type=float)

parser.add_argument("--psi0arr",nargs='+',type=float)
parser.add_argument("--psi1arr",nargs='+',type=float)
# parser.add_argument("--psi2arr",nargs='+',type=float)
parser.add_argument("--num_gamma",type=int)

parser.add_argument("--hXarr",nargs='+',type=float)
parser.add_argument("--Xminarr",nargs='+',type=float)
parser.add_argument("--Xmaxarr",nargs='+',type=float)

parser.add_argument("--auto",type=int)
parser.add_argument("--IntPeriod",type=int)
parser.add_argument("--plot_year_gamma",type=int)
parser.add_argument("--plot_year_theta",type=int)

parser.add_argument("--scheme",type=str)
parser.add_argument("--HJB_solution",type=str)

# parser.add_argument("--Update",type=int)


args = parser.parse_args()

dataname = args.dataname

# Update = args.Update
IntPeriod = args.IntPeriod
timespan = 1/12
plot_year_gamma=args.plot_year_gamma
plot_year_theta=args.plot_year_theta

psi0arr = args.psi0arr
psi1arr = args.psi1arr
# psi2arr = args.psi2arr
xiaarr = args.xiaarr
xikarr = args.xikarr 
xicarr = args.xicarr 
xijarr = args.xijarr 
xidarr = args.xidarr 
xigarr = args.xigarr 
varrhoarr = args.varrhoarr
rho = args.rhoarr
phi_0 = args.phi_0

# if len(xicarr)==4:
#     labellist = ['Capital Aversion', 'Climate Aversion', 'Technology Aversion', 'Damage Aversion']
#     Filename = 'Uncertainty Channels'
#     colors = ['blue','red', 'green', 'cyan', 'purple']
    
# if len(xicarr)==4 and min(xikarr)==0.050 and phi_0==0.5:
#     labellist = ['Climate Uncertainty', 'Damage Uncertainty', 'Productivity Uncertainty', ' No Technology Uncertainty']
#     Filename = 'Uncertainty Channels Less High Cost'
#     colors = ['blue','red', 'green', 'cyan', 'purple']
        
    
# if len(xicarr)==4 and min(xikarr)==0.025 and phi_0==0.5:
#     labellist = ['Climate Uncertainty', 'Damage Uncertainty', 'Productivity Uncertainty', 'No Technology Uncertainty']
#     Filename = 'Uncertainty Channels More High Cost'
#     colors = ['blue','red', 'green', 'cyan', 'purple']
        
        
        
# if len(xicarr)==4 and min(xikarr)==0.050:
#     labellist = ['Climate Uncertainty', 'Damage Uncertainty', 'Productivity Uncertainty', 'Technology Uncertainty']
#     Filename = 'Uncertainty Channels'
#     colors = ['blue','red', 'green', 'cyan', 'purple']
    

# if len(xicarr)==4 and min(xikarr)==0.075:
#     labellist = ['Climate Uncertainty', 'Damage Uncertainty', 'Productivity Uncertainty', 'Technology Uncertainty']
#     Filename = 'Uncertainty Channels More'
#     colors = ['blue','red', 'green', 'cyan', 'purple']
    
if len(xicarr)==4 and min(xikarr)==0.150:
    labellist = ['Climate Uncertainty', 'Damage Uncertainty', 'Productivity Uncertainty', 'Technology Uncertainty']
    labellist2 = ['Climate Uncertainty', 'Damage Uncertainty', 'Productivity Uncertainty', 'Technology Uncertainty']
    Filename = 'Uncertainty Channels Less'
    colors = ['blue','red', 'green', 'cyan', 'purple']

    
if len(xicarr)==2 and min(xikarr)==0.150:
    labellist = ['Productivity Uncertainty', 'Neutrality']
    Filename = 'DoubleComparison'
    colors = ['blue','red', 'green', 'cyan', 'purple']
    
    
# if len(xicarr)==3 and min(xikarr)==0.150:
#     labellist = ['Productivity Uncertainty', 'Technology Uncertainty', 'Neutrality']
#     Filename = 'TripleComparison'
#     colors = ['blue','green', 'red', 'cyan', 'purple']
    
# # if len(xicarr)==5:
# #     labellist = ['Capital Aversion', 'Climate Aversion', 'Technology Aversion', 'Damage Aversion', 'Full Aversion']
# #     Filename = 'Uncertainty Channels'
# #     # if rho==0.66:
# #     #     Filename = 'Uncertainty Channels_Rho<1'
# #     # if rho==1:
# #     #     Filename = 'Uncertainty Channels_Rho=1'    
# #     # if rho==1.5:
# #     #     Filename = 'Uncertainty Channels_Rho>1'

# #     colors =['blue','red', 'green', 'cyan', 'purple']
    
# if len(xicarr)==3 and min(xikarr)==0.025:
#     labellist = ['More Aversion', 'Less Aversion', 'Neutrality']
#     Filename = 'Aversion Intensity Intense'
#     # Filename = 'Aversion Intensity_old'
#     # Filename = 'Aversion Intensity_onlyj'
#     # Filename = 'Aversion Intensity_onlyk'
#     colors = ['blue','red', 'green', 'cyan', 'purple']
#     colors2 = ['blue','red', 'green', 'cyan', 'purple']


# if len(xicarr)==3:
#     labellist = ['More Aversion', 'Less Aversion', 'Neutrality']
#     # Filename = 'Uncertainty Channel Intense'
#     Filename = 'Aversion Intensity'
#     # Filename = 'Aversion Intensity_old'
#     # Filename = 'Aversion Intensity_onlyj'
#     # Filename = 'Aversion Intensity_onlyk'
#     colors = ['blue','red', 'green', 'cyan', 'purple', 'orange']
#     colors2 = ['blue','red', 'green', 'cyan', 'purple', 'orange']

# if len(xicarr)==3:
#     labellist = ['More Aversion', 'Less Aversion', 'Baseline']
#     # Filename = 'Uncertainty Channel Intense'
#     Filename = 'Aversion Intensity Test'
#     # Filename = 'Aversion Intensity_old'
#     # Filename = 'Aversion Intensity_onlyj'
#     # Filename = 'Aversion Intensity_onlyk'
#     colors = ['blue','red', 'green', 'cyan', 'purple', 'orange']
#     colors2 = ['blue','red', 'green', 'cyan', 'purple', 'orange']

# if len(xicarr)==7:
#     labellist = ['More Aversion', '2', '3', '4', '5', '6', '7']
#     # Filename = 'Uncertainty Channel Intense'
#     Filename = 'Aversion Intensity Test'
#     # Filename = 'Aversion Intensity_old'
#     # Filename = 'Aversion Intensity_onlyj'
#     # Filename = 'Aversion Intensity_onlyk'
#     colors = ['blue','red', 'green', 'cyan', 'purple', 'orange', 'black']
#     colors2 = ['blue','red', 'green', 'cyan', 'purple', 'orange', 'black']

# # if len(xicarr)==7:
# #     labellist = ['More Aversion', 'Less Aversion', 'Neutrality', 'Climate Uncertainty', 'Damage Uncertainty', 'Productivity Uncertainty', 'Technology Uncertainty']
# #     # Filename = 'Uncertainty Channel Intense'
# #     Filename = 'Aversion Intensity Computation'
# #     # Filename = 'Aversion Intensity_old'
# #     # Filename = 'Aversion Intensity_onlyj'
# #     # Filename = 'Aversion Intensity_onlyk'
# #     colors = ['blue','red', 'green', 'cyan', 'purple', 'orange', 'black']
# #     colors2 = ['blue','red', 'green', 'cyan', 'purple', 'orange', 'black']


if len(xicarr)==1:
    labellist = ['More Aversion']
    # Filename = 'Uncertainty Channel Intense'
    Filename = 'Aversion Intensity New'
    # Filename = 'Aversion Intensity_old'
    # Filename = 'Aversion Intensity_onlyj'
    # Filename = 'Aversion Intensity_onlyk'
    colors = ['blue','red', 'green', 'cyan', 'purple', 'orange', 'black']
    colors2 = ['blue','red', 'green', 'cyan', 'purple', 'orange', 'black']


# if len(xicarr)==3 and min(xikarr)==0.005:
#     labellist = ['More Aversion', 'Less Aversion', 'Neutrality']
#     Filename = 'Aversion Intensity Extreme'
#     # Filename = 'Aversion Intensity_old'
#     # Filename = 'Aversion Intensity_onlyj'
#     # Filename = 'Aversion Intensity_onlyk'
#     colors = ['blue','red', 'green', 'cyan', 'purple']
#     colors2 = ['blue','red', 'green', 'cyan', 'purple']


# if len(xicarr)==4 and min(xikarr)==0.005:
#     labellist = [r'$\gamma=201$', r'$\gamma=101$', r'$\gamma=41$', 'Neutrality']
#     Filename = 'Aversion Intensity 4'
#     # Filename = 'Aversion Intensity_old'
#     # Filename = 'Aversion Intensity_onlyj'
#     # Filename = 'Aversion Intensity_onlyk'
#     colors = ['blue','red', 'green', 'cyan', 'purple']
#     colors2 = ['blue','red', 'green', 'cyan', 'purple']

# if len(xicarr)==4:
#     labellist = [r'$\gamma=201$', r'$\gamma=41$', r'$\gamma=14$', 'Neutrality']
#     Filename = 'Aversion Intensity Current'
#     # Filename = 'Aversion Intensity_old'
#     # Filename = 'Aversion Intensity_onlyj'
#     # Filename = 'Aversion Intensity_onlyk'
#     colors = ['purple', 'cyan', 'blue', 'green']
#     colors2 = ['purple', 'cyan', 'blue', 'green']
    

# if len(xicarr)==4:
#     labellist = ['Less Aversion', '2', '3', '4']
#     Filename = 'Aversion Intensity Test'
#     # Filename = 'Aversion Intensity_old'
#     # Filename = 'Aversion Intensity_onlyj'
#     # Filename = 'Aversion Intensity_onlyk'
#     colors = ['purple', 'cyan', 'blue', 'green']
#     colors2 = ['purple', 'cyan', 'blue', 'green']
    
# if len(xicarr)==5 and min(xijarr)==0.005:
#     labellist = [r'$\gamma=201$', r'$\gamma=101$', r'$\gamma=41$', 'Neutrality']
#     Filename = 'R&D Channel on 4'
#     # Filename = 'Aversion Intensity_old'
#     # Filename = 'Aversion Intensity_onlyj'
#     # Filename = 'Aversion Intensity_onlyk'
#     colors = ['blue','red', 'green', 'cyan', 'purple']
#     colors2 = ['blue','red', 'green', 'cyan', 'purple']

if len(xicarr)==6 and min(xijarr)==0.005:
    labellist = [r'$\gamma=201$', r'$\gamma=51$', r'$\gamma=26$', r'$\gamma=13.5$', r'$\gamma=9.3$', r'$\gamma=7.6$']
    Filename = 'RD Channel on 6'
    # Filename = 'Aversion Intensity_old'
    # Filename = 'Aversion Intensity_onlyj'
    # Filename = 'Aversion Intensity_onlyk'
    colors = ['blue','red', 'green', 'cyan', 'purple', 'black']
    colors2 = ['blue','red', 'green', 'cyan', 'purple', 'black']


if len(xicarr)==7 and min(xicarr)==0.005:
    labellist = [r'$\gamma=201$', r'$\gamma=101$', r'$\gamma=51$', r'$\gamma=26$', r'$\gamma=13.5$', r'$\gamma=9.3$', r'$\gamma=7.6$']
    labellist2 = ['$gamma=201$', '$gamma=101$', '$gamma=51$', '$gamma=26$', '$gamma=13.5$', '$gamma=9.3$', '$gamma=7.6$']
    Filename = 'Only RD Channel off 7'
    # Filename = 'Aversion Intensity_old'
    # Filename = 'Aversion Intensity_onlyj'
    # Filename = 'Aversion Intensity_onlyk'
    colors = ['blue','red', 'green', 'cyan', 'purple', 'black', 'orange']
    colors2 = ['blue','red', 'green', 'cyan', 'purple', 'black', 'orange']


if len(xicarr)==7 and min(xijarr)==0.005:
    labellist = [r'$\gamma=201$', r'$\gamma=101$', r'$\gamma=51$', r'$\gamma=26$', r'$\gamma=13.5$', r'$\gamma=9.3$', r'$\gamma=7.6$']
    labellist2 = ['$gamma=201$', '$gamma=101$', '$gamma=51$', '$gamma=26$', '$gamma=13.5$', '$gamma=9.3$', '$gamma=7.6$']
    Filename = 'Only RD Channel on 7'
    # Filename = 'Aversion Intensity_old'
    # Filename = 'Aversion Intensity_onlyj'
    # Filename = 'Aversion Intensity_onlyk'
    colors = ['blue','red', 'green', 'cyan', 'purple', 'black', 'orange']
    colors2 = ['blue','red', 'green', 'cyan', 'purple', 'black', 'orange']

# if len(xicarr)==5 and min(xikarr)==0.005:
#     labellist = [r'$\gamma=201$', r'$\gamma=101$', r'$\gamma=41$', r'$\gamma=21$','Neutrality']
#     Filename = 'Aversion Intensity2 5'
#     # Filename = 'Aversion Intensity_old'
#     # Filename = 'Aversion Intensity_onlyj'
#     # Filename = 'Aversion Intensity_onlyk'
#     colors = ['blue','red', 'green', 'cyan', 'purple']
#     colors2 = ['blue','red', 'green', 'cyan', 'purple']


if len(xicarr)==4 and min(xikarr)==0.005:
    labellist = [r'$\xi=0.005$', r'$\xi=0.025$', r'$\xi=0.075$','Neutrality']
    labellist2 = [r'$xi=0.005$', r'$xi=0.025$', r'$xi=0.075$','Neutrality']
    Filename = 'Extreme Case'
    colors = ['purple','cyan', 'blue', 'green', 'purple']
    
# # if len(xicarr)==3 and min(xikarr)==0.075:
# #     labellist = ['Even Less Aversion', 'Much Less Aversion', 'Neutrality']
# #     Filename = 'Aversion Intensity'
# #     # Filename = 'Aversion Intensity_old'
# #     # Filename = 'Aversion Intensity_onlyj'
# #     # Filename = 'Aversion Intensity_onlyk'
# #     colors = ['blue','red', 'green', 'cyan', 'purple']
# #     colors2 = ['blue','red', 'green', 'cyan', 'purple']

if len(xicarr)==3 and min(xikarr)==0.075:
    labellist = ['More Aversion', 'Less Aversion', 'Neutrality']
    labellist2 = ['More Aversion', 'Less Aversion', 'Neutrality']
    Filename = 'Aversion Intensity'
    # Filename = 'Aversion Intensity_old'
    # Filename = 'Aversion Intensity_onlyj'
    # Filename = 'Aversion Intensity_onlyk'
    colors = ['blue','red', 'green', 'cyan', 'purple']
    colors = ['blue','red', 'green', 'cyan', 'purple']
    print("define success")

# if len(xicarr)==5:
#     labellist = ['Climate Uncertainty', 'Damage Uncertainty', 'Productivity Uncertainty', 'Technology Uncertainty', 'No Technology Uncertainty']
#     Filename = 'Uncertainty Channel5'
#     # Filename = 'Aversion Intensity_old'
#     # Filename = 'Aversion Intensity_onlyj'
#     # Filename = 'Aversion Intensity_onlyk'
#     colors = ['blue','red', 'green', 'cyan', 'purple']
#     colors2 = ['blue','red', 'green', 'cyan', 'purple']

# if len(xicarr)==1:
#     labellist = ['No Technology Uncertainty']
#     Filename = 'Uncertainty Channel1'
#     # Filename = 'Aversion Intensity_old'
#     # Filename = 'Aversion Intensity_onlyj'
#     # Filename = 'Aversion Intensity_onlyk'
#     colors = ['blue','red', 'green', 'cyan', 'purple']
#     colors2 = ['blue','red', 'green', 'cyan', 'purple']

    
# colors = ['blue','green', 'red', 'cyan']

Xminarr = args.Xminarr
Xmaxarr = args.Xmaxarr
hXarr = args.hXarr
auto = args.auto

num_gamma = args.num_gamma
gamma_3_list = np.linspace(0,1./3.,num_gamma)

scheme = args.scheme
HJB_solution = args.HJB_solution


delta = args.delta
alpha = 0.115
kappa = 6.667
mu_k  = -0.045
sigma_k = 0.01
beta_f = 1.86/1000
sigma_y = 1.2 * 1.86 / 1000
zeta = 0.0
# psi_0 = 0.00025
# psi_1 = 1/2
sigma_g = 0.0078
gamma_1 = 1.7675 / 1000
gamma_2 = 0.0022 * 2


y_bar = 2.
y_bar_lower = 1.5

# Tech
theta = 3
lambda_bar = 0.1206
# vartheta_bar = 0.0453
# vartheta_bar = 0.05
# vartheta_bar = 0.056
# vartheta_bar = 0.5
vartheta_bar = args.phi_0

lambda_bar_first = lambda_bar / 2.
vartheta_bar_first = vartheta_bar / 2.

lambda_bar_second = 1e-3
vartheta_bar_second = 0.


# print(plt.rcParamsDefault)
# print("Before, figure default size is: ", plt.rcParams["figure.figsize"])
# print("Before, figure default dpi is: ", plt.rcParams["figure.dpi"])
# print("Before, figure default size is: ", plt.rcParams["font.size"])
# print("Before, legend.frameon is: ", plt.rcParams["legend.frameon"])
# print("Before, lines.linewidth is: ", plt.rcParams["lines.linewidth"])

plt.style.use('classic')
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["figure.figsize"] = (10,8)
plt.rcParams["figure.dpi"] = 500
plt.rcParams["font.size"] = 12
plt.rcParams["legend.frameon"] = False
plt.rcParams["lines.linewidth"] = 5
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors) 
    # colors = ['blue','green', 'red', 'cyan', 'purple']

print("After, figure default size is: ", plt.rcParams["savefig.bbox"])
print("After, figure default size is: ", plt.rcParams["figure.figsize"])
print("After, figure default dpi is: ", plt.rcParams["figure.dpi"])
print("After, figure default size is: ", plt.rcParams["font.size"])
print("After, legend.frameon is: ", plt.rcParams["legend.frameon"])
print("After, lines.linewidth is: ", plt.rcParams["lines.linewidth"])


os.makedirs("./figure/"+args.dataname+"/", exist_ok=True)

Plot_Dir = "./figure/"+args.dataname+"/"

def model_simulation_generate(xi_a,xi_k,xi_c,xi_j,xi_d,xi_g,psi_0,psi_1,varrho,rho):

    # Output_Dir = "/scratch/bincheng/"
    Output_Dir = args.outputname
    Data_Dir = Output_Dir+"abatement/data_2tech/"+args.dataname+"/"
    # File_Dir = "xi_a_{}_xi_k_{}_xi_c_{}_xi_j_{}_xi_d_{}_xi_g_{}_psi_0_{}_psi_1_{}_varrho_{}_rho_{}_" .format(xi_a,xi_k,xi_c,xi_j,xi_d,xi_g,psi_0,psi_1,varrho,rho)
    File_Dir = "xi_a_{}_xi_k_{}_xi_c_{}_xi_j_{}_xi_d_{}_xi_g_{}_psi_0_{}_psi_1_{}_varrho_{}_rho_{}_delta_{}_" .format(xi_a, xi_k, xi_c, xi_j, xi_d, xi_g, psi_0,psi_1, varrho, rho, delta)
    

    with open(Data_Dir + File_Dir+"model_tech1_pre_damage"+"_UD_simul_{}".format(IntPeriod)+ scheme + "_" +HJB_solution, "rb") as f:
        res = pickle.load(f)


    
    return res

print("RD")
for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], ((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], ((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                plt.xlabel('Years')
                plt.ylabel('$\%$ of GDP')
                # plt.title("R&D investment as percentage of  GDP")
                # if auto==0:   
                # if vartheta_bar==0.1:
                #     plt.ylim(0,4)
                if vartheta_bar==0.5:
                    plt.ylim(0,6)
                plt.ylim(0,15)
                plt.xlim(0,40)
                # plt.xlim(0,IntPeriod)

                plt.legend(loc='upper left')        
                
                if id_xiag==2:
                    RD_Neu_0 = (((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[0]
                    RD_Neu_1 = (((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[-1]
                if id_xiag==1:
                    RD_Less_0 = (((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[0]
                    RD_Less_1 = (((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[-1]
                if id_xiag==0:
                    RD_More_0 = (((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[0]
                    RD_More_1 = (((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[-1]
                                    
                if id_xiag>2:
                    
                    print(labellist[id_xiag]+": Neutrality Ratio Left = {:.1f}\%".format((((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[0]/RD_Neu_0*100)+": Neutrality Ratio Right = {:.1f}\%".format((((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[-1]/RD_Neu_1*100))
                    print(labellist[id_xiag]+": Less Ratio Left = {:.1f}\%".format((((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[0]/RD_Less_0*100)+": Neutrality Ratio Right = {:.1f}\%".format((((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[-1]/RD_Less_1*100))
                    print(labellist[id_xiag]+": More Ratio Left = {:.1f}\%".format((((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[0]/RD_More_0*100)+": Neutrality Ratio Right = {:.1f}\%".format((((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[-1]/RD_More_1*100))
                print("delta = {}".format(delta)+"xi = {}, number = {:.1f}\%".format(xikarr[id_xiag],(((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[0]))


print(res.keys())
plt.savefig(Plot_Dir+"/RD_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/RD_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

print("RD_log")


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)


                if id_xiag==2:
                    RD_Neu_0 = np.log(((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[0]
                    RD_Neu_1 = np.log(((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[-1]
                if id_xiag==1:
                    RD_Less_0 = np.log(((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[0]
                    RD_Less_1 = np.log(((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[-1]
                    print(RD_Less_0)
                if id_xiag==0:
                    RD_More_0 = np.log(((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[0]
                    RD_More_1 = np.log(((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[-1]

                if id_xiag>2:
                    
                    temp_0 = np.log(((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[0]
                    temp_1 = np.log(((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5])[0-1]
                    print(temp_0,temp_1)
                    print(labellist[id_xiag]+": Neutrality Ratio Left = {:.2f}".format(RD_Neu_0-temp_0)+": Neutrality Ratio Right = {:.2f}".format(RD_Neu_1-temp_1))
                    print(labellist[id_xiag]+": Less Ratio Left = {:.2f}".format(RD_Less_0-temp_0)+": Neutrality Ratio Right = {:.2f}".format(RD_Less_1-temp_1))
                    print(labellist[id_xiag]+": More Ratio Left = {:.2f}".format(RD_More_0-temp_0)+": Neutrality Ratio Right = {:.2f}".format(RD_More_1-temp_1))
                                    


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], ((res["c"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], ((res["c"]/(alpha*np.exp(res["states"][:,0])))*100)[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                plt.xlabel('Years')
                plt.ylabel('$\%$ of GDP')
                # plt.title("R&D investment as percentage of  GDP")
                # if auto==0:   
                # if vartheta_bar==0.1:
                #     plt.ylim(0,4)
                if vartheta_bar==0.5:
                    plt.ylim(20,40)
                # plt.xlim(0,30)
                plt.xlim(0,IntPeriod)

                plt.legend(loc='upper left')        
print(res.keys())
plt.savefig(Plot_Dir+"/ConOutput_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/ConOutput_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["x"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["x"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                plt.xlabel('Years')
                # plt.ylabel('')
                plt.title("Raw R&D investment")
                # if auto==0:   
                # plt.ylim(0,1)
                plt.xlim(0,30)

                plt.legend(loc='upper left')        
print(res.keys())
plt.savefig(Plot_Dir+"/RD_Raw"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/RD_Raw"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()



for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)*(1-res["true_tech_prob"]))[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (((res["x"]/(alpha*np.exp(res["states"][:,0])))*100)*(1-res["true_tech_prob"]))[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                plt.xlabel('Years')
                plt.ylabel('$\%$ of GDP')
                plt.title("R&D investment as percentage of GDP")
                # if auto==0:   
                plt.ylim(0,15)
                plt.xlim(0,40)

                plt.legend(loc='upper left')        
print(res.keys())
plt.savefig(Plot_Dir+"/RD_Expected_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/RD_Expected_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()



for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["i"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["i"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                plt.xlabel('Years')
                # plt.title("Capital investment")
                # if auto==0:   
                plt.ylim(30,160)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

# plt.savefig(Plot_Dir+"/CapI_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/CapI_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()



for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], (res["i"]/(np.exp(res["states"][:,0]))*100)[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (res["i"]/(np.exp(res["states"][:,0]))*100)[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                plt.xlabel('Years')
                # plt.title("Capital investment as percentage of GDP")
                # if auto==0:   
                plt.ylim(5,15)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper left')

# plt.savefig(Plot_Dir+"/CapI_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/CapIRatio_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

print("CapIRatioOutput_")
for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], (res["i"]/(alpha*np.exp(res["states"][:,0]))*100)[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (res["i"]/(alpha*np.exp(res["states"][:,0]))*100)[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                plt.xlabel('Years')
                plt.ylabel('$\%$ of GDP')
                # plt.title("Capital investment as percentage of GDP")
                # if auto==0:   
                plt.ylim(50,80)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper right')


                print("delta = {}".format(delta)+"xi = {}, number = {:.1f}\%".format(xikarr[id_xiag],((res["i"]/(alpha*np.exp(res["states"][:,0]))*100)[res["states"][:, 1]<1.5])[0]))
# plt.savefig(Plot_Dir+"/CapI_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/CapIRatioOutput_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

data = pd.read_csv("./data/emission.csv")

print("Emission")
for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["e"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["e"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                # plt.plot(res2["years"][res2["states"][:, 1]<1.5], res2["e"][res2["states"][:, 1]<1.5],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"][res3["states"][:, 1]<1.5], res3["e"][res3["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                # plt.title("Carbon Emissions")
                # if auto==0:   
                plt.ylim(0.0,20.0)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

                
                if id_xiag==2:
                    RD_Neu_0 = (res["e"][res["states"][:, 1]<1.5])[0]
                    RD_Neu_1 = (res["e"][res["states"][:, 1]<1.5])[-1]
                if id_xiag==1:
                    RD_Less_0 = (res["e"][res["states"][:, 1]<1.5])[0]
                    RD_Less_1 = (res["e"][res["states"][:, 1]<1.5])[-1]
                if id_xiag==0:
                    RD_More_0 = (res["e"][res["states"][:, 1]<1.5])[0]
                    RD_More_1 = (res["e"][res["states"][:, 1]<1.5])[-1]
                                    
                if id_xiag>2:
                    
                    print(labellist[id_xiag]+": Neutrality Ratio Left = {:.1f}\%".format(RD_Neu_0/(res["e"][res["states"][:, 1]<1.5])[0]*100)+": Neutrality Ratio Right = {:.1f}\%".format(RD_Neu_1/(res["e"][res["states"][:, 1]<1.5])[-1]*100))
                    print(labellist[id_xiag]+": Less Ratio Left = {:.1f}\%".format(RD_Less_0/(res["e"][res["states"][:, 1]<1.5])[0]*100)+": Neutrality Ratio Right = {:.1f}\%".format(RD_Less_1/(res["e"][res["states"][:, 1]<1.5])[-1]*100))
                    print(labellist[id_xiag]+": More Ratio Left = {:.1f}\%".format(RD_More_0/(res["e"][res["states"][:, 1]<1.5])[0]*100)+": Neutrality Ratio Right = {:.1f}\%".format(RD_More_1/(res["e"][res["states"][:, 1]<1.5])[-1]*100))
                                    


plt.savefig(Plot_Dir+"/E_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/E_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

print("Emission_Log")

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (res["e"]*(1-res["true_tech_prob"]))[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (res["e"]*(1-res["true_tech_prob"]))[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                # plt.plot(res2["years"][res2["states"][:, 1]<1.5], res2["e"][res2["states"][:, 1]<1.5],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"][res3["states"][:, 1]<1.5], res3["e"][res3["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                # plt.title("Carbon Emissions")
                # if auto==0:   
                plt.ylim(0.0,20.0)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

                
                if id_xiag==2:
                    RD_Neu_0 = (res["e"][res["states"][:, 1]<1.5])[0]
                    RD_Neu_1 = (res["e"][res["states"][:, 1]<1.5])[-1]
                if id_xiag==1:
                    RD_Less_0 = (res["e"][res["states"][:, 1]<1.5])[0]
                    RD_Less_1 = (res["e"][res["states"][:, 1]<1.5])[-1]
                if id_xiag==0:
                    RD_More_0 = (res["e"][res["states"][:, 1]<1.5])[0]
                    RD_More_1 = (res["e"][res["states"][:, 1]<1.5])[-1]
                                    
                if id_xiag>2:
                    
                    print(labellist[id_xiag]+": Neutrality Ratio Left = {:.1f}\%".format(RD_Neu_0/(res["e"][res["states"][:, 1]<1.5])[0]*100)+": Neutrality Ratio Right = {:.1f}\%".format(RD_Neu_1/(res["e"][res["states"][:, 1]<1.5])[-1]*100))
                    print(labellist[id_xiag]+": Less Ratio Left = {:.1f}\%".format(RD_Less_0/(res["e"][res["states"][:, 1]<1.5])[0]*100)+": Neutrality Ratio Right = {:.1f}\%".format(RD_Less_1/(res["e"][res["states"][:, 1]<1.5])[-1]*100))
                    print(labellist[id_xiag]+": More Ratio Left = {:.1f}\%".format(RD_More_0/(res["e"][res["states"][:, 1]<1.5])[0]*100)+": Neutrality Ratio Right = {:.1f}\%".format(RD_More_1/(res["e"][res["states"][:, 1]<1.5])[-1]*100))
                                    


plt.savefig(Plot_Dir+"/E_Expected"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/E_Expected"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()



for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                # if xiaarr[id_xiag]>10:
                #     plt.plot(res["years"][res["states"][:, 1]<1.5], res["e"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                # else:
                #     plt.plot(res["years"][res["states"][:, 1]<1.5], res["e"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                # # plt.plot(res2["years"][res2["states"][:, 1]<1.5], res2["e"][res2["states"][:, 1]<1.5],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # # plt.plot(res3["years"][res3["states"][:, 1]<1.5], res3["e"][res3["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=7.0)
                # plt.xlabel('Years')
                # # plt.title("Carbon Emissions")
                # # if auto==0:   
                # plt.ylim(3.0,15.0)
                # plt.xlim(0,30)
                # plt.legend(loc='upper left')


                if id_xiag==2:
                    RD_Neu_0 = np.log(res["e"][res["states"][:, 1]<1.5])[0]
                    RD_Neu_1 = np.log(res["e"][res["states"][:, 1]<1.5])[-1]
                if id_xiag==1:
                    RD_Less_0 = np.log(res["e"][res["states"][:, 1]<1.5])[0]
                    RD_Less_1 = np.log(res["e"][res["states"][:, 1]<1.5])[-1]
                if id_xiag==0:
                    RD_More_0 = np.log(res["e"][res["states"][:, 1]<1.5])[0]
                    RD_More_1 = np.log(res["e"][res["states"][:, 1]<1.5])[-1]
                                    
                if id_xiag>2:
                    
                    print(labellist[id_xiag]+": Neutrality Ratio Left = {:.2f}".format(RD_Neu_0-np.log(res["e"][res["states"][:, 1]<1.5])[0])+": Neutrality Ratio Right = {:.2f}".format(RD_Neu_1-np.log(res["e"][res["states"][:, 1]<1.5])[-1]))
                    print(labellist[id_xiag]+": Less Ratio Left = {:.2f}".format(RD_Less_0-np.log(res["e"][res["states"][:, 1]<1.5])[0])+": Neutrality Ratio Right = {:.2f}".format(RD_Less_1-np.log(res["e"][res["states"][:, 1]<1.5])[-1]))
                    print(labellist[id_xiag]+": More Ratio Left = {:.2f}".format(RD_More_0-np.log(res["e"][res["states"][:, 1]<1.5])[0])+": Neutrality Ratio Right = {:.2f}".format(RD_More_1-np.log(res["e"][res["states"][:, 1]<1.5])[-1]))
                                    

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:
                    plt.plot(res["years"], res["e"],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"], res["e"],label=labellist[id_xiag],linewidth=5.0)
                # plt.plot(res2["years"][res2["states"][:, 1]<1.5], res2["e"][res2["states"][:, 1]<1.5],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"][res3["states"][:, 1]<1.5], res3["e"][res3["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                # plt.title("Carbon Emissions")
                # if auto==0:   
                # plt.ylim(6.0,15.0)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper left')


txt="""Annual emissions, initialized at the year 2020, in units of GtC. \n The three solid lines correspond to trajectories of annual emissions under three model misspecification scenarios: \n neutrality, less aversion and more aversion from Barnett, Brock, Hansen and Zhang (2023). \n The four dashed lines correspond to trajectories of emissions under four projected RCP scenarios: \n RCP 2.6, RCP 4.5, RCP 6.0 and RCP 8.5 from the 2017 Climate Science Special Report by the U.S. Global Change Research Program"""
plt.figtext(0.5, -0.2, txt, wrap=True, horizontalalignment='center', fontsize=18)

plt.plot(data.loc[13:17,"YEARS"]-2020,data.loc[13:17,"RCP3"],label="RCP2.6",linestyle="--", color="palegreen")
plt.plot(data.loc[13:17,"YEARS"]-2020,data.loc[13:17,"RCP4.5"],label="RCP4.5",linestyle="--", color="lightskyblue")
plt.plot(data.loc[13:17,"YEARS"]-2020,data.loc[13:17,"RCP6"],label="RCP6.0",linestyle="--", color = "orange")
plt.plot(data.loc[13:17,"YEARS"]-2020,data.loc[13:17,"RCP8.5"],label="RCP8.5",linestyle="--", color="salmon")
plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/E2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/E2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:
                    plt.plot(res["years"], res["e"]* (1-res["true_tech_prob"]),label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"], res["e"]* (1-res["true_tech_prob"]),label=labellist[id_xiag],linewidth=5.0)
                # plt.plot(res2["years"][res2["states"][:, 1]<1.5], res2["e"][res2["states"][:, 1]<1.5],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"][res3["states"][:, 1]<1.5], res3["e"][res3["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                # plt.title("Carbon Emissions")
                # if auto==0:   
                # plt.ylim(6.0,15.0)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper left')


txt="""Annual emissions, initialized at the year 2020, in units of GtC. \n The three solid lines correspond to trajectories of annual expected emissions under three model misspecification scenarios: \n neutrality, less aversion and more aversion from Barnett, Brock, Hansen and Zhang (2023). \n The four dashed lines correspond to trajectories of emissions under four projected RCP scenarios: \n RCP 2.6, RCP 4.5, RCP 6.0 and RCP 8.5 from the 2017 Climate Science Special Report by the U.S. Global Change Research Program"""
plt.figtext(0.5, -0.2, txt, wrap=True, horizontalalignment='center', fontsize=18)

plt.plot(data.loc[13:17,"YEARS"]-2020,data.loc[13:17,"RCP3"],label="RCP2.6",linestyle="--", color="palegreen")
plt.plot(data.loc[13:17,"YEARS"]-2020,data.loc[13:17,"RCP4.5"],label="RCP4.5",linestyle="--", color="lightskyblue")
plt.plot(data.loc[13:17,"YEARS"]-2020,data.loc[13:17,"RCP6"],label="RCP6.0",linestyle="--", color = "orange")
plt.plot(data.loc[13:17,"YEARS"]-2020,data.loc[13:17,"RCP8.5"],label="RCP8.5",linestyle="--", color="salmon")
plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/EExpec2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/EExpec2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], 886+np.cumsum(res["e"][res["states"][:, 1]<1.5])*timespan,label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], 886+np.cumsum(res["e"][res["states"][:, 1]<1.5])*timespan,label=labellist[id_xiag],linewidth=5.0)
                # plt.plot(res2["years"][res2["states"][:, 1]<1.5], res2["e"][res2["states"][:, 1]<1.5],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"][res3["states"][:, 1]<1.5], res3["e"][res3["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                # plt.title("Carbon Emissions")
                # if auto==0:   
                # plt.ylim(6.0,15.0)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/ECum_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/ECum_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], ( res["e"] * (1-res["true_tech_prob"]) )[res["states"][:, 1]<1.5] ,label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], ( res["e"] * (1-res["true_tech_prob"]) )[res["states"][:, 1]<1.5] ,label=labellist[id_xiag],linewidth=5.0)
                # plt.plot(res2["years"][res2["states"][:, 1]<1.5], res2["e"][res2["states"][:, 1]<1.5],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"][res3["states"][:, 1]<1.5], res3["e"][res3["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                plt.title("Carbon Emissions")
                # if auto==0:   
                plt.ylim(6.0,12.0)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/ETrue_Expected_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/ETrue_Expected_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["states"][:, 1][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["states"][:, 1][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                # plt.plot(res2["years"][res2["states"][:, 1]<1.5], res2["states"][:, 1][res2["states"][:, 1]<1.5],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"][res3["states"][:, 1]<1.5], res3["states"][:, 1][res3["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                # plt.title("Temperature Anomaly")
                # if auto==0:   
                plt.ylim(1.1,1.5)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/TA_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/TA_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()


# data_cmip5= pd.read_csv("./data/cmip5.csv")
data_back= pd.read_csv("./data/BackData.csv")


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"], res["states"][:, 1],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"], res["states"][:, 1],label=labellist[id_xiag],linewidth=5.0)
                # plt.plot(res2["years"][res2["states"][:, 1]<1.5], res2["states"][:, 1][res2["states"][:, 1]<1.5],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"][res3["states"][:, 1]<1.5], res3["states"][:, 1][res3["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                # plt.title("Temperature Anomaly")
                # if auto==0:   
                plt.ylim(1,3)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper left')

plt.plot(data_back.loc[0:8,"YEARS"]-2020,data_back.loc[0:8,"RCP3"]*5/9,label="RCP2.6",linestyle="--")
plt.plot(data_back.loc[0:8,"YEARS"]-2020,data_back.loc[0:8,"RCP4.5"]*5/9,label="RCP4.5",linestyle="--")
plt.plot(data_back.loc[0:8,"YEARS"]-2020,data_back.loc[0:8,"RCP8"]*5/9,label="RCP8.5",linestyle="--")
plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/TAF_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/TAF_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()




for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"], res["states"][:, 1],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"], res["states"][:, 1],label=labellist[id_xiag],linewidth=5.0)
                # plt.plot(res2["years"][res2["states"][:, 1]<1.5], res2["states"][:, 1][res2["states"][:, 1]<1.5],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"][res3["states"][:, 1]<1.5], res3["states"][:, 1][res3["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                plt.title("Temperature Anomaly")
                # if auto==0:   
                plt.ylim(1,3)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper left')

plt.plot(data_back.loc[0:8,"YEARS"]-2020,data_back.loc[0:8,"RCP3"]*5/9,label="RCP2.6",linestyle="--", color="palegreen")
plt.plot(data_back.loc[0:8,"YEARS"]-2020,data_back.loc[0:8,"RCP4.5"]*5/9,label="RCP4.5",linestyle="--", color = "lightskyblue")
plt.plot(data_back.loc[0:8,"YEARS"]-2020,data_back.loc[0:8,"RCP8"]*5/9,label="RCP8.5",linestyle="--", color = "salmon")
plt.legend(loc='upper left')

txt="""Temperature anomaly, initialized at the year 2020, in units of degrees Celsius. \n The three solid lines correspond to trajectories of temperature anomaly under three model misspecification scenarios: \n neutrality, less aversion and more aversion from Barnett, Brock, Hansen and Zhang (2023). \n The three dashed lines correspond to trajectories of temperature anomaly under three projected RCP scenarios: \n RCP 2.6, RCP 4.5 and RCP 8.5 from the 2017 Climate Science Special Report by the U.S. Global Change Research Program"""
plt.figtext(0.5, -0.2, txt, wrap=True, horizontalalignment='center', fontsize=18)

plt.savefig(Plot_Dir+"/TAF2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/TAF2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()


    
for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                e_annual = np.zeros([IntPeriod+1])
                year_annual = np.zeros([IntPeriod+1])


                for t in range(IntPeriod):
    
                    e_annual[t] = (res["e"]* (1-res["true_tech_prob"]))[t*12]
                    year_annual[t] = t 
                    
                e_annual [-1] = (res["e"]* (1-res["true_tech_prob"]))[-1]
                year_annual[-1] = IntPeriod

    
                if xiaarr[id_xiag]>10:
                    plt.plot(year_annual, np.cumsum(e_annual)*0.00186+1.1,label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(year_annual, np.cumsum(e_annual)*0.00186+1.1,label=labellist[id_xiag],linewidth=5.0)
                # plt.plot(res2["years"][res2["states"][:, 1]<1.5], res2["e"][res2["states"][:, 1]<1.5],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"][res3["states"][:, 1]<1.5], res3["e"][res3["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                # plt.title("Carbon Emissions")
                # if auto==0:   
                # plt.ylim(6.0,15.0)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper left')


txt="""Temperature anomaly, initialized at the year 2020, in units of degrees Celsius. \n The three solid lines correspond to trajectories of temperature anomaly with expectation value under three model misspecification scenarios: \n neutrality, less aversion and more aversion from Barnett, Brock, Hansen and Zhang (2023). \n The three dashed lines correspond to trajectories of temperature anomaly under three projected RCP scenarios: \n RCP 2.6, RCP 4.5 and RCP 8.5 from the 2017 Climate Science Special Report by the U.S. Global Change Research Program"""
plt.figtext(0.5, -0.2, txt, wrap=True, horizontalalignment='center', fontsize=18)

plt.plot(data_back.loc[0:8,"YEARS"]-2020,data_back.loc[0:8,"RCP3"]*5/9,label="RCP2.6",linestyle="--", color="palegreen")
plt.plot(data_back.loc[0:8,"YEARS"]-2020,data_back.loc[0:8,"RCP4.5"]*5/9,label="RCP4.5",linestyle="--", color = "lightskyblue")
plt.plot(data_back.loc[0:8,"YEARS"]-2020,data_back.loc[0:8,"RCP8"]*5/9,label="RCP8.5",linestyle="--", color = "salmon")
plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/TA_TCREExpec2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/TA_TCREExpec2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()


plt.plot(data.loc[12:17,"YEARS"]-2010,(np.cumsum(data.loc[:,"RCP3"])*.00186*10)[12:18],label="RCP3-TCRE",color="blue")
plt.plot(data.loc[12:17,"YEARS"]-2010,(np.cumsum(data.loc[:,"RCP4.5"])*.00186*10)[12:18],label="RCP4.5-TCRE",color="red")
# plt.plot(data.loc[13:17,"YEARS"],(np.cumsum(data.loc[:,"RCP6"])*.00186*10)[13:18],label="RCP6")
plt.plot(data.loc[12:17,"YEARS"]-2010,(np.cumsum(data.loc[:,"RCP8.5"])*.00186*10)[12:18],label="RCP8.5-TCRE",color="green")
plt.plot(data_back.loc[0:8,"YEARS"]-2020,data_back.loc[0:8,"RCP3"]*5/9,label="RCP3",linestyle="--",color="blue")
plt.plot(data_back.loc[0:8,"YEARS"]-2020,data_back.loc[0:8,"RCP4.5"]*5/9,label="RCP4.5",linestyle="--",color="red")
plt.plot(data_back.loc[0:8,"YEARS"]-2020,data_back.loc[0:8,"RCP8"]*5/9,label="RCP8.5",linestyle="--",color="green")
plt.legend(loc='upper left')
plt.xlim(0,IntPeriod)
plt.ylim(1,3)
plt.savefig(Plot_Dir+"/TACompare_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.exp(res["states"][:, 2])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.exp(res["states"][:, 2])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                # plt.plot(res2["years"][res2["states"][:, 1]<1.5], np.exp(res2["states"][:, 2])[res2["states"][:, 1]<1.5],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"][res3["states"][:, 1]<1.5], np.exp(res3["states"][:, 2])[res3["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                # plt.title("Knowledge Stock $J_g$")
                # if auto==0:   
                plt.ylim(11.0,45.0)
                plt.xlim(0,30)
                plt.legend(loc='upper left')


plt.savefig(Plot_Dir+"/Ig_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/Ig_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()



for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"], res["distorted_tech_prob"],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"], res["distorted_tech_prob"],label=labellist[id_xiag],linewidth=5.0)
                # plt.plot(res2["years"], res2["distorted_tech_prob"],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"], res3["distorted_tech_prob"],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                # plt.title("Distorted Probability of a Technology Jump")
                plt.ylim(0,1)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/ProbTechJump_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/ProbTechJump_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"], res["distorted_damage_prob"],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"], res["distorted_damage_prob"],label=labellist[id_xiag],linewidth=5.0)
                # plt.plot(res2["years"], res2["distorted_damage_prob"],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"], res3["distorted_damage_prob"],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                # plt.title("Distorted Probability of Damage Changes")
                plt.ylim(0,1)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/ProbDamageChange_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/ProbDamageChange_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                color_one = colors[id_xiag % len(xiaarr)]   

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"], res["true_tech_prob"],label=labellist[id_xiag],linewidth=5.0,linestyle = 'dashed',color=color_one)
                else:
                    plt.plot(res["years"], res["true_tech_prob"],label=labellist[id_xiag] ,linewidth=5.0,linestyle = 'dashed',color=color_one)
                    
                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"], res["distorted_tech_prob"],label=labellist[id_xiag],linewidth=5.0,color=color_one)
                else:
                    plt.plot(res["years"], res["distorted_tech_prob"],label=labellist[id_xiag] ,linewidth=5.0,color=color_one)
                # plt.plot(res2["years"], res2["distorted_tech_prob"],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"], res3["distorted_tech_prob"],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                plt.title("Distorted(Solid) and True(Dashed) Probability of a Technology Jump")
                plt.ylim(0,1)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/CombProbTechJump_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/CombProbTechJump_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                color_one = colors[id_xiag % len(xiaarr)]   

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"], res["true_damage_prob"],label=labellist[id_xiag],linewidth=5.0,linestyle = 'dashed',color=color_one)
                else:
                    plt.plot(res["years"], res["true_damage_prob"],label=labellist[id_xiag] ,linewidth=5.0,linestyle = 'dashed',color=color_one)
                    
                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"], res["distorted_damage_prob"],label=labellist[id_xiag],linewidth=5.0,color=color_one)
                else:
                    plt.plot(res["years"], res["distorted_damage_prob"],label=labellist[id_xiag] ,linewidth=5.0,color=color_one)
                # plt.plot(res2["years"], res2["distorted_tech_prob"],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"], res3["distorted_tech_prob"],label=labellist[id_xiag],linewidth=7.0)
                # plt.plot(res2["years"], res2["distorted_damage_prob"],label=r'$\xi_a=\\xi_g=0.050$',linewidth=7.0)
                # plt.plot(res3["years"], res3["distorted_damage_prob"],label=labellist[id_xiag],linewidth=7.0)
                plt.xlabel('Years')
                plt.title("Distorted Probability of Damage Changes")
                plt.ylim(0,1)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/CombProbDamageChange_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/CombProbDamageChange_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"], res["true_tech_prob"],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"], res["true_tech_prob"],label=labellist[id_xiag],linewidth=5.0)

                plt.xlabel("Years")
                plt.title("True Probability of a Technology Jump")
                plt.ylim(0.0,1.0)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/TPIg_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/TPIg_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"], res["true_damage_prob"],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"], res["true_damage_prob"],label=labellist[id_xiag],linewidth=5.0)

                plt.xlabel("Years")
                plt.title("True Probability of Damage Changes")
                plt.ylim(0,1)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/TPId_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/TPId_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.log(res["scc"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.log(res["scc"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)

                plt.xlabel("Years")
                plt.title("Log of Social Cost of Carbon")
                # if auto==0:   
                plt.ylim(4.0,8.0)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/logSCC_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/logSCC_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

print("SVRD")
for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.log(res["scrd"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.log(res["scrd"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)

                plt.xlabel("Years")
                plt.ticklabel_format(useOffset=False)

                # plt.title("Log of Social Value of R&D")
                # if auto==0:   
                plt.ylim(3.0,12.0)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

                if id_xiag==2:
                    RD_Neu_0 = (np.log(res["scrd"])[res["states"][:, 1]<1.5])[0]
                    RD_Neu_1 = (np.log(res["scrd"])[res["states"][:, 1]<1.5])[-1]
                if id_xiag==1:
                    RD_Less_0 = (np.log(res["scrd"])[res["states"][:, 1]<1.5])[0]
                    RD_Less_1 = (np.log(res["scrd"])[res["states"][:, 1]<1.5])[-1]
                if id_xiag==0:
                    RD_More_0 = (np.log(res["scrd"])[res["states"][:, 1]<1.5])[0]
                    RD_More_1 = (np.log(res["scrd"])[res["states"][:, 1]<1.5])[-1]
                                    
                if id_xiag>2:
                    
                    print(labellist[id_xiag]+": Neutrality Ratio Left = {:.2f}".format(RD_Neu_0-(np.log(res["scrd"])[res["states"][:, 1]<1.5])[0])+": Neutrality Ratio Right = {:.2f}".format(RD_Neu_1-(np.log(res["scrd"])[res["states"][:, 1]<1.5])[-1]))
                    print(labellist[id_xiag]+": Less Ratio Left = {:.2f}".format(RD_Less_0-(np.log(res["scrd"])[res["states"][:, 1]<1.5])[0])+": Neutrality Ratio Right = {:.2f}".format(RD_Less_1-(np.log(res["scrd"])[res["states"][:, 1]<1.5])[-1]))
                    print(labellist[id_xiag]+": More Ratio Left = {:.2f}".format(RD_More_0-(np.log(res["scrd"])[res["states"][:, 1]<1.5])[0])+": Neutrality Ratio Right = {:.2f}".format(RD_More_1-(np.log(res["scrd"])[res["states"][:, 1]<1.5])[-1]))
                                    


plt.savefig(Plot_Dir+"/logSVRD_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/logSVRD_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

print("SCGW")


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.log(res["scgw"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.log(res["scgw"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)

                plt.xlabel("Years")
                plt.ticklabel_format(useOffset=False)

                # plt.title("Log of Social Cost of Global Warming")
                # if auto==0:   
                plt.ylim(11.0,14.0)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

                if id_xiag==2:
                    RD_Neu_0 = (np.log(res["scgw"])[res["states"][:, 1]<1.5])[0]
                    RD_Neu_1 = (np.log(res["scgw"])[res["states"][:, 1]<1.5])[-1]
                if id_xiag==1:
                    RD_Less_0 = (np.log(res["scgw"])[res["states"][:, 1]<1.5])[0]
                    RD_Less_1 = (np.log(res["scgw"])[res["states"][:, 1]<1.5])[-1]
                if id_xiag==0:
                    RD_More_0 = (np.log(res["scgw"])[res["states"][:, 1]<1.5])[0]
                    RD_More_1 = (np.log(res["scgw"])[res["states"][:, 1]<1.5])[-1]
                                    
                if id_xiag>2:
                    
                    print(labellist[id_xiag]+": Neutrality Ratio Left = {:.2f}".format(RD_Neu_0-(np.log(res["scgw"])[res["states"][:, 1]<1.5])[0])+": Neutrality Ratio Right = {:.2f}".format(RD_Neu_1-(np.log(res["scgw"])[res["states"][:, 1]<1.5])[-1]))
                    print(labellist[id_xiag]+": Less Ratio Left = {:.2f}".format(RD_Less_0-(np.log(res["scgw"])[res["states"][:, 1]<1.5])[0])+": Neutrality Ratio Right = {:.2f}".format(RD_Less_1-(np.log(res["scgw"])[res["states"][:, 1]<1.5])[-1]))
                    print(labellist[id_xiag]+": More Ratio Left = {:.2f}".format(RD_More_0-(np.log(res["scgw"])[res["states"][:, 1]<1.5])[0])+": Neutrality Ratio Right = {:.2f}".format(RD_More_1-(np.log(res["scgw"])[res["states"][:, 1]<1.5])[-1]))
                                    

plt.savefig(Plot_Dir+"/logSCGW_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/logSCGW_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.log(res["scgw2"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.log(res["scgw2"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)

                plt.xlabel("Years")
                plt.ticklabel_format(useOffset=False)

                # plt.title("Log of Social Cost of Global Warming")
                # if auto==0:   
                plt.ylim(11.0,14.0)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/logSCGW2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/logSCGW2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.log(res["scrd_2"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.log(res["scrd_2"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)

                plt.xlabel("Years")
                plt.ticklabel_format(useOffset=False)

                plt.title("Log of Term (ii)")
                # if auto==0:   
                # plt.ylim(4.0,8.0)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/logSVRD2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/logSVRD2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.log(res["spo"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.log(res["spo"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)

                plt.xlabel("Years")
                plt.ticklabel_format(useOffset=False)

                plt.title("Log of Term (iii)")
                # if auto==0:   
                # plt.ylim(4.0,8.0)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/spo_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/spo_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.log(res["spo2"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], np.log(res["spo2"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)

                plt.xlabel("Years")
                plt.ticklabel_format(useOffset=False)

                plt.title("Log of Social Payoff")
                # if auto==0:   
                # plt.ylim(4.0,8.0)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/spo2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/spo2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["uncertainty_adjusted_diff"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["uncertainty_adjusted_diff"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)

                plt.xlabel("Years")
                plt.ticklabel_format(useOffset=False)

                plt.title("Uncertainty Adjusted Difference")
                # if auto==0:   
                plt.ylim(-0.1, 0)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/UAD_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/UAD_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()



for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["uncertainty_adjusted_diff2"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=1.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["uncertainty_adjusted_diff2"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=1.0)

                plt.xlabel("Years")
                plt.ticklabel_format(useOffset=False)

                plt.title("Uncertainty Adjusted Difference")
                # if auto==0:   
                # plt.ylim(-0.1, 0)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/UAD2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/UAD2_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()



for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["v_post_techt"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5],res["v_post_techt"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)

                plt.xlabel("Years")
                plt.ticklabel_format(useOffset=False)

                plt.title("$V_g$")
                # if auto==0:   
                plt.ylim(4.0,5.3)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/Vgt_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/Vgt_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["vt"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], res["vt"][res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)

                plt.xlabel("Years")
                plt.ticklabel_format(useOffset=False)

                plt.title("$V$")
                # if auto==0:   
                plt.ylim(4.0,5.3)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/Vt_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/Vt_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], (res["v_post_techt"]-res["vt"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (res["v_post_techt"]-res["vt"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)

                plt.xlabel("Years")
                plt.ticklabel_format(useOffset=False)

                plt.title("$V_g-V$")
                # if auto==0:   
                plt.ylim(0.0,1)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/Vg-Vt_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/Vg-Vt_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):


                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                if xiaarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], (res["gt_tech"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (res["gt_tech"])[res["states"][:, 1]<1.5],label=labellist[id_xiag],linewidth=5.0)

                plt.xlabel("Years")
                plt.ticklabel_format(useOffset=False)

                plt.title("$g_t$")
                # if auto==0:   
                # plt.ylim(0.0,1)
                plt.xlim(0,30)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/gt_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".pdf")
plt.savefig(Plot_Dir+"/gt_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png")
plt.close()


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], ( sigma_y * res["e"] * res["ht"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], ( sigma_y * res["e"] * res["ht"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                plt.title(r"$\sigma_y e h$")
                plt.ylim(0,0.005)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/h_{}_"+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/h_{}_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], -(sigma_k * res["hkt"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], -(sigma_k * res["hkt"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                plt.title(r"-$\sigma_k h_k$")
                plt.ylim(0,0.005)
                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/hk_{}_"+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/hk_{}_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], -(sigma_g * res["hjt"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], -(sigma_g * res["hjt"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                plt.title(r"-$\sigma_j h_j$")
                plt.ylim(0,0.005)

                plt.legend(loc='upper left')

plt.savefig(Plot_Dir+"/hj_{}_"+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/hj_{}_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()

print("h")

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["ht"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["ht"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                plt.title(r"$ h$")
                plt.ylim(0,0.45)
                plt.legend(loc='upper left')

                if id_xiag==2:
                    RD_Neu_0 = (res["ht"][res["states"][:, 1]<1.5])[0]
                    RD_Neu_1 = (res["ht"][res["states"][:, 1]<1.5])[-1]
                    print("Neutrality Left={:.4f}".format(RD_Neu_0))
                if id_xiag==1:
                    RD_Less_0 = (res["ht"][res["states"][:, 1]<1.5])[0]
                    RD_Less_1 = (res["ht"][res["states"][:, 1]<1.5])[-1]
                    print("Less Left={:.4f}".format(RD_Less_0))

                if id_xiag==0:
                    RD_More_0 = (res["ht"][res["states"][:, 1]<1.5])[0]
                    RD_More_1 = (res["ht"][res["states"][:, 1]<1.5])[-1]
                    print("More Left={:.4f}".format(RD_More_0))

                    
plt.savefig(Plot_Dir+"/h2_{}_"+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/h2_{}_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()

print("hk")
for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], -( res["hkt"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], -( res["hkt"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                plt.title(r"-$ h_k$")
                plt.ylim(0,0.45)
                plt.legend(loc='upper left')


                if id_xiag==2:
                    RD_Neu_0 = (res["hkt"][res["states"][:, 1]<1.5])[0]
                    RD_Neu_1 = (res["hkt"][res["states"][:, 1]<1.5])[-1]
                    print("Neutrality Left={:.4f}".format(RD_Neu_0))
                if id_xiag==1:
                    RD_Less_0 = (res["hkt"][res["states"][:, 1]<1.5])[0]
                    RD_Less_1 = (res["hkt"][res["states"][:, 1]<1.5])[-1]
                    print("Less Left={:.4f}".format(RD_Less_0))

                if id_xiag==0:
                    RD_More_0 = (res["hkt"][res["states"][:, 1]<1.5])[0]
                    RD_More_1 = (res["hkt"][res["states"][:, 1]<1.5])[-1]
                    print("More Left={:.4f}".format(RD_More_0))


plt.savefig(Plot_Dir+"/hk2_{}_"+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/hk2_{}_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()

print("hj")

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], -( res["hjt"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], -( res["hjt"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                plt.title(r"-$h_j$")
                plt.ylim(0,0.02)

                plt.legend(loc='upper left')


                if id_xiag==2:
                    RD_Neu_0 = (res["hjt"][res["states"][:, 1]<1.5])[0]
                    RD_Neu_1 = (res["hjt"][res["states"][:, 1]<1.5])[-1]
                    print("Neutrality Left={:.4f}".format(RD_Neu_0))
                if id_xiag==1:
                    RD_Less_0 = (res["hjt"][res["states"][:, 1]<1.5])[0]
                    RD_Less_1 = (res["hjt"][res["states"][:, 1]<1.5])[-1]
                    print("Less Left={:.4f}".format(RD_Less_0))

                if id_xiag==0:
                    RD_More_0 = (res["hjt"][res["states"][:, 1]<1.5])[0]
                    RD_More_1 = (res["hjt"][res["states"][:, 1]<1.5])[-1]
                    print("More Left={:.4f}".format(RD_More_0))


plt.savefig(Plot_Dir+"/hj2_{}_"+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/hj2_{}_"+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["RelativeEntropy_hk"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["RelativeEntropy_hk"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                # plt.title("Relative Entropy for Capital")
                # plt.ylim(0,0.45)
                plt.legend(loc='upper left')

# plt.savefig(Plot_Dir+"/h2_{}_".format(IntPeriod)+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/RE_K_{}_".format(IntPeriod)+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["RelativeEntropy_h"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["RelativeEntropy_h"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                # plt.title(" Relative Entropy for Temperature Anomaly")
                # plt.ylim(0,0.45)
                plt.legend(loc='upper left')

# plt.savefig(Plot_Dir+"/h2_{}_".format(IntPeriod)+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/RE_TA_{}_".format(IntPeriod)+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["RelativeEntropy_hj"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["RelativeEntropy_hj"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                # plt.title("Relative Entropy for Knowledge Capital")
                plt.ticklabel_format(style='plain')
                # plt.ylim(0,0.45)
                plt.legend(loc='upper left')

# plt.savefig(Plot_Dir+"/h2_{}_".format(IntPeriod)+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/RE_J_{}_".format(IntPeriod)+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["dL"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["dL"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                # plt.title("Relative Entropy for Knowledge Capital")
                plt.ticklabel_format(style='plain')
                # plt.ylim(0,0.2)
                plt.legend(loc='upper left')

# plt.savefig(Plot_Dir+"/h2_{}_".format(IntPeriod)+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/dL_{}_".format(IntPeriod)+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["dK"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["dK"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                # plt.title("Relative Entropy for Knowledge Capital")
                plt.ticklabel_format(style='plain')
                # plt.ylim(0.8,1)
                plt.legend(loc='upper left')

# plt.savefig(Plot_Dir+"/h2_{}_".format(IntPeriod)+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/dK_{}_".format(IntPeriod)+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  np.exp(res["dL"]))[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  np.exp(res["dL"]))[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                # plt.title("Relative Entropy for Knowledge Capital")
                plt.ticklabel_format(style='plain')
                plt.ylim(1,1.2)
                plt.legend(loc='upper left')

# plt.savefig(Plot_Dir+"/h2_{}_".format(IntPeriod)+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/dL_exp_{}_".format(IntPeriod)+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  np.exp(res["dK"]))[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  np.exp(res["dK"]))[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                # plt.title("Relative Entropy for Knowledge Capital")
                plt.ticklabel_format(style='plain')
                plt.ylim(2.4,2.7)
                plt.legend(loc='upper left')

# plt.savefig(Plot_Dir+"/h2_{}_".format(IntPeriod)+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/dK_exp_{}_".format(IntPeriod)+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["RelativeEntropy_TechJump"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["RelativeEntropy_TechJump"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                plt.title(r"$ h$")
                # plt.ylim(0,0.45)
                plt.legend(loc='upper left')

# plt.savefig(Plot_Dir+"/h2_{}_".format(IntPeriod)+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/RE_TechJump_{}_".format(IntPeriod)+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"], (  res["RelativeEntropy_TechJump"]),label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"], (  res["RelativeEntropy_TechJump"]),label=labellist[id_xiag]  )

                plt.xlabel("Years")
                # plt.title("Relative Entropy for Technology Jump")
                # plt.ylim(0,0.45)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper left')

# plt.savefig(Plot_Dir+"/h2_{}_".format(IntPeriod)+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/RE_TechJump2_{}_".format(IntPeriod)+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()



for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["RelativeEntropy_DamageJump"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"][res["states"][:, 1]<1.5], (  res["RelativeEntropy_DamageJump"])[res["states"][:, 1]<1.5],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                plt.title(r"$ h$")
                # plt.ylim(0,0.45)
                plt.legend(loc='upper left')

# plt.savefig(Plot_Dir+"/h2_{}_".format(IntPeriod)+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/RE_DamageJump_{}_".format(IntPeriod)+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)
                
                
                if xigarr[id_xiag]>10:

                    plt.plot(res["years"], res["RelativeEntropy_DamageJump"],label=labellist[id_xiag]  )
                else:
                    plt.plot(res["years"], res["RelativeEntropy_DamageJump"],label=labellist[id_xiag]  )

                plt.xlabel("Years")
                # plt.title("Relative Entropy for Damage Jump")
                # plt.ylim(0,0.45)
                plt.xlim(0,IntPeriod)
                plt.legend(loc='upper left')

# plt.savefig(Plot_Dir+"/h2_{}_".format(IntPeriod)+Filename+".pdf".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.savefig(Plot_Dir+"/RE_DamageJump2_{}_".format(IntPeriod)+Filename+"_rho={}_delta={}_phi0={}".format(rho,delta,phi_0)+".png".format(IntPeriod, xiaarr,xicarr,xidarr,xigarr,psi0arr,psi1arr,varrhoarr))
plt.close()

plt.style.use('default')
plt.rcParams["lines.linewidth"] = 20
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["figure.figsize"] = (16,10)
plt.rcParams["font.size"] = 25
plt.rcParams["legend.frameon"] = False


print("After, figure default size is: ", plt.rcParams["savefig.bbox"])
print("After, figure default size is: ", plt.rcParams["figure.figsize"])
print("After, figure default dpi is: ", plt.rcParams["figure.dpi"])
print("After, figure default size is: ", plt.rcParams["font.size"])
print("After, legend.frameon is: ", plt.rcParams["legend.frameon"])
print("After, lines.linewidth is: ", plt.rcParams["lines.linewidth"])


for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

                
                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                NUM_DAMAGE = res["gt_dmg"].shape[0]
                gamma_3_list = np.linspace(0., 1./3., NUM_DAMAGE)

                # γ3_distort = np.load("γ3_5.npy")
                print(NUM_DAMAGE)
                
                time_frame = int(plot_year_gamma/timespan)
                γ3_distort = res["gt_dmg"][:, time_frame] 
                # plt.figure(figsize=(16,10))
                plt.hist(gamma_3_list, weights=np.ones(len(gamma_3_list)) / len(gamma_3_list), 
                        alpha=0.5, color="C3", ec="darkgray",label='Baseline' .format(xicarr[id_xiag], xidarr[id_xiag], xigarr[id_xiag]), bins=NUM_DAMAGE)
                plt.hist(gamma_3_list, weights= γ3_distort / np.sum(γ3_distort), 
                        alpha=0.5, color="C0", ec="darkgray",label=labellist[id_xiag], bins=NUM_DAMAGE)
                if vartheta_bar==0.1:
                    plt.ylim(0, 0.15)
                if vartheta_bar==0.5:
                    plt.ylim(0, 0.15)
                plt.ylim(0, 0.3)
                # plt.title("Distorted Probability of Damage Models")
                # plt.xlabel("Damage Curvature")
                plt.legend(loc='upper left',frameon=False)

                    
                # plt.savefig(Plot_Dir+"/Gamma3_{},xia={:.5f},xik={:.3f},xic={:.3f},xij={:.3f},xid={:.3f},xig={:.3f},psi0={:.3f},psi1={:.3f},varrho={:.1f}\%_rho={}_delta={}_phi0={}.png".format(plot_year_gamma,xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho,delta,phi_0))
                plt.savefig(Plot_Dir+"/Gamma3_{},".format(plot_year_gamma)+Filename+labellist2[id_xiag] + "_rho={}_delta={}_phi0={}.png".format(rho,delta,phi_0))
                plt.savefig(Plot_Dir+"/Gamma3_{},".format(plot_year_gamma)+Filename+labellist2[id_xiag] + "_rho={}_delta={}_phi0={}.pdf".format(rho,delta,phi_0))
                plt.close()

for id_xiag in range(len(xiaarr)): 
    for id_psi0 in range(len(psi0arr)):
        for id_psi1 in range(len(psi1arr)):
            for id_varrho in range(len(varrhoarr)):

        
                res = model_simulation_generate(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho)

                time_frame = int(plot_year_theta/timespan)

                # histogram of beta_f
                theta_ell = pd.read_csv("./data/model144.csv", header=None).to_numpy()[:, 0]
                # print("theta_ell")
                # print(theta_ell)
                # print("theta_ell_new")
                # print(theta_ell_new)
                pi_c_o = np.ones(len(theta_ell)) / len(theta_ell)
                # pi_c = np.load("πc_5.npy")
                theta_ell_new = res["theta_ell_new"][:,time_frame]

                pi_c = res["pic_t"][:, time_frame]
                # plt.figure(figsize=(16,10))

                # plt.hist(theta_ell, weights=pi_c_o, bins=np.linspace(0.8, 3., 16), density=True, 
                #         alpha=0.5, ec="darkgrey", color="C3",label='Baseline')
                # plt.hist(theta_ell, weights=pi_c, bins=np.linspace(0.8, 3., 16), density=True, 
                #         alpha=0.5, ec="darkgrey", color="C0",label='$\\xi_a={:.4f}$,$\\xi_g=\\xi_d=\\xi_r={:.3f}$'.format(xicarr[id_xiag], xidarr[id_xiag], xigarr[id_xiag])  )
                # plt.legend(loc='upper left')
                # plt.title("Distorted probability of Climate Models")

                # plt.ylim(0, 1.4)
                # plt.xlabel("Climate Sensitivity")
                # plt.savefig("./abatement/pdf_2tech/"+args.dataname+"/ClimateSensitivity_25,xia={:.4f},xic={:.3f},xid={:.3f},xig={:.3f},psi0={:.3f},psi1={:.3f},BC.pdf".format(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1]))
                # plt.savefig("./abatement/pdf_2tech/"+args.dataname+"/ClimateSensitivity_25,xia={:.4f},xic={:.3f},xid={:.3f},xig={:.3f},psi0={:.3f},psi1={:.3f},BC.png".format(xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1]))
                # plt.close()

                plt.figure(figsize=(16,10))

                plt.hist(theta_ell, weights=pi_c_o, bins=np.arange(0.9, 2.85, 0.15), density=True, 
                        alpha=0.5, ec="darkgrey", color="C3",label='Baseline')
                plt.hist(theta_ell_new*1000, weights=pi_c, bins=np.arange(0.9, 2.85, 0.15), density=True, 
                        alpha=0.5, ec="darkgrey", color="C0",label=labellist[id_xiag]  )
                # plt.hist(theta_ell, weights=pi_c, bins=np.linspace(0.8, 3., 16), density=True, 
                #         alpha=0.5, ec="darkgrey", color="C0",label=labellist[id_xiag]  )
                plt.legend(loc='upper left')
                # plt.title("Distorted Probability of Climate Models")


                print("mean of uncondition = {}" .format(np.average(theta_ell,weights = pi_c_o)))
                print("mean of condition = {}" .format(np.average(theta_ell,weights = pi_c)))
                print("mean of condition mean shift = {}" .format(np.average(theta_ell_new*1000,weights = pi_c)))
                    

                plt.ylim(0, 1.4)
                # plt.xlabel("Climate Sensitivity")
                # plt.savefig(Plot_Dir+"/ClimateSensitivity_pmean_{},xia={:.5f},xik={:.3f},xic={:.3f},xij={:.3f},xid={:.3f},xig={:.3f},psi0={:.3f},psi1={:.3f},varrho={:.1f}\%_rho={}_delta={}_phi0={}.pdf".format(plot_year_theta, xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho,delta,phi_0))
                # plt.savefig(Plot_Dir+"/ClimateSensitivity_pmean_{},xia={:.5f},xik={:.3f},xic={:.3f},xij={:.3f},xid={:.3f},xig={:.3f},psi0={:.3f},psi1={:.3f},varrho={:.1f}\%_rho={}_delta={}_phi0={}.png".format(plot_year_theta, xiaarr[id_xiag],xikarr[id_xiag],xicarr[id_xiag],xijarr[id_xiag],xidarr[id_xiag],xigarr[id_xiag],psi0arr[id_psi0],psi1arr[id_psi1], varrhoarr[id_varrho],rho,delta,phi_0))
                plt.savefig(Plot_Dir+"/ClimateSensitivity_pmean_{},xic={:.3f}".format(plot_year_theta,xicarr[id_xiag]) +Filename+labellist2[id_xiag]+"_rho={}_delta={}_phi0={}.png".format(rho,delta,phi_0))
                plt.savefig(Plot_Dir+"/ClimateSensitivity_pmean_{},xic={:.3f}".format(plot_year_theta,xicarr[id_xiag]) +Filename+labellist2[id_xiag]+"_rho={}_delta={}_phi0={}.pdf".format(rho,delta,phi_0))
                plt.close()
