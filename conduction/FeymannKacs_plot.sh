#! /bin/bash

# actiontime=1
epsilonarraypost=(0.1) 

python_name_unit="Result_2jump_UD_simulate_CRS_long_rhodelta3_FKPlot_Withuncertainty_addpartialg_improveM_adddamage.py"
python_dir="/home/bcheng4/TwoCapital_Shrink/abatement_UD"
output_dir="/scratch/bincheng/"

NUM_DAMAGE=20

declare -A hXarr1=([0]=0.2 [1]=0.2 [2]=0.2)
declare -A hXarr2=([0]=0.1 [1]=0.1 [2]=0.1)
hXarrays=(hXarr1)
# hXarrays=(hXarr2)

Xminarr=(4.00 0.0 1.0 0.0)
Xmaxarr=(9.00 4.0 6.0 3.0)


# xi_a=(100000. 100000. 100000.)
# xi_k=(0.025  0.050 100000.)
# xi_c=(0.025  0.050 100000.)
# xi_j=(0.025  0.050 100000.)
# xi_d=(0.025  0.050 100000.)
# xi_g=(0.025  0.050 100000.)


xi_a=(100000. 100000. 100000.)
xi_k=(0.075  0.150 100000.)
xi_c=(0.075  0.150 100000.)
xi_j=(0.075  0.150 100000.)
xi_d=(0.075  0.150 100000.)
xi_g=(0.075  0.150 100000.)



# xi_a=(100000. 100000. 100000. 100000.)
# xi_k=(0.075 100000. 100000. 100000.)
# xi_c=(100000. 0.075 100000. 100000.)
# xi_j=(100000. 100000. 0.075 100000.)
# xi_d=(100000. 100000. 100000. 0.075)
# xi_g=(100000. 100000. 0.075 100000.)

# xi_a=(100000. 100000. 100000. 100000.)
# xi_k=(0.150 100000. 100000. 100000.)
# xi_c=(100000. 0.150 100000. 100000.)
# xi_j=(100000. 100000. 0.150 100000.)
# xi_d=(100000. 100000. 100000. 0.150)
# xi_g=(100000. 100000. 0.150 100000.)

varrhoarr=(1120)


psi0arr=(0.105830)
psi1arr=(0.5)



# rhoarr=(0.66 1 1.5)
# deltaarr=(0.010 0.010 0.010)

rhoarr=(1)
deltaarr=(0.010)


# rhoarr=(1 1 1)
# deltaarr=(0.010 0.020 0.030)




LENGTH_rho=$((${#rhoarr[@]} - 1))


phi0arr=(0.5)
# phi0arr=(0.1)
LENGTH_phi0=$((${#phi0arr[@]} - 1))


server_name="mercury"

LENGTH_psi=$((${#psi0arr[@]} - 1))
LENGTH_xi=$((${#xi_a[@]} - 1))

hXarr_SG=(0.2 0.2 0.2)
Xminarr_SG=(4.00 0.0 -5.5 0.0)
Xmaxarr_SG=(9.00 4.0 0.0 3.0)
interp_action_name="2jump_step_0.2_0.2_0.2_LR_0.01"
fstr_SG="NearestNDInterpolator"


scheme_array=("direct")
HJBsolution_array=("direct")


LENGTH_scheme=$((${#scheme_array[@]} - 1))


auto=1
# year=25
# year=40
# year=45
year=50
# year=55
# year=130
# year=500
# year=1500
# year=500
# m0_array=("Temp")
m0_array="Temperature"
# m0_array="Technology"

LENGTH_m0_array=$((${#m0_array[@]} - 1))
for epsilonpost in ${epsilonarraypost[@]}; do
    for hXarri in "${hXarrays[@]}"; do
        for phi0index in $(seq 0 $LENGTH_phi0); do

        count=0
        declare -n hXarr="$hXarri"
        action_name="2jump_step_${Xminarr[0]},${Xmaxarr[0]}_${Xminarr[1]},${Xmaxarr[1]}_${Xminarr[2]},${Xmaxarr[2]}_${Xminarr[3]},${Xmaxarr[3]}_SS_${hXarr[0]},${hXarr[1]},${hXarr[2]}_LR_${epsilonpost}_Current_phi0_${phi0arr[$phi0index]}"

        for PSI_0 in ${psi0arr[@]}; do
            for PSI_1 in ${psi1arr[@]}; do
                for varrho in ${varrhoarr[@]}; do
                    for j in $(seq 0 $LENGTH_xi); do
                        for k in $(seq 0 $LENGTH_scheme); do
							for kk in $(seq 0 $LENGTH_rho); do

                    mkdir -p ./job-outs4/${action_name}/Graph_Simulate_plot/scheme_${scheme_array[$k]}_HJB_${HJBsolution_array[$k]}/

                    if [ -f ./bash/${action_name}/hX_${hXarr[0]}_PSI0_${PSI_0}_PSI1_${PSI_1}_varrho_${varrho}_rho_${rhoarr[$kk]}_delta_${deltaarr[$kk]}_Graph.sh ]; then
                        rm ./bash/${action_name}/hX_${hXarr[0]}_PSI0_${PSI_0}_PSI1_${PSI_1}_varrho_${varrho}_rho_${rhoarr[$kk]}_delta_${deltaarr[$kk]}_Graph.sh
                    fi
                    mkdir -p ./bash/${action_name}/

                    touch ./bash/${action_name}/hX_${hXarr[0]}_PSI0_${PSI_0}_PSI1_${PSI_1}_varrho_${varrho}_rho_${rhoarr[$kk]}_delta_${deltaarr[$kk]}_Graph.sh

                    tee -a ./bash/${action_name}/hX_${hXarr[0]}_PSI0_${PSI_0}_PSI1_${PSI_1}_varrho_${varrho}_rho_${rhoarr[$kk]}_delta_${deltaarr[$kk]}_Graph.sh <<EOF
#! /bin/bash


######## login 
#SBATCH --job-name=plot_${year}_${m0_array}
#SBATCH --output=./job-outs4/${action_name}/Graph_Simulate_plot/scheme_${scheme_array[$k]}_HJB_${HJBsolution_array[$k]}/graph_${python_name_unit}_${year}_${m0_array}_${xi_k[$j]}.out
#SBATCH --error=./job-outs4/${action_name}/Graph_Simulate_plot/scheme_${scheme_array[$k]}_HJB_${HJBsolution_array[$k]}/graph_${python_name_unit}_${year}_${m0_array}_${xi_k[$j]}.err

#SBATCH --account=pi-lhansen
#SBATCH --partition=highmem
#SBATCH --cpus-per-task=1
#SBATCH --mem=52G
#SBATCH --time=0-02:00:00
#SBATCH --exclude=mcn53,mcn55,mcn57,mcn08

####### load modules
module load python/booth/3.8  gcc/9.2.0


echo "\$SLURM_JOB_NAME"
echo "Program starts \$(date)"
start_time=\$(date +%s)

python3 ${python_dir}/${python_name_unit} --dataname  ${action_name}  --outputname ${output_dir} --pdfname ${server_name} --psi0 ${PSI_0} --psi1 ${PSI_1} --xiaarr ${xi_a[$j]} --xikarr ${xi_k[$j]}  --xicarr ${xi_c[$j]}  --xijarr ${xi_j[$j]} --xidarr ${xi_d[$j]} --xigarr ${xi_g[$j]}     --hXarr ${hXarr[@]} --Xminarr ${Xminarr[@]} --Xmaxarr ${Xmaxarr[@]} --auto $auto --IntPeriod ${year} --scheme ${scheme_array[$k]}  --HJB_solution ${HJBsolution_array[$k]}  --varrho ${varrho}   --phi_0 ${phi0arr[$phi0index]}  --rhoarr ${rhoarr[$kk]} --delta ${deltaarr[$kk]} --m0 ${m0_array}

echo "Program ends \$(date)"
end_time=\$(date +%s)

# elapsed time with second resolution
elapsed=\$((end_time - start_time))

eval "echo Elapsed time: \$(date -ud "@\$elapsed" +'\$((%s/3600/24)) days %H hr %M min %S sec')"

EOF

                    sbatch ./bash/${action_name}/hX_${hXarr[0]}_PSI0_${PSI_0}_PSI1_${PSI_1}_varrho_${varrho}_rho_${rhoarr[$kk]}_delta_${deltaarr[$kk]}_Graph.sh

                                    done
                                done
                            done
                        done
                    done
                done
            done
        done
    done
# done
# done