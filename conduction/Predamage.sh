#! /bin/bash



current_dir=$(pwd)
epsilonarraypost=(0.1) 
epsilonarraypre=(0.1)

python_name="Predamage.py"
python_dir="./python"
output_dir="/scratch/pengyu/"

NUM_DAMAGE=20
ID_MAX_DAMAGE=$((NUM_DAMAGE - 1))
maxiterarr=(100000 800000)

declare -A hXarr1=([0]=0.2 [1]=0.2 [2]=0.2)
declare -A hXarr2=([0]=0.1 [1]=0.1 [2]=0.1)
hXarrays=(hXarr1)
# hXarrays=(hXarr2)


Xminarr=(4.00 0.0 1.0 0.0)
Xmaxarr=(9.00 4.0 6.0 3.0)


all_channel=$1
single_channel_more_aversion=$2
single_channel_less_aversion=$3
figure4_search_xia=$4
appendix=$5

# Setting xi arrays based on flags
if [ "$all_channel" = true ]; then
    xi_a=(100000. 100000. 100000.)
    xi_k=(0.075 0.150 100000.)
    xi_c=(0.075 0.150 100000.)
    xi_j=(0.075 0.150 100000.)
    xi_d=(0.075 0.150 100000.)
    xi_g=(0.075 0.150 100000.)
elif [ "$single_channel_more_aversion" = true ]; then
    xi_a=(100000. 100000. 100000. 100000.)
    xi_k=(0.075 100000. 100000. 100000.)
    xi_c=(100000. 0.075 100000. 100000.)
    xi_j=(100000. 100000. 0.075 100000.)
    xi_d=(100000. 100000. 100000. 0.075)
    xi_g=(100000. 100000. 0.075 100000.)
elif [ "$single_channel_less_aversion" = true ]; then
    xi_a=(100000. 100000. 100000. 100000.)
    xi_k=(0.150 100000. 100000. 100000.)
    xi_c=(100000. 0.150 100000. 100000.)
    xi_j=(100000. 100000. 0.150 100000.)
    xi_d=(100000. 100000. 100000. 0.150)
    xi_g=(100000. 100000. 0.150 100000.)
elif [ "$figure4_search_xia" = true ]; then
    xi_a=(0.0025 0.0030 100000.)
    xi_k=(0.075 0.150 100000.)
    xi_c=(100000. 100000. 100000.)
    xi_j=(0.075 0.150 100000.)
    xi_d=(0.075 0.150 100000.)
    xi_g=(0.075 0.150 100000.)
elif [ "$appendix" = true ]; then
    xi_a=(100000. 100000. 100000. 100000.)
    xi_k=(0.005 100000. 100000. 100000.)
    xi_c=(0.005 100000. 100000. 100000.)
    xi_j=(0.005 0.005 0.075 0.150)
    xi_d=(0.005 100000. 100000. 100000.)
    xi_g=(0.005 0.005 0.075 0.150)
else
    echo "No valid condition set"
    exit 1
fi


varrhoarr=(1120)


rhoarr=(1)
deltaarr=(0.010)

LENGTH_rho=$((${#rhoarr[@]} - 1))

psi0arr=(0.105830)
psi1arr=(0.5)
LENGTH_psi=$((${#psi0arr[@]} - 1))
LENGTH_xi=$((${#xi_a[@]} - 1))

phi0arr=(0.5)

LENGTH_phi0=$((${#phi0arr[@]} - 1))


for epsilon in ${epsilonarraypre[@]}; do
	for epsilonpost in ${epsilonarraypost[@]}; do
		for hXarri in "${hXarrays[@]}"; do
        for phi0index in $(seq 0 $LENGTH_phi0); do
			count=0
			declare -n hXarr="$hXarri"

			action_name="2jump_step_${Xminarr[0]},${Xmaxarr[0]}_${Xminarr[1]},${Xmaxarr[1]}_${Xminarr[2]},${Xmaxarr[2]}_${Xminarr[3]},${Xmaxarr[3]}_SS_${hXarr[0]},${hXarr[1]},${hXarr[2]}_LR_${epsilonpost}_Current_phi0_${phi0arr[$phi0index]}"
			epsilonarr=(0.1 ${epsilon})
			fractionarr=(0.1 ${epsilon})


			for PSI_0 in ${psi0arr[@]}; do
				for PSI_1 in ${psi1arr[@]}; do
					for varrho in ${varrhoarr[@]}; do
					for j in $(seq 0 $LENGTH_xi); do
							for k in $(seq 0 $LENGTH_rho); do

							mkdir -p ./job-outs/${action_name}/Pre/xia_${xi_a[$j]}_xik_${xi_k[$j]}_xic_${xi_c[$j]}_xij_${xi_j[$j]}_xid_${xi_d[$j]}_xig_${xi_g[$j]}_PSI0_${PSI_0}_PSI1_${PSI_1}_varrho_${varrho}_rho_${rhoarr[$k]}_delta_${deltaarr[$k]}/

							if [ -f ./bash/${action_name}/hX_${hXarr[0]}_xia_${xi_a[$j]}_xik_${xi_k[$j]}_xic_${xi_c[$j]}_xij_${xi_j[$j]}_xid_${xi_d[$j]}_xig_${xi_g[$j]}_PSI0_${PSI_0}_PSI1_${PSI_1}_varrho_${varrho}_rho_${rhoarr[$k]}_delta_${deltaarr[$k]}_Eps_${epsilon}.sh ]; then
								rm ./bash/${action_name}/hX_${hXarr[0]}_xia_${xi_a[$j]}_xik_${xi_k[$j]}_xic_${xi_c[$j]}_xij_${xi_j[$j]}_xid_${xi_d[$j]}_xig_${xi_g[$j]}_PSI0_${PSI_0}_PSI1_${PSI_1}_varrho_${varrho}_rho_${rhoarr[$k]}_delta_${deltaarr[$k]}_Eps_${epsilon}.sh
							fi

							mkdir -p ./bash/${action_name}/

							touch ./bash/${action_name}/hX_${hXarr[0]}_xia_${xi_a[$j]}_xik_${xi_k[$j]}_xic_${xi_c[$j]}_xij_${xi_j[$j]}_xid_${xi_d[$j]}_xig_${xi_g[$j]}_PSI0_${PSI_0}_PSI1_${PSI_1}_varrho_${varrho}_rho_${rhoarr[$k]}_delta_${deltaarr[$k]}_Eps_${epsilon}.sh

							tee -a ./bash/${action_name}/hX_${hXarr[0]}_xia_${xi_a[$j]}_xik_${xi_k[$j]}_xic_${xi_c[$j]}_xij_${xi_j[$j]}_xid_${xi_d[$j]}_xig_${xi_g[$j]}_PSI0_${PSI_0}_PSI1_${PSI_1}_varrho_${varrho}_rho_${rhoarr[$k]}_delta_${deltaarr[$k]}_Eps_${epsilon}.sh <<EOF
#! /bin/bash

######## login
#SBATCH --job-name=${Xminarr[1]}_${hXarr[0]}_xia_${xi_a[$j]}_xik_${xi_k[$j]}_xic_${xi_c[$j]}_xij_${xi_j[$j]}_xid_${xi_d[$j]}_xig_${xi_g[$j]}_${rhoarr[$k]}_phi0_${phi0arr[$phi0index]}_${epsilon}
#SBATCH --output=./job-outs/${action_name}/Pre/xia_${xi_a[$j]}_xik_${xi_k[$j]}_xic_${xi_c[$j]}_xij_${xi_j[$j]}_xid_${xi_d[$j]}_xig_${xi_g[$j]}_PSI0_${PSI_0}_PSI1_${PSI_1}_varrho_${varrho}_rho_${rhoarr[$k]}_delta_${deltaarr[$k]}/mercury_pre_${epsilon}.out
#SBATCH --error=./job-outs/${action_name}/Pre/xia_${xi_a[$j]}_xik_${xi_k[$j]}_xic_${xi_c[$j]}_xij_${xi_j[$j]}_xid_${xi_d[$j]}_xig_${xi_g[$j]}_PSI0_${PSI_0}_PSI1_${PSI_1}_varrho_${varrho}_rho_${rhoarr[$k]}_delta_${deltaarr[$k]}/mercury_pre_${epsilon}.err

#SBATCH --account=pi-lhansen
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=7-00:00:00
#SBATCH --exclude=mcn53,mcn55,mcn57,mcn08


module purge
####### load modules
module load python/booth/3.10  gcc/9.2.0

echo "\$SLURM_JOB_NAME"

echo "Program starts \$(date)"
start_time=\$(date +%s)
# perform a task

srun python3 ${python_dir}/$python_name  --outputname ${output_dir} --num_gamma $NUM_DAMAGE --xi_a ${xi_a[$j]} --xi_k ${xi_k[$j]} --xi_c ${xi_c[$j]} --xi_j ${xi_j[$j]} --xi_d ${xi_d[$j]} --xi_g ${xi_g[$j]}   --epsilonarr ${epsilonarr[@]}  --fractionarr ${fractionarr[@]}   --maxiterarr ${maxiterarr[@]}  --psi_0 $PSI_0 --psi_1 $PSI_1    --name ${action_name} --hXarr ${hXarr[@]} --Xminarr ${Xminarr[@]} --Xmaxarr ${Xmaxarr[@]} --varrho ${varrho}  --phi_0 ${phi0arr[$phi0index]} --rho ${rhoarr[$k]} --delta ${deltaarr[$k]}

echo "Program ends \$(date)"
end_time=\$(date +%s)

# elapsed time with second resolution
elapsed=\$((end_time - start_time))

eval "echo Elapsed time: \$(date -ud "@\$elapsed" +'\$((%s/3600/24)) days %H hr %M min %S sec')"

EOF
									count=$(($count + 1))
									sbatch ./bash/${action_name}/hX_${hXarr[0]}_xia_${xi_a[$j]}_xik_${xi_k[$j]}_xic_${xi_c[$j]}_xij_${xi_j[$j]}_xid_${xi_d[$j]}_xig_${xi_g[$j]}_PSI0_${PSI_0}_PSI1_${PSI_1}_varrho_${varrho}_rho_${rhoarr[$k]}_delta_${deltaarr[$k]}_Eps_${epsilon}.sh
								done
							done
						done
					done
				done
			done
		done
	done
done