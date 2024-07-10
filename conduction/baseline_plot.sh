#!/bin/bash

temp_script=$(mktemp)
python_name_unit="frequency_plot.py"
python_dir="./python/"

# Create the temporary script
echo "#!/bin/bash" > $temp_script
cat ./setting/server_settings.sh >> $temp_script
cat << EOF >> $temp_script
python3 ${python_dir}/${python_name_unit}
EOF

# Submit the SLURM job
sbatch $temp_script

# Remove the temporary script after submission
rm $temp_script
