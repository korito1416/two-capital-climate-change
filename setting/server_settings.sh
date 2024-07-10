# server_settings.sh
#SBATCH --account=pi-lhansen
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --exclude=mcn53,mcn55,mcn57,mcn08

module load python/booth/3.10  
module load gcc/9.2.0

