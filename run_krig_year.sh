#!/bin/bash
#SBATCH --partition=postproc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-0

# Activate the user environment
#bash
source /users/antoumos/.local/bin/activate-uenv
uenv start --view=climana climana/24.10:rc1

# Add required module paths and load necessary modules
module use -a /user-environment/modules
module use -a /mch-environment/v6/modules
module use -a /usr/share/modules
module use -a /usr/share/Modules/3.2.10/modulefiles
module use -a /store_new/mch/msclim/share/modulefiles/modules/all
module use -a /oprusers/osm/opr.emme/modules/modulefiles
module load udunits
module load proj
module load sqlite
module load r
module load gdal
module load geos
module load hdf5
module load cats

# Set environment variables
export TEXMFHOME=/store_new/mch/msclim/share/CATs/TinyTex
export PATH=/store_new/mch/msclim/share/CATs/TinyTex/bin/x86_64-linux/:${PATH}
export PATH=./:${PATH}

#set -euo pipefail
mkdir -p logs

#cd /store_new/mch/msclim/antoumos/R/develop/CPC/new_project/

# List of years corresponding to array indices
YEARS=(2019)
YEAR=${YEARS[$SLURM_ARRAY_TASK_ID]}

#MODE="relunc"
#MU_MIN="0.1"

echo "Running year: $YEAR (task $SLURM_ARRAY_TASK_ID)"

srun Rscript run_krig_year.r "$YEAR"