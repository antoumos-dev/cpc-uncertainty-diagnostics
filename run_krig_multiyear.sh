#!/bin/bash
#SBATCH --partition=postproc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00

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

set -euo pipefail
mkdir -p logs

# Years
START_YEAR=2016
END_YEAR=2025

BANDS="0.01-0.1"

srun Rscript run_krig_multiyear.r \
    $START_YEAR \
    $END_YEAR \
    $BANDS
