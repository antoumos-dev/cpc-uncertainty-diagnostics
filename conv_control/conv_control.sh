#!/bin/bash
#SBATCH --partition=postproc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=14:00:00
#SBATCH --job-name=conv_control


module use -a /user-environment/modules
module use -a /mch-environment/v6/modules
module use -a /usr/share/modules
module use -a /usr/share/Modules/3.2.10/modulefiles
module use -a /store_new/mch/msclim/share/modulefiles/modules/all
module load udunits proj sqlite r gdal geos hdf5 cats

export TEXMFHOME=/store_new/mch/msclim/share/CATs/TinyTex
export PATH=/store_new/mch/msclim/share/CATs/TinyTex/bin/x86_64-linux/:${PATH}

#mkdir -p logs

YEAR=${1:-25}

cd /store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats/conv_control
srun Rscript conv_control.r "$YEAR"
