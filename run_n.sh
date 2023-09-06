#!/bin/bash
#PBS -P n74
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -lmem=128GB
#PBS -lncpus=4
#PBS -l wd


module load intel-mkl/2020.3.304
module load geos/3.8.0
module load hdf5/1.10.7
module load openmpi/4.1.4
module load python3/3.9.2
# Load gdal after python to avoid conflict
module load gdal/3.5.0

# Local pythonpaths
export PYTHONPATH=.:/home/547/jdg547/.local/lib/python3.9/site-packages/:${PYTHONPATH}
export PYTHONPATH=.::/home/547/jdg547/modelling/EQIAT/:${PYTHONPATH}                                                     
export PYTHONPATH=.::/home/547/jdg547/modelling/oq-engine/:${PYTHONPATH}
export PYTHONPATH=.::/scratch/w84/jdg547/:${PYTHONPATH}

# Script to submit several single cpu jobs at once
counter=0
one=1
# List all subdirectories at the level of individual tsunami runs
all_param_files=$(ls data/*1852*mmi_params.txt)

#mybasedir=$(pwd)

# Loop over all subdirectories
for i in $all_param_files; do
    log_file=$i'.log'
    python estimate_magnitude.py -param_file $i > $log_file &
    counter=$(($counter+$one));
    # Once we have submitted n jobs, break from this loop.
    if [ $counter = 15 ];
    then
        break
    fi
done      

wait
