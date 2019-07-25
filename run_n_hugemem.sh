#!/bin/bash
#PBS -P n74
#PBS -q hugemem
#PBS -l walltime=32:00:00
#PBS -lmem=1024GB
#PBS -lncpus=7
#PBS -l wd

module load intel-cc/12.1.9.293
module load intel-fc/12.1.9.293
module load gdal
module load geos/3.5.0
module load hdf5/1.8.10
#module load openmpi/1.6.3
module unload python
module load python/2.7.11
module load python/2.7.11-matplotlib

# To get rtree to run
export SPATIALINDEX_C_LIBRARY=/short/n74/jdg547/spatialindex-src-1.8.5/lib/libspatialindex_c.so.4
export LD_LIBRARY_PATH=/short/n74/jdg547/spatialindex-src-1.8.5/lib:$LD_LIBRARY_PATH
# Python paths for local openquake installs and dependencies
export PYTHONPATH=.:/home/547/jdg547/.local/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=.:/short/w84/NSHA18/sandpit/jdg547/oq-hazardlib:${PYTHONPATH}
export PYTHONPATH=.:/short/w84/NSHA18/sandpit/jdg547/oq-engine:${PYTHONPATH}
export PYTHONPATH=.:/short/n74/src/lib/python/:${PYTHONPATH}

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
    if [ $counter = 17 ];
    then
        break
    fi
done      

wait
