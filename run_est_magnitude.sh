#PBS -P n74
#PBS -q normal
#PBS -l walltime=18:00:00
#PBS -l ncpus=1
#PBS -l mem=128GB
#PBS -l wd

#New
module load intel-mkl/2020.3.304
module load geos/3.8.0
module load hdf5/1.10.7
module load openmpi/4.1.4
module load python3/3.9.2
# Load gdal after python to avoid conflict                                                                                                                             
module load gdal/3.5.0

#NEW
# Local pythonpaths                                                                                                                                                                        
export PYTHONPATH=.:/home/547/jdg547/.local/lib/python3.9/site-packages/:${PYTHONPATH}
export PYTHONPATH=.::/home/547/jdg547/modelling/EQIAT/:${PYTHONPATH}                                                                                     
export PYTHONPATH=.::/home/547/jdg547/modelling/oq-engine/:${PYTHONPATH}
export PYTHONPATH=.::/scratch/w84/jdg547/:${PYTHONPATH}

#python estimate_magnitude.py -param_file data/1699slab_params.txt >& data/1699slab_params.txt.log
python3 estimate_magnitude.py -param_file data/1699megathrust_params.txt >& data/1699megathrust_params.txt.log 
