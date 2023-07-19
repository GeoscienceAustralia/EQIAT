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
#python estimate_magnitude.py -param_file data/1867slab_params.txt >& data/1867slab_params.txt.log
#python estimate_magnitude.py -param_file data/1867megathrust_params.txt >& data/1867megathrust_params.txt.log
#python estimate_magnitude.py -param_file data/1867_params.txt >& data/1867_params.txt.log
#python estimate_magnitude.py -param_file data/1847_params.txt >& data/1847_params.txt.log
#python estimate_magnitude.py -param_file data/1840_params.txt >& data/1840_params.txt.log
#python estimate_magnitude.py -param_file data/1834_params.txt >& data/1834_params.txt.log
#python estimate_magnitude.py -param_file data/1780megathrust_params.txt >& data/1780megathrust_params.txt.log
#python estimate_magnitude.py -param_file data/1780_params.txt >& data/1780_params.txt.log
#python3 estimate_magnitude.py -param_file data/1852Banda_domain_ryan_mmi_params.txt >& data/1852Banda_domain_ryan_mmi_params.txt.log
#python estimate_magnitude.py -param_file data/1852Banda_area_params.txt >& data/1852_area_params.txt.log
#python estimate_magnitude.py -param_file data/1699megathrust_params.txt >& data/1699megathrust_params.txt.log
#python estimate_magnitude.py -param_file data/1820_M8.4_params.txt >& data/1820_M8.4_params.txt.log
#python estimate_magnitude.py -param_file data/1820_params.txt >& data/1820_params.txt.log
#python estimate_magnitude.py -param_file data/1815_params.txt >& data/1815_params.txt.log
#python estimate_magnitude.py -param_file data/1818_params.txt >& data/1818_params.txt.log
#python estimate_magnitude.py -param_file data/2006_params.txt >& data/2006_params.txt.log
#python estimate_magnitude.py -param_file data/2017slab_params.txt >& data/2017slab_params.txt.log
#python estimate_magnitude.py -param_file data/2018_params.txt >& data/2018_params.txt.log
