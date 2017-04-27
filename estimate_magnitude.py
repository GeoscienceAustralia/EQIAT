"""Use ground motion simulations to find the best fitting earthquake
source model based on historical intensity data

Jonathan Griffin
Geoscience Australia, December 2016
"""

import os, sys
import numpy as np
import argparse

from gmf_calculator import get_pt_sources, get_sources, RuptureGmf
from get_site_model import get_site_collection
from openquake.hazardlib.site import SiteCollection
from openquake.hazardlib import nrml
from openquake.hazardlib.gsim.campbell_bozorgnia_2008 import CampbellBozorgnia2008
from openquake.hazardlib.gsim.chiou_youngs_2008 import ChiouYoungs2008
from openquake.hazardlib.gsim.boore_atkinson_2008 import BooreAtkinson2008
from openquake.hazardlib.gsim.youngs_1997 import YoungsEtAl1997SInter, YoungsEtAl1997SSlab
from openquake.hazardlib.gsim.atkinson_boore_2003 import AtkinsonBoore2003SSlab, \
    AtkinsonBoore2003SSlabCascadia, AtkinsonBoore2003SInter
from openquake.hazardlib.gsim.zhao_2006 import ZhaoEtAl2006SSlab, ZhaoEtAl2006SInter
from openquake.hazardlib.geo.point import Point

parser = argparse.ArgumentParser(
        description='Estimate magnitude and location based on historical intensity data')
parser.add_argument('-param_file', type=str, default="",
                        help='Filename containing model parameters')
args = parser.parse_args()
f_in = open(args.param_file, 'r')
event_name = f_in.readline().rstrip().split(',')[1]
area_source_file = f_in.readline().rstrip().split(',')[1]
trt = f_in.readline().rstrip().split(',')[1]
site_model_file = f_in.readline().rstrip().split(',')[1]
site_file = f_in.readline().rstrip().split(',')[1]
site_file_pts = f_in.readline().rstrip().split(',')[1]
# Area source model is used to define the parameter space to be searched
# Note this should be made more generic to use faults as well
#area_source_file = 'data/java_slab_source_model.xml'
#event_name = '1699'
#paramfile = 'params.txt'
#gsim_list = [ChiouYoungs2008(), BooreAtkinson2008()]
#gsim = AtkinsonBoore2003SSlab()
gsim_list = [ZhaoEtAl2006SSlab(), AtkinsonBoore2003SSlab()]
#trt = 'Active'
if trt == 'Subduction Intraslab':
    gsim_list = [ZhaoEtAl2006SSlab(), AtkinsonBoore2003SSlab(), AtkinsonBoore2003SSlabCascadia()]
if trt == 'Active':
    gsim_list = [ChiouYoungs2008(), BooreAtkinson2008(), CampbellBozorgnia2008()]
if trt == 'Subduction Interface':
    gsim_list = [YoungsEtAl1997SInter(), AtkinsonBoore2003SInter(), ZhaoEtAl2006SInter()]

def build_site_col(sites_data, site_model_file):
    """Interpolate vs30 values to sites of interest
    """
#    sites_data = np.genfromtxt(site_file)
    site_points = []
    for i in range(len(sites_data[:,0])):
        site_pt = Point(sites_data[i,0], sites_data[i,1])
        site_points.append(site_pt)
    sitecol = get_site_collection(site_model_file, site_points)
    return sitecol

# Background site class model
#Using inferred vs30, z1.0 and z2.5 values                                                               
#site_model_file = 'data/jawa_site_model.xml'

#Build sites for observations
#site_file = 'data/1699HMMI_weighted_mod.txt'
sites_data = np.genfromtxt(site_file)
sitecol = build_site_col(sites_data, site_model_file)
# Build sites for scenarios
#site_file_pts = 'data/jawa_sites_thin.csv'
sites_data_pts = np.genfromtxt(site_file_pts, delimiter=',')
site_col_scenario = build_site_col(sites_data_pts, site_model_file)


"""
class SiteModel(object):
    def __init__(self, reference_vs30_type,
                 reference_vs30_value,
                 reference_depth_to_2pt5km_per_sec,
                 reference_depth_to_1pt0km_per_sec):
        self.reference_vs30_type = reference_vs30_type
        self.reference_vs30_value = reference_vs30_value
        self.reference_depth_to_2pt5km_per_sec = reference_depth_to_2pt5km_per_sec
        self.reference_depth_to_1pt0km_per_sec = reference_depth_to_1pt0km_per_sec
        self.reference_backarc = None

sitemodel = SiteModel(reference_vs30_type,
                      reference_vs30_value,
                      reference_depth_to_2pt5km_per_sec,
                      reference_depth_to_1pt0km_per_sec)

sitecol = SiteCollection.from_points(sites_data[:,0], sites_data[:,1], sitemodel)             
"""

output_dir = 'outputs'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
#Generate pt sources
pt_sources = get_sources(area_source_file, discretisation = 10)
#Loop through source types
for key,item in pt_sources.iteritems():
#    print key
    # Different key for different tectonic regions
    for gsim in gsim_list:
        rupture_gmfs = RuptureGmf(pt_sources[key], gsim, sitecol)#, 
                              #imts = ['0.0', '1.0'])
        rupture_gmfs.calculate()
        #    print rupture_gmfs.rupture_list
        print rupture_gmfs.gmf_list[-1]
        # Convert to MMI
        rupture_gmfs.rsa2mmi()
        print rupture_gmfs.mmi_list[-1]
    #rupture_gmfs.calc_sum_squares_mmi_weighted(sites_data[:,2]) 
    #    rupture_gmfs.calc_sum_squares_mmi(sites_data[:,2])
        rupture_gmfs.calc_rmse(sites_data[:,2], weights = sites_data[:,3])
    #    print rupture_gmfs.sum_squares_list
        rupture_gmfs.find_best_fit()
        print 'Results using %s GMM' % gsim
        print 'Best magnitude', rupture_gmfs.best_rupture.mag
        print rupture_gmfs.best_rupture.hypocenter
    #    print rupture_gmfs.best_rupture.surface
        print 'Best dip', rupture_gmfs.best_rupture.surface.get_dip()
        print 'Best strike', rupture_gmfs.best_rupture.surface.get_strike()
        print 'RMSE', rupture_gmfs.min_rmse
    # Calculate scenario for best rupture
        rupture_gmfs.calculate_from_rupture(rupture_gmfs.best_rupture, site_col_scenario)
#        print rupture_gmfs.rupture_gmf
        #np.savetxt(('scenario_gmf_mmi_%s_%s.csv' % (event_name, gsim)), rupture_gmfs.rupture_gmf_mmi, delimiter =',')
        #np.savetxt(('scenario_gmf_%s_%s.csv' % (event_name, gsim)), rupture_gmfs.rupture_gmf, delimiter =',')
        # Save best fit rupture parameters
        rupture_filename = os.path.join(output_dir, ('rupture_scenario_%s_%s' % (event_name, gsim)))
        scenario_loc = np.vstack([sites_data_pts[:,0], sites_data_pts[:,1], rupture_gmfs.rupture_gmf]).transpose()
        scenario_loc_mmi = np.vstack([sites_data_pts[:,0], sites_data_pts[:,1], rupture_gmfs.rupture_gmf_mmi]).transpose()
        np.savetxt((os.path.join(output_dir, ('scenario_gmf_loc_rsa1p0_%s_%s.csv' % (event_name, gsim)))), scenario_loc, delimiter = ',')
        np.savetxt((os.path.join(output_dir, ('scenario_gmf_loc_mmi_%s_%s.csv' % (event_name, gsim)))), scenario_loc_mmi, delimiter = ',')
