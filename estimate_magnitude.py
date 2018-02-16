"""Use ground motion simulations to find the best fitting earthquake
source model based on historical intensity data

Jonathan Griffin
Geoscience Australia, December 2016
"""

import os, sys
import numpy as np
import argparse

from gmf_calculator import get_pt_sources, get_sources, RuptureGmf
from get_site_model import get_site_collection, read_site_col
from write_fault_shp import fault2shp
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
from openquake.hazardlib.geo.surface.simple_fault import SimpleFaultSurface
from openquake.hazardlib.geo.surface.complex_fault import ComplexFaultSurface

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
#gsim_list = [ZhaoEtAl2006SSlab(), AtkinsonBoore2003SSlab()]
#trt = 'Active'
if trt == 'Subduction Intraslab':
    gsim_list = [ZhaoEtAl2006SSlab(), AtkinsonBoore2003SSlab(), AtkinsonBoore2003SSlabCascadia()]
if trt == 'Active':
    gsim_list = [ChiouYoungs2008(), BooreAtkinson2008(), CampbellBozorgnia2008()]
if trt == 'Subduction Interface':
    gsim_list = [YoungsEtAl1997SInter(), AtkinsonBoore2003SInter(), ZhaoEtAl2006SInter()]

def build_site_col(sites_data, site_model_file, filename=None):
    """Interpolate vs30 values to sites of interest
    """
#    sites_data = np.genfromtxt(site_file)
    site_points = []
    for i in range(len(sites_data[:,0])):
        site_pt = Point(sites_data[i,0], sites_data[i,1])
        site_points.append(site_pt)
    sitecol = get_site_collection(site_model_file, site_points, None,filename)
    return sitecol
    

# Background site class model
#Using inferred vs30, z1.0 and z2.5 values                                                               
#site_model_file = 'data/jawa_site_model.xml'

#Build sites for observations
#site_file = 'data/1699HMMI_weighted_mod.txt'
sites_data = np.genfromtxt(site_file)
sitecol = build_site_col(sites_data, site_model_file)
# Build sites for scenarios (or use pre-calculated site file
#site_file_pts = 'data/jawa_sites_thin.csv'
if site_file_pts.endswith('.csv'):
    sites_data_pts = np.genfromtxt(site_file_pts, delimiter=',')
    site_col_scenario = build_site_col(sites_data_pts, site_model_file)
elif site_file_pts.endswith('.xml'):
    site_col_scenario = read_site_col(site_file_pts)
    sites_data_pts = np.vstack([site_col_scenario.lons, site_col_scenario.lats]).transpose()
else:
    msg = 'Invalid site model file %s ' % site_file_pts
    raise ValueError(msg)


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
# write best-fit rupture to shapefile                                                                                      
        output_shp_filename = 'rupture_scenario_%s_%s_Mw%.2f_%.3f_%.3f_%.2fkm.shp' % (event_name, gsim,
                                                                                      rupture_gmfs.best_rupture.mag,
                                                                                      rupture_gmfs.best_rupture.hypocenter.longitude,
                                                                                      rupture_gmfs.best_rupture.hypocenter.latitude,
                                                                                      rupture_gmfs.best_rupture.hypocenter.depth)
        output_shp = os.path.join(output_dir, output_shp_filename)
        print type(rupture_gmfs.best_rupture.surface)
        if isinstance(rupture_gmfs.best_rupture.surface, SimpleFaultSurface) or \
                isinstance(rupture_gmfs.best_rupture.surface, ComplexFaultSurface):
            # Dealing with simple faults
            #print rupture_gmfs.best_rupture.surface.mesh.__dict__
            #print rupture_gmfs.best_rupture.surface.__dict__
            #print rupture_gmfs.best_rupture.__dict__
            #print type(rupture_gmfs.best_rupture.surface.mesh)
            #print type(rupture_gmfs.best_rupture.surface)
            #print type(rupture_gmfs.best_rupture)
            alllons = rupture_gmfs.best_rupture.surface.mesh.lons
            #print alllons
            #print type(alllons)
            lons = alllons #np.concatenate(alllons[0], alllons[-1][::-1])
            alllats = rupture_gmfs.best_rupture.surface.mesh.lats
            #print alllats
            lats = alllats #np.concatenate(alllats[0], alllats[-1][::-1])
            alldepths = rupture_gmfs.best_rupture.surface.mesh.depths
            depths = alldepths #np.concatenate(alldepths[0], alldepths[-1][::-1])
            #print lons, lats, depths
            #lons,lats = rupture_gmfs.best_rupture.surface.surface_projection_from_fault_data(rupture_gmfs.best_rupture.surface.fault_trace,
            #                                                                                 rupture_gmfs.best_rupture.surface.upper_seismogenic_depth,
            #                                                                                 rupture_gmfs.best_rupture.surface.lower_seismogenic_depth,
            #                                                                                 rupture_gmfs.best_rupture.surface.dip)
            #print lons,lats
            try:
                fault2shp(lons,
                          lats,
                          output_shp,
                          depths)
            except:
                print 'Could not make shapefile, writing to text'
                fault_mesh = np.vstack([lons,lats,depths])
                mesh_filename = 'rupture_mesh_%s_%s_Mw%.2f_%.3f_%.3f_%.2fkm.txt' % (
                    event_name, gsim,
                    rupture_gmfs.best_rupture.mag,
                    rupture_gmfs.best_rupture.hypocenter.longitude,
                    rupture_gmfs.best_rupture.hypocenter.latitude,
                    rupture_gmfs.best_rupture.hypocenter.depth)
                mesh_filename = os.path.join(output_dir, mesh_filename)
                np.savetxt(mesh_filename, fault_mesh, delimiter=',')
        else:
            print 'Corner lons', rupture_gmfs.best_rupture.surface.corner_lons
            print 'Corner lats', rupture_gmfs.best_rupture.surface.corner_lats
            print 'Corner depths', rupture_gmfs.best_rupture.surface.corner_depths
            fault2shp(rupture_gmfs.best_rupture.surface.corner_lons, 
                      rupture_gmfs.best_rupture.surface.corner_lats,
                      output_shp,
                      rupture_gmfs.best_rupture.surface.corner_depths)
        print 'RMSE', rupture_gmfs.min_rmse
        # Calculate uncertainty model
        rupture_gmfs.uncertainty_model()
        try:
            print 'Magnitude range is %.2f - %.2f' % (rupture_gmfs.min_mag, rupture_gmfs.max_mag)
            print 'Longitude range is %.2f - %.2f' % (rupture_gmfs.min_lon,rupture_gmfs.max_lon)
            print 'Latitude range is %.2f - %.2f' % (rupture_gmfs.min_lat,rupture_gmfs.max_lat)
            print 'Depth range is %.2f - %.2f' % (rupture_gmfs.min_depth,rupture_gmfs.max_depth)
            print 'Strike range is %.2f - %.2f' % (rupture_gmfs.min_strike,rupture_gmfs.max_strike)
            print 'Dip range is %.2f - %.2f' % (rupture_gmfs.min_dip,rupture_gmfs.max_dip)
        except:
            # Not enough data to calcualte uncertainties
            pass
        # Get slices
        fig_comment = event_name
        lons,lats,rmse = rupture_gmfs.uncertainty_slice('longitude', 'latitude', 'mag', 
                                                        rupture_gmfs.best_rupture.mag,
                                                        fig_comment=fig_comment)
        print lons, lats, rmse
        unc_array = np.vstack([lons,lats,rmse])
        np.savetxt('rmse_%s_%s.csv' % (event_name, gsim), unc_array.T, delimiter=',')
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
