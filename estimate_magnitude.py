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
#from openquake.hazardlib.contexts import ContextMaker
from openquake.hazardlib.gsim.campbell_bozorgnia_2008 import CampbellBozorgnia2008
from openquake.hazardlib.gsim.campbell_bozorgnia_2014 import CampbellBozorgnia2014
from openquake.hazardlib.gsim.abrahamson_2015 import AbrahamsonEtAl2015SSlab, AbrahamsonEtAl2015SInter
from openquake.hazardlib.gsim.chiou_youngs_2008 import ChiouYoungs2008
from openquake.hazardlib.gsim.chiou_youngs_2014 import ChiouYoungs2014
from openquake.hazardlib.gsim.boore_atkinson_2008 import BooreAtkinson2008
from openquake.hazardlib.gsim.boore_2014 import BooreEtAl2014
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
try:
    limits_filename = f_in.readline().rstrip().split(',')[1]
except:
    limits_filename = None
if event_name == '1852Banda_area' or event_name == '1852Banda_domain_ryan_mmi' or \
        event_name == '1852Banda_domain_FH_mmi' or event_name == '1852Banda_exclude_20min_ryan_mmi' or\
        event_name == '1852Banda_exclude_20min_FH_mmi' or event_name == '1852Banda_exclude_15min_FH_mmi' or \
        event_name == '1852Banda_exclude_15min_ryan_mmi':
    disc = 10
else:
    disc = 5
# Area or fault source model is used to define the parameter space to be searched
if trt == 'Subduction Intraslab':
    gsim_list = [ZhaoEtAl2006SSlab(), AtkinsonBoore2003SSlab(), AtkinsonBoore2003SSlabCascadia(), AbrahamsonEtAl2015SSlab()]
if trt == 'Active':
#    gsim_list = [ChiouYoungs2008(), ChiouYoungs2014(), BooreAtkinson2008(), BooreEtAl2014(), CampbellBozorgnia2008(), CampbellBozorgnia2014() ]
    gsim_list = [ChiouYoungs2014(), BooreEtAl2014(), CampbellBozorgnia2014()]
if trt == 'Subduction Interface':
    gsim_list = [AbrahamsonEtAl2015SInter(), ZhaoEtAl2006SInter(), YoungsEtAl1997SInter(), AtkinsonBoore2003SInter()]
if trt == 'Mixed':
    gsim_list = [AbrahamsonEtAl2015SInter(), BooreEtAl2014(), ChiouYoungs2014(), ZhaoEtAl2006SInter()]
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
#Build sites for observations
sites_data = np.genfromtxt(site_file)
sitecol = build_site_col(sites_data, site_model_file)
print('Sitecol', sitecol, type(sitecol))
# Yield to list from generator
##def yield_sites():
##    yield from sitecol
##sitecol = list(yield_sites())
#print('Sitecol', sitecol, type(sitecol))
# Build sites for scenarios (or use pre-calculated site file
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
print('sitecol', sitecol, type(sitecol))

output_dir = 'outputs'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
#Generate pt sources
print('Generating source model')
pt_sources = get_sources(area_source_file, discretisation = disc)
#Loop through source types
for key,item in pt_sources.items():
    # Different key for different tectonic regions
    for gsim in gsim_list:
        print('Running for %s' % gsim)
        rupture_gmfs = RuptureGmf(pt_sources[key], gsim, sitecol, trt)
        rupture_gmfs.calculate()
        print(rupture_gmfs.gmf_list[-1])
        # Convert to MMI
        rupture_gmfs.rsa2mmi()
        print(rupture_gmfs.mmi_list[-1])
        rupture_gmfs.calc_rmse(sites_data[:,2], weights = sites_data[:,3])
        rupture_gmfs.find_best_fit()
        print('Results using %s GMM' % gsim)
        print('Best magnitude', rupture_gmfs.best_rupture.mag)
        print(rupture_gmfs.best_rupture.hypocenter)
        # Print rupture_gmfs.best_rupture.surface
        print('Best dip', rupture_gmfs.best_rupture.surface.get_dip())
        print('Best strike', rupture_gmfs.best_rupture.surface.get_strike())
        # Write best-fit rupture to shapefile                                                                                      
        output_shp_filename = 'rupture_scenario_%s_%s_Mw%.2f_%.3f_%.3f_%.2fkm.shp' % (event_name, gsim,
                                                                                      rupture_gmfs.best_rupture.mag,
                                                                                      rupture_gmfs.best_rupture.hypocenter.longitude,
                                                                                      rupture_gmfs.best_rupture.hypocenter.latitude,
                                                                                      rupture_gmfs.best_rupture.hypocenter.depth)
        output_shp = os.path.join(output_dir, output_shp_filename)
        #print type(rupture_gmfs.best_rupture.surface)
        if isinstance(rupture_gmfs.best_rupture.surface, SimpleFaultSurface) or \
                isinstance(rupture_gmfs.best_rupture.surface, ComplexFaultSurface):
            # Dealing with simple faults
            alllons = rupture_gmfs.best_rupture.surface.mesh.lons
            lons = alllons 
            alllats = rupture_gmfs.best_rupture.surface.mesh.lats
            lats = alllats 
            alldepths = rupture_gmfs.best_rupture.surface.mesh.depths
            depths = alldepths 

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
                print('Could not make shapefile, writing to text')
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
            print('Corner lons', rupture_gmfs.best_rupture.surface.corner_lons)
            print('Corner lats', rupture_gmfs.best_rupture.surface.corner_lats)
            print('Corner depths', rupture_gmfs.best_rupture.surface.corner_depths)
            fault2shp(rupture_gmfs.best_rupture.surface.corner_lons, 
                      rupture_gmfs.best_rupture.surface.corner_lats,
                      output_shp,
                      rupture_gmfs.best_rupture.surface.corner_depths)
        print('RMSE', rupture_gmfs.min_rmse)
        # Calculate uncertainty model
        parameter_llh_filename = os.path.join(output_dir, event_name + ('_%s_parameter_llh.csv' % gsim))
        parameter_llh_filename = parameter_llh_filename.replace('()', '')
        rupture_gmfs.uncertainty_model(filename=parameter_llh_filename)
        try:
            print('Magnitude range is %.2f - %.2f' % (rupture_gmfs.min_mag, rupture_gmfs.max_mag))
            print('Longitude range is %.2f - %.2f' % (rupture_gmfs.min_lon,rupture_gmfs.max_lon))
            print('Latitude range is %.2f - %.2f' % (rupture_gmfs.min_lat,rupture_gmfs.max_lat))
            print('Depth range is %.2f - %.2f' % (rupture_gmfs.min_depth,rupture_gmfs.max_depth))
            print('Strike range is %.2f - %.2f' % (rupture_gmfs.min_strike,rupture_gmfs.max_strike))
            print('Dip range is %.2f - %.2f' % (rupture_gmfs.min_dip,rupture_gmfs.max_dip))
        except:
            # Not enough data to calcualte uncertainties
            pass
        # Get parameter pdfs
        fig_comment = 'figures/' + event_name
        try:
            rupture_gmfs.parameter_pdf(fig_comment=fig_comment,limits_filename=limits_filename)
        # Get slices
            fig_comment = 'rmse_png/' + event_name
            lons,lats,rmse = rupture_gmfs.uncertainty_slice2D('longitude', 'latitude', 'mag', 
                                                              rupture_gmfs.best_rupture.mag,
                                                              fig_comment=fig_comment,
                                                              limits_filename=limits_filename)
            unc_array = np.vstack([lons,lats,rmse])
            np.savetxt('rmse_csv/rmse_%s_%s.csv' % (event_name, gsim), unc_array.T, delimiter=',')
        except:
            # Not enough data to calcualte uncertainties                                                                             
            pass
        try:
            min_mag, max_mag =  rupture_gmfs.uncertainty_slice1D('mag', 'longitude', 'latitude',
                                                                 rupture_gmfs.best_rupture.hypocenter.longitude,
                                                                 rupture_gmfs.best_rupture.hypocenter.latitude)
            print('95 percent confidence range on magnitude at best-fit location is %.2f - %.2f' % (
                min_mag, max_mag))
        except AttributeError:
            print('Not enough data to calculate magnitude uncertainty on best-fit location')
        # Calculate scenario for best rupture
        rupture_gmfs.calculate_from_rupture(rupture_gmfs.best_rupture, site_col_scenario)
        # Save best fit rupture parameters
        rupture_filename = os.path.join(output_dir, ('rupture_scenario_%s_%s' % (event_name, gsim)))
        scenario_loc = np.vstack([sites_data_pts[:,0], sites_data_pts[:,1], rupture_gmfs.rupture_gmf]).transpose()
        scenario_loc_mmi = np.vstack([sites_data_pts[:,0], sites_data_pts[:,1], rupture_gmfs.rupture_gmf_mmi]).transpose()
        np.savetxt((os.path.join(output_dir, ('scenario_gmf_loc_rsa1p0_%s_%s.csv' % (event_name, gsim)))), scenario_loc, delimiter = ',')
        np.savetxt((os.path.join(output_dir, ('scenario_gmf_loc_mmi_%s_%s.csv' % (event_name, gsim)))), scenario_loc_mmi, delimiter = ',')
