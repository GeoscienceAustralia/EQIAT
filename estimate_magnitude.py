"""Use ground motion simulations to find the best fitting earthquake
source model based on historical intensity data

Jonathan Griffin
Geoscience Australia, December 2016
"""

import os, sys
import numpy as np

from gmf_calculator import get_pt_sources, RuptureGmf
from get_site_model import get_site_collection
from openquake.hazardlib.site import SiteCollection
from openquake.commonlib import nrml
from openquake.hazardlib.gsim.chiou_youngs_2008 import ChiouYoungs2008
from openquake.hazardlib.geo.point import Point
# Area source model is used to define the parameter space to be searched
# Note this should be made more generic to use faults as well
area_source_file = 'source_model_1840.xml'
#paramfile = 'params.txt'
gsim = ChiouYoungs2008()
trt = 'Active'
site_file = '1840HMMI.txt'
sites_data = np.genfromtxt(site_file)
site_points = []
for i in range(len(sites_data[:,0])):
    site_pt = Point(sites_data[i,0], sites_data[i,1])
    site_points.append(site_pt)
##[site_params]
##reference_vs30_type = 'measured'
##reference_vs30_value = 760.0
##reference_depth_to_2pt5km_per_sec = 1400.0
##reference_depth_to_1pt0km_per_sec = 300.0

#Using inferred vs30, z1.0 and z2.5 values
site_model_file = 'jawa_site_model.xml'

sitecol = get_site_collection(site_model_file, site_points) 
print sitecol

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

#Generate pt sources
pt_sources = get_pt_sources(area_source_file, discretisation = 20)
#Loop through source types
for key,item in pt_sources.iteritems():
#    print key
    # Different key for different tectonic regions
    rupture_gmfs = RuptureGmf(pt_sources[key], gsim, sitecol)#, 
                              #imts = ['0.0', '1.0'])
    rupture_gmfs.calculate()
    #    print rupture_gmfs.rupture_list
    print rupture_gmfs.gmf_list[-1]
    # Convert to MMI
    rupture_gmfs.rsa2mmi()
    print rupture_gmfs.mmi_list[-1]
    rupture_gmfs.calc_sum_squares_mmi(sites_data[:,2]) 
    #    print rupture_gmfs.sum_squares_list
    rupture_gmfs.find_best_fit()
    print 'Best magnitude', rupture_gmfs.best_rupture.mag
    print rupture_gmfs.best_rupture.hypocenter
#    print rupture_gmfs.best_rupture.surface
    print 'Best dip', rupture_gmfs.best_rupture.surface.get_dip()
    print 'Best strike', rupture_gmfs.best_rupture.surface.get_strike()
