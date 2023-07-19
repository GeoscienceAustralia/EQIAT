from get_site_model import get_site_collection
import numpy as np
from openquake.hazardlib.geo.point import Point

site_file = 'data/1699HMMI_weighted_mod.txt'
site_model_file= 'data/jawa_site_model.xml'


sites_data = np.genfromtxt(site_file)

site_points = []                           
for i in range(len(sites_data[:,0])):
    site_pt = Point(sites_data[i,0], sites_data[i,1])
    site_points.append(site_pt)

filename = None
sitecol = get_site_collection(site_model_file, site_points, None, filename)
print(sitecol)
