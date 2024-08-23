"""My hacked version from openquake.commonlib.readinput

Jonathan Griffin
Geoscienec Australia
December 2016
"""

from openquake.hazardlib import nrml
from openquake.hazardlib.nrml import NAMESPACE
from openquake.hazardlib.sourcewriter import obj_to_node
from openquake.baselib.node import Node
from openquake.hazardlib import geo, site, imt
from openquake.hazardlib.site import SiteCollection, Site
#from openquake.risklib import valid
from openquake.hazardlib import valid
import numpy as np

MAX_SITE_MODEL_DISTANCE = 5  # km, given by Graeme Weatherill

@obj_to_node.add('Site')
def build_site(site):
    return Node('site',
                {'lon': site.location.longitude, 'lat': site.location.latitude,
                 'vs30': site.vs30, 'vs30Type': site.vs30measured,
                 'z1pt0': site.z1pt0, 'z2pt5': site.z2pt5,
                 'backarc': site.backarc})

@obj_to_node.add('SiteModel')
def build_site_model(siteModel):
    site_nodes = []
    for site in siteModel:
        site_nodes.append(build_site(site))
    return Node('siteModel', nodes=site_nodes)

#@obj_to_node.add('siteModel')
#def build_site_model(siteModel):
#    site_nodes = []
#    for site in siteModel:
#        site_nodes.append(build_site(site))
#    return Node('siteModel', nodes=site_nodes)    

def build_site(site):
    return Node('site',
                {'lon': site.longitude, 'lat': site.latitude,
                 'vs30': site.vs30, 'vs30Type': site.measured,
                 'z1pt0': site.z1pt0, 'z2pt5': param.z2pt5,
                 'backarc': site.backarc})

def get_site_model(site_model_file):
    """
    Convert the NRML file into an iterator over 6-tuple of the form
    (z1pt0, z2pt5, measured, vs30, lon, lat)
    :param oqparam:
        an :class:`openquake.commonlib.oqvalidation.OqParam` instance
    """
#    print(site_model_file)
    try:
        for node in nrml.read(site_model_file).siteModel:
#        print(node, type(node))
#        print(node.attrib)
#        yield valid.site_param(**node.attrib) # Orginal
#        new = valid.site_param(node.attrib)
#        print(new)
#        print(valid.site_param(node.attrib))
            yield valid.site_param(node.attrib)
    except AttributeError:
        for node in nrml.read(site_model_file):
            yield valid.site_param(node.attrib)
            

def read_site_col(site_model_file):
    """Directly read a site colection from nrml
    """
    sitecol = []
    for param in get_site_model(site_model_file):
#        print(param)
#        pt = geo.Point(param.lon, param.lat, 0.)
        pt = geo.Point(param['lon'], param['lat'], 0.)
#        print(pt)
#        sitecol.append(site.Site(
#                pt, param.vs30, param.measured,
#                param.z1pt0, param.z2pt5, param.backarc))
        sitecol.append(site.Site(
            pt, float(param['vs30']),
            float(param['z1pt0']), float(param['z2pt5']),
            vs30measured = param['vs30measured']))#, param['backarc']))
    return site.SiteCollection(sitecol)

def read_site_col_csv(site_model_file):
    """Read site collection from a csv file"""
    sitecol = []
    data = np.genfromtxt(site_model_file, delimiter=',', skip_header=1)
    for i in range(len(data)):
        pt = geo.Point(float(data[i,0]), float(data[i,1]))
        sitecol.append(site.Site(pt, float(data[i,2]),
                                 float(data[i,3]), float(data[i,4]),
                                 vs30measured = float(data[i,5]),
                                 backarc = float(data[i,6])))
    return site.SiteCollection(sitecol)                                 
                                           
        
def yield_site_mod(site_mod):
    yield from site_mod

def get_site_collection(site_model_file, sites, site_model_params=None,
                        filename=None):
    """
    Returns a SiteCollection instance by looking at the points and the
    site model defined by the configuration parameters.
    :param site_model_file:
        path to the site_model nrml file
    :param sites:
        Locations of hazard points that 
    :param site_model_params:
        object with a method .get_closest returning the closest site
        model parameters
    :param filename:
        path to output nrml file where new site model will be written,
        if this parameter is specified.
    """
    
    if site_model_params is None:
        # read the parameters directly from their file
        if site_model_file.endswith('.csv'):
            site_mod_read = read_site_col_csv(site_model_file)
            print('filename', filename)
#            if filename is not None:
#                name = filename.split('/')[-1][:-4]
#                site_nodes = list(map(obj_to_node, site_mod_read))#sorted(sitecol)))   
#                site_model = Node("siteModel", nodes=site_nodes)
##                print(site_model)
#               with open(filename, 'wb') as f:
#                    nrml.write(site_model, f, '%s', xmlns = NAMESPACE)
#            return site_mod_read
        else:
            site_mod_read = read_site_col(site_model_file)
        #### Try associate function
##        target_sites = SiteCollection(sites)
##        print('target_sites', target_sites)
##        print(sites, type(sites))
##        assoc_dist = 100.0
 ##       site_model = site_mod_read.assoc(sites, assoc_dist, ignore=())
##        print('site model', site_model)
##        sys.exit()
        ####
#        site_mod = get_site_model_v2(site_model_file)
#        site_mod = list(yield_site_mod(site_mod))
#        print('Here')
#        print(site_mod)
        site_model_params = geo.utils._GeographicObjects(
            site_mod_read)
#        print(site_model_params)
#        site_model_params = geo.utils._GeographicObjects(
#            site_mod)
#    site_model_params = site_mod_read #Test!
    sitecol = []
    for pt in sites:
        # NB: the mesh, when read from the datastore, is a 32 bit array;
        # however, the underlying C library expects 64 bit floats, thus
        # we have to cast float(pt.longitude), float(pt.latitude);
        # we should change the geodetic speedups instead
        param, dist = site_model_params.\
                      get_closest(float(pt.longitude), float(pt.latitude))
        print('param',param)
        if dist >= MAX_SITE_MODEL_DISTANCE:
            #logging.warn('The site parameter associated to %s came from a '
            #             'distance of %d km!' % (pt, dist))
            print('WARNING:The site parameter associated to %s came from a ' \
                'distance of %d km!' % (pt, dist))
#        sitecol.append(
#            site.Site(pt, param.vs30, param.measured,
#                      param.z1pt0, param.z2pt5, param.backarc))
#        print(param)
#        sitecol.append(                                                                                                                 
#            site.Site(pt, param[4],                                                                                                       
#                      param[6], param[7], vs30measured=1, backarc = 0)) ## FIXME -currentl backarc is hardcoded
        sitecol.append(
            site.Site(pt, param[5],
                      param[7], param[8], vs30measured=1, backarc=0))
    if filename is not None:
        name = filename.split('/')[-1][:-4]
        site_nodes = list(map(obj_to_node, sitecol))#sorted(sitecol)))
        site_model = Node("siteModel", nodes=site_nodes)
        with open(filename, 'wb') as f:
            nrml.write(site_model, f, '%s', xmlns = NAMESPACE)
    return site.SiteCollection(sitecol)


