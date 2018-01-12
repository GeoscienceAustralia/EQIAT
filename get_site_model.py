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
from openquake.risklib import valid
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
    for node in nrml.read(site_model_file).siteModel:
        yield valid.site_param(**node.attrib)

def read_site_col(site_model_file):
    """Directly read a site colection from nrml
    """
    sitecol = []
    for param in sorted(get_site_model(site_model_file)):
        pt = geo.Point(param.lon, param.lat, 0.)
        sitecol.append(site.Site(
                pt, param.vs30, param.measured,
                param.z1pt0, param.z2pt5, param.backarc))
    return site.SiteCollection(sitecol)

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
        site_model_params = geo.geodetic.GeographicObjects(
            get_site_model(site_model_file))
    sitecol = []
    for pt in sites:
        # NB: the mesh, when read from the datastore, is a 32 bit array;
        # however, the underlying C library expects 64 bit floats, thus
        # we have to cast float(pt.longitude), float(pt.latitude);
        # we should change the geodetic speedups instead
        param, dist = site_model_params.\
                      get_closest(float(pt.longitude), float(pt.latitude))
        if dist >= MAX_SITE_MODEL_DISTANCE:
            #logging.warn('The site parameter associated to %s came from a '
            #             'distance of %d km!' % (pt, dist))
            print 'WARNING:The site parameter associated to %s came from a ' \
                'distance of %d km!' % (pt, dist)
        sitecol.append(
            site.Site(pt, param.vs30, param.measured,
                      param.z1pt0, param.z2pt5, param.backarc))
    if filename is not None:
        name = filename.split('/')[-1][:-4]
        site_nodes = list(map(obj_to_node, sitecol))#sorted(sitecol)))
        #print site_nodes
        site_model = Node("siteModel", nodes=site_nodes)
        #site_model = site_nodes
        #print site_model
        with open(filename, 'wb') as f:
            nrml.write(site_model, f, '%s', xmlns = NAMESPACE)
        #    nrml.write([site_model], f, '%s', xmlns = NAMESPACE)
    return site.SiteCollection(sitecol)

