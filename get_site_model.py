"""My hacked version from openquake.commonlib.readinput

Jonathan Griffin
Geoscienec Australia
December 2016
"""

from openquake.hazardlib import nrml
from openquake.hazardlib import geo, site, imt
from openquake.risklib import valid
MAX_SITE_MODEL_DISTANCE = 5  # km, given by Graeme Weatherill

def get_site_model(site_model_file):
    """
    Convert the NRML file into an iterator over 6-tuple of the form
    (z1pt0, z2pt5, measured, vs30, lon, lat)
    :param oqparam:
        an :class:`openquake.commonlib.oqvalidation.OqParam` instance
    """
    for node in nrml.read(site_model_file).siteModel:
        yield valid.site_param(**node.attrib)


def get_site_collection(site_model_file, sites, site_model_params=None):
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
    return site.SiteCollection(sitecol)

