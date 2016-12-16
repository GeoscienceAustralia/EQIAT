"""Generate ground motion fields using OpenQuake functions for locations 
where we have MMI observations. This builds a dataset by which to find 
the best fit event

Jonathan Griffin
Geoscience Australia, December 2016
"""

import os, sys
import numpy as np

from openquake.hazardlib.calc.gmf import GmfComputer
from openquake.commonlib.source import SourceModelParser
from openquake.commonlib.sourceconverter import SourceConverter, \
    area_to_point_sources, SourceGroup

paramfile = 'params.txt'

def get_pt_sources(area_source_file, discretisation=200.):
    """Calls OpenQuake parsers to read area source model
    from source_mode.xml type file, convert to point sources
    and return to calculate ruptures.
    :params area_source_file:
        nrml format file of the area source
    :params discretisation:
        Grid size (km) for the area source discretisation, 
        which defines the distance between resulting point
        sources.
    :returns new_pt_sources
        Point source for the area source model
    """
    converter = SourceConverter(50, 10, width_of_mfd_bin=0.1,
                                area_source_discretization=discretisation)
    parser = SourceModelParser(converter)
    print [method for method in dir(parser)]# if callable(getattr(parser, method))]
    try:
        sources = parser.parse_sources(area_source_file)
    except AttributeError: # Handle version 2.1 and above
        sources = []
        groups = parser.parse_src_groups(area_source_file)
        for group in groups:
            for source in group:
                sources.append(source)
    name = 'test_point_model'
    new_pt_sources = {}
    for source in sources:
        pt_sources = area_to_point_sources(source)
        for pt in pt_sources:
            pt.source_id = pt.source_id.replace(':','')
            pt.name = pt.name.replace(':','_')
            try:
                new_pt_sources[pt.tectonic_region_type].append(pt)
            except KeyError:
                new_pt_sources[pt.tectonic_region_type] = [pt]
    return new_pt_sources

class rupture_gmf(object):
    """Class for storing ruptures and associated
    ground motion fields for later analysis

def gmf_calculate(pt_sources, gsim):
    """Generates ruptures for each pt source and calculates ground motion
    field.
    :params pt_sources:
        Point source objects derived from original area source model
    :params gsim:
        GSIM instance (i.e. subclass of openquake.hazardlib.gsim.base.GMPE)
    :returns gmfs:
        Set of ruptures and associated parameters for ground motion
        calculations
    """
    
    for pt in pt_sources:
#        rupture_mags = []
#        rupture_hypocenter = []
        ruptures = pt.iter_ruptures()
        for rupture in ruptures:
#            rupture_mags.append(rupture.mag)
#            rupture_hypocenter.append(rupture.hypocenter)
            computer = GmfComputer(
                rupture, self.sitecol, ['0.0', '1.0'], [gsim])
            gmf = computer.compute(gsim, 1)
            
 #       rupture_mags = np.array(rupture_mags).flatten()
                # make the same length as the corners
#        rupture_mags = np.repeat(rupture_mags, 4)
#        rupture_lons = np.array(rupture_lons).flatten()

