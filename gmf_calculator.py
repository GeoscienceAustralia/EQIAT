"""Generate ground motion fields using OpenQuake functions for locations 
where we have MMI observations. This builds a dataset by which to find 
the best fit event

Jonathan Griffin
Geoscience Australia, December 2016
"""

import os, sys
import numpy as np

from openquake.hazardlib.calc.gmf import GmfComputer
from openquake.hazardlib.nrml import SourceModelParser
from openquake.hazardlib.sourceconverter import SourceConverter, \
    area_to_point_sources, SourceGroup
from RSA2MMI import rsa2mmi9

def get_pt_sources(area_source_file, discretisation=50.):
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

def get_sources(source_model_file, discretisation=50.):
    """Calls OpenQuake parsers to read  source model
    from source_mode.xml type file, 
    and return to calculate ruptures. A more generic verions
    than that above.
    :params source_model_file:
        nrml format file of source model
    :params discretisation:
        Grid size (km) for the area source discretisation, 
        which defines the distance between resulting point
        sources.
    :returns sources
        Source for the source model
    """
    converter = SourceConverter(50, 10, width_of_mfd_bin=0.1,
                                area_source_discretization=discretisation)
    parser = SourceModelParser(converter)
    try:
        sources = parser.parse_sources(source_model_file)
    except AttributeError: # Handle version 2.1 and above
        sources = []
        groups = parser.parse_src_groups(source_model_file)
        for group in groups:
            for source in group:
                sources.append(source)
#    name = 'test_point_model'
    new_sources = {}
    for source in sources:
        #pt_sources = area_to_point_sources(source)
        #for pt in pt_sources:
        #    pt.source_id = pt.source_id.replace(':','')
        #    pt.name = pt.name.replace(':','_')
            try:
                new_sources[source.tectonic_region_type].append(source)
            except KeyError:
                new_sources[source.tectonic_region_type] = [source]
    return new_sources

class RuptureGmf(object):
    """Class for storing ruptures and associated
    ground motion fields for later analysis
    """

    def __init__(self, sources, gsim, sitecol, imts = ['SA(1.0)']):
        """
        :params pt_sources:
            Point source objects derived from original area source model
        :params gsim:
            GSIM instance (i.e. subclass of openquake.hazardlib.gsim.base.GMPE)
        """
        self.sources = sources
        self.gsim  = gsim
        self.imts = imts
        self.sitecol = sitecol
        self.rupture_list = [] # list for storing ruptures
        self.gmf_list = [] # list for storing associated gmfs
        self.mmi_list = []
    
    def calculate_from_pts(self):
        """Generates ruptures for each pt source and calculates ground motion
        field.
        :returns gmfs:
            Set of ruptures and associated parameters for ground motion
            calculations  
        """
        for pt in self.sources:
            #        rupture_mags = []
            #        rupture_hypocenter = []
            ruptures = pt.iter_ruptures()
            for rupture in ruptures:
                computer = GmfComputer(rupture, self.sitecol,
                                       self.imts, [self.gsim],
                                       truncation_level=0)
                gmf = computer.compute(self.gsim, 1)
                gmf = gmf.flatten()
                self.rupture_list.append(rupture)
                self.gmf_list.append(gmf)

    def calculate(self):
        """Generates ruptures for each source and calculates ground motion
        field.
        :returns gmfs:
            Set of ruptures and associated parameters for ground motion
            calculations  
        """
        for source in self.sources:
            #        rupture_mags = []
            #        rupture_hypocenter = []
            ruptures = source.iter_ruptures()
            for rupture in ruptures:
               # print 'Calculating rupture', rupture.hypocenter
                computer = GmfComputer(rupture, self.sitecol,
                                       self.imts, [self.gsim],
                                       truncation_level=0)
                gmf = computer.compute(self.gsim, 1)
                gmf = gmf.flatten()
                self.rupture_list.append(rupture)
                self.gmf_list.append(gmf)

    def rsa2mmi(self):
        """Convert ground motion fields to MMI intensity
        """
        for gmf in self.gmf_list:
            mmi = rsa2mmi9(gmf, period = 1.0)
            self.mmi_list.append(mmi)

    def calc_sum_squares_mmi(self, mmi_obs):
        """Calculates sum of squares for each rupture gmf compared 
        with historical observations
        """
        self.sum_squares_list = []
        for mmi in self.mmi_list:
            sum_squares = np.sum((mmi - mmi_obs)**2)
            self.sum_squares_list.append(sum_squares)

    def rmse(self, mmi_obs):
        """Calculates root-mean-square error of each rupture gmt
        compared with historical observations
        """
        try:
            self.sum_squares_list
        except AttributeError:
            self.calc_sum_squares_mmi(mmi_obs)
        self.rmse = np.sqrt(np.array(self.sum_squares_list)/float(len(mmi_obs)))

    def calc_sum_squares_mmi_weighted(self, mmi_obs):
        """Calculates sum of squares for each rupture gmf compared 
        with historical observations
        """
        self.sum_squares_list = []
        weights = np.where(mmi_obs < 5, 2, 1) # Increase weight for low MMI events
        for mmi in self.mmi_list:
            sum_squares = np.sum(np.dot(weights,(mmi - mmi_obs))**2)/(np.sum(weights**2))
            self.sum_squares_list.append(sum_squares)

    def find_best_fit(self):
        """Find rupture with minimm sum of squares
        """
        index = np.argmin(self.rmse)
        self.best_rupture = self.rupture_list[index]
        self.min_rmse = self.rmse[index]
