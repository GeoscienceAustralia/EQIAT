"""Generate ground motion fields using OpenQuake functions for locations 
where we have MMI observations. This builds a dataset by which to find 
the best fit event

Jonathan Griffin
Geoscience Australia, December 2016
"""

import os, sys
import numpy as np
from scipy.stats import norm, chi2
from scipy import interpolate
from mpl_toolkits.basemap import Basemap, maskoceans
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Polygon
from collections import defaultdict
from openquake.hazardlib.calc.gmf import GmfComputer
from openquake.hazardlib.nrml import SourceModelParser
from openquake.hazardlib.sourceconverter import SourceConverter, \
    area_to_point_sources, SourceGroup
from RSA2MMI import rsa2mmi8


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
        :params sources:
            Source objects derived from original area source model
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
                #print type(rupture)
                #print 'Calculating rupture', rupture.hypocenter
                computer = GmfComputer(rupture, self.sitecol,
                                       self.imts, [self.gsim],
                                       truncation_level=0)
                gmf = computer.compute(self.gsim, 1)
                gmf = gmf.flatten()
                self.rupture_list.append(rupture)
                self.gmf_list.append(gmf)

    def calculate_from_rupture(self, rupture, rup_sitecol=None):
        """Method to generate scenario ground motion
        for a specific rupture
        """
        if rup_sitecol is None:
            rup_sitecol = self.sitecol
        computer = GmfComputer(rupture, rup_sitecol,
                               self.imts, [self.gsim],
                               truncation_level=0)
        gmf = computer.compute(self.gsim, 1)
        gmf = gmf.flatten()
        print 'gmf', gmf
        self.rupture_scenario = rupture
        self.rupture_gmf = gmf
        self.rupture_gmf_mmi = rsa2mmi8(gmf, period = 1.0)
        

    def rsa2mmi(self):
        """Convert ground motion fields to MMI intensity
        """
        for gmf in self.gmf_list:
            mmi = rsa2mmi8(gmf, period = 1.0)
            self.mmi_list.append(mmi)

    def calc_sum_squares_mmi(self, mmi_obs):
        """Calculates sum of squares for each rupture gmf compared 
        with historical observations
        """
        self.sum_squares_list = []
        for mmi in self.mmi_list:
            sum_squares = np.sum((mmi - mmi_obs)**2)
            self.sum_squares_list.append(sum_squares)

    def calc_rmse(self, mmi_obs, weights = None):
        """Calculates root-mean-square error of each rupture gmf
        compared with historical observations
        """
        self.mmi_obs = mmi_obs
        try:
            self.sum_squares_list
        except AttributeError:
            if weights is not None:
                self.calc_sum_squares_mmi_weighted(self.mmi_obs, weights)
                self.rmse = np.sqrt(np.array(self.sum_squares_list))
            else:
                self.calc_sum_squares_mmi(self.mmi_obs)
                self.rmse = np.sqrt(np.array(self.sum_squares_list)/float(len(self.mmi_obs)))

    def calc_sum_squares_mmi_weighted(self, mmi_obs, weights):
        """Calculates sum of squares for each rupture gmf compared 
        with historical observations
        """
        self.sum_squares_list = []
#        weights = np.where(mmi_obs < 5, 2, 1) # Increase weight for low MMI events
        for mmi in self.mmi_list:
            #sum_squares = np.sum(np.dot(weights,(mmi - mmi_obs))**2)/(np.sum(weights**2))
#            print weights
            weights = np.array(weights)
            weights = weights*(1./sum(weights)) # Normalise weights to sum to 1
#            print weights
#            print sum(weights)
            sum_squares = np.dot(weights,(mmi - mmi_obs)**2)
#            sum_squares = np.sum(np.dot(weights,(mmi - mmi_obs))**2)
            self.sum_squares_list.append(sum_squares)

    def find_best_fit(self):
        """Find rupture with minimm sum of squares
        """
        index = np.argmin(self.rmse)
        self.best_rupture = self.rupture_list[index]
        self.min_rmse = self.rmse[index]

    def rupture_params_2_array(self):
        """Create an array of all rupture parameters for later slicing                                                                      
        """
        # Create an array of the full rupture parameter space                                                                               
        mags = np.array([rup.mag for rup in self.rupture_list])
        lons = np.array([rup.hypocenter.longitude for rup in self.rupture_list])
        lats = np.array([rup.hypocenter.latitude for rup in self.rupture_list])
        depths = np.array([rup.hypocenter.depth for rup in self.rupture_list])
        strikes = np.array([rup.surface.get_strike() for rup in self.rupture_list])
        dips = np.array([rup.surface.get_dip() for rup in self.rupture_list])
        self.parameter_space = np.vstack([mags,lons,lats,depths,strikes,dips])
        self.parameter_dict = {'mag': 0, 'longitude': 1,
                          'latitude': 2, 'depth': 3,
                          'strike': 4, 'dip':5}

    def uncertainty_model(self, min_rmse=None, filename='param_model.csv'):
        """Estimate parameter uncertainties by 
        assuming the minimum rmse is maximum likelihood.
        We use the residuals of the best-fit parameters 
        to estimate the standard deviation of the model fit. 
        Can optionally pass min rmse,
        e.g. if you wish to apply uncertainties from another model 
        (i.e. using a different GMPE).
        """
        if len(self.mmi_obs) <= 6:
            print 'Not enough data points to calculate uncertainties'
            #indices = np.where(self.rmse < 1e24)[0]
            # Hacky estimate of sigma
            self.sigma = 0.3
            index = np.argmin(self.rmse)
            self.uncert_fun = norm(min(self.rmse),self.sigma)
            self.max_lik = np.power((1/(self.sigma*np.sqrt(2*np.pi))), len(self.mmi_obs)) * \
                np.exp((-1/2)*((self.sum_squares_list[index]/self.sigma**2)))
        else:
            if min_rmse is not None:
                index = np.argmin(self.rmse)
                self.sigma=np.sqrt((1./(len(self.mmi_obs)-6))*self.sum_squares_list[index])
#                self.sigma=(1./(len(self.mmi_obs)-6))*np.power(min_rmse,2)
                print 'sigma', self.sigma
                self.uncert_fun = norm(min_rmse,self.sigma)
            else:
                index = np.argmin(self.rmse)
                self.sigma=np.sqrt((1./(len(self.mmi_obs)-6))*self.sum_squares_list[index])
                #self.sigma=(1./(len(self.mmi_obs)-6))*np.power(min(self.rmse),2)
                print 'sigma', self.sigma
                self.uncert_fun = norm(min(self.rmse),self.sigma)
                #self.log_lik = (-1/2)*self.sum_squares_list[index]/(self.sigma**2) - \
                  #  len(self.mmi_obs)*np.log(self.sigma*np.sqrt(2*np.pi))
                #print 'Log likelihood', self.log_lik
                self.max_lik = np.power((1/(self.sigma*np.sqrt(2*np.pi))), len(self.mmi_obs)) * \
                    np.exp((-1/2)*((self.sum_squares_list[index]/self.sigma**2)))
#                self.max_lik = np.power(10, self.log_lik)
                print 'Maxmium likelihood = ', self.max_lik
                #df = len(self.mmi_obs)-6
                #chi_bv = chi2(df).ppf(0.95)/df
                #self.uncert_fun = min(self.rmse)+ chi2(df).ppf(0.95)/df
            print self.uncert_fun.ppf(0.975)
            indices = np.where(self.rmse < self.uncert_fun.ppf(0.975))[0]
        # Dump all data to a file
        try:
            self.parameter_space
        except AttributeError:
            self.rupture_params_2_array() 
        rmses = np.array(self.rmse)
        llhs = []
        for i, rup in enumerate(self.rupture_list):
            likelihood = np.power((1/(self.sigma*np.sqrt(2*np.pi))), len(self.mmi_obs)) * \
                np.exp((-1/2)*((self.sum_squares_list[i]/self.sigma**2)))
            llhs.append(likelihood)
        llhs = np.array(llhs)
        self.parameter_space_llh = np.vstack([self.parameter_space, rmses, llhs])
        self.parameter_space_llh = self.parameter_space_llh.T
        header = 'mag,lon,lat,depth,strike,dip,rmse,likelihood'
        np.savetxt(filename, self.parameter_space_llh, delimiter=',', header = header)
        # all ruptures within 95% uncertainty bounds 
        self.fitted_ruptures = []
        self.fitted_rmse = []
        for index in indices:
            self.fitted_ruptures.append(self.rupture_list[index])
            self.fitted_rmse.append(self.rmse[index])
        # Create lists for storing ranges of each parameter
        self.fitted_mags = []
        self.fitted_lons = []
        self.fitted_lats = []
        self.fitted_depths = []
        self.fitted_strikes = []
        self.fitted_dips = []
        for rup in self.fitted_ruptures:
            self.fitted_mags.append(rup.mag)
            self.fitted_lons.append(rup.hypocenter.longitude)
            self.fitted_lats.append(rup.hypocenter.latitude)
            self.fitted_depths.append(rup.hypocenter.depth)
            self.fitted_strikes.append(rup.surface.get_strike())
            self.fitted_dips.append(rup.surface.get_dip())
        self.min_mag = min(self.fitted_mags)
        self.max_mag = max(self.fitted_mags)
        self.min_lon = min(self.fitted_lons)
        self.max_lon = max(self.fitted_lons)
        self.min_lat = min(self.fitted_lats)
        self.max_lat = max(self.fitted_lats)
        self.min_depth = min(self.fitted_depths)
        self.max_depth = max(self.fitted_depths)
        self.min_strike = min(self.fitted_strikes)
        self.max_strike = max(self.fitted_strikes)
        self.min_dip = min(self.fitted_dips)
        self.max_dip = max(self.fitted_dips)
        print 'strikes',  np.unique(np.array(self.fitted_strikes))
        print 'dips', np.unique(np.array(self.fitted_dips))

    def parameter_pdf(self, fig_comment='', limits_filename=None):
        """Calculate a pdf for parameter values based on the uncertainty model
        """
        try:
            self.parameter_space
        except AttributeError:
            self.rupture_params_2_array()
        xlabel_dict = {'mag': 'Magnitude ($M_w$)', 'longitude': 'Longitude',
                          'latitude': 'Latitude', 'depth': 'Depth (km)',
                          'strike': 'Strike ($^\circ$)', 'dip':'Dip ($^\circ$)'}    
        #self.parameter_pdfs = self.parameter_space*self.uncert_fun.pdf(self.rmse)
        
        # Now sum equal values to build pdfs
        self.parameter_pdf_sums = {}
        self.parameter_pdf_values = {}
        plt.clf()
        fig = plt.figure(figsize=(16,8))#, tight_layout=True)
        gs = plt.GridSpec(2,4)
        for key, value in self.parameter_dict.iteritems():
            unique_vals = np.unique(self.parameter_space[value])
            pdf_sums = []
            bin=False
            # Determine if we need to bin values of strike and dip (i.e.
            # for non-area sources
            if key == 'strike':
                if len(unique_vals) > 24:
                    bin=True
                    bins =  np.arange(0, 360, 15.)
            if key == 'dip' or key == 'depth':
                if len(unique_vals) > 4:
                    bin=True
                    if key == 'dip':
                        bins = np.arange(min(self.parameter_space[value]), max(self.parameter_space[value])+5, 5.)
                    elif key == 'depth':
                        bins = np.arange(min(self.parameter_space[value]), max(self.parameter_space[value])+20, 20.)
            if bin: # Calculate as histogram
                hist, bins = np.histogram(self.parameter_space[value], bins)
                # align to bin centre for plotting
                #bin_width = bins[1] - bins[0]
                unique_vals = []
                for i, edge in enumerate(bins):
                    try:
                        ind = np.intersect1d(np.where(self.parameter_space[value] >= edge), 
                                             np.where(self.parameter_space[value] < bins[i+1]))
                    except IndexError:
                        ind = np.where(self.parameter_space[value] >= edge)
                    pdf_sum = 0
                    if len(ind) < 1:
                        print ind
                    else:
                        for index in ind:
                            #                        print index
                            try:
                                likelihood = np.power((1/(self.sigma*np.sqrt(2*np.pi))), len(self.mmi_obs)) * \
                                    np.exp((-1/2)*((self.sum_squares_list[index]/self.sigma**2)))
                                pdf_sum += likelihood
                            except TypeError:
                                print 'ind', ind
                                print 'index',index
                    pdf_sums.append(pdf_sum)                              
                    unique_vals.append(edge)# + bin_width)
            else: # Use raw values
                for val in unique_vals:
                    # Just unique values, as pdf sums are repeated
                    ind = np.argwhere(self.parameter_space[value]==val)
                    pdf_sum = 0
#                    print ind
                    for index in ind:
 #                       print index
                        likelihood = np.power((1/(self.sigma*np.sqrt(2*np.pi))), len(self.mmi_obs)) * \
                            np.exp((-1/2)*((self.sum_squares_list[index]/self.sigma**2)))
                        pdf_sum += likelihood
                    #pdf_sum = np.sum(self.uncert_fun.pdf(self.rmse[ind]))
                    pdf_sums.append(pdf_sum)
            # Normalise pdf sums
            pdf_sums = np.array(pdf_sums)
            pdf_sums = pdf_sums/np.sum(pdf_sums)
            print 'pdf_sums', pdf_sums
            self.parameter_pdf_sums[key] = pdf_sums
            self.parameter_pdf_values[key] = unique_vals
            
            # Get the best-fit value for plotting on top
            index = np.argmin(self.rmse)
            best_fit_x =  self.parameter_space[value][index]
            #y_index = np.where(unique_vals == best_fit_x)[0]
            try:
                y_index = (np.abs(unique_vals - best_fit_x)).argmin()[0]
            except IndexError:
                y_index = (np.abs(unique_vals - best_fit_x)).argmin()
            best_fit_y = pdf_sums[y_index]
            # Now plot the results
            try:
                width = unique_vals[1] - unique_vals[0] # Scale width by discretisation
            except IndexError: # if only one value
                width = 1.0
            if key=='strike': 
                # Plot as rose diagram
                ax = fig.add_subplot(gs[0,2],projection='polar')
                ax.bar(np.deg2rad(unique_vals), pdf_sums, width=np.deg2rad(width), bottom=0.0,
                       align='center', color='0.5', edgecolor='k')
                ax.scatter(np.deg2rad(best_fit_x), best_fit_y, marker = '*', color='0.5', edgecolor='k', s=200,zorder=10 )
                ax.set_theta_zero_location('N')
                ax.set_theta_direction(-1)
                ax.set_thetagrids(np.arange(0, 360, 15))
                # Define grids intervals for radial axis
                if max(pdf_sums) < 0.11:
                    r_int = 0.02
                else:
                    r_int = 0.1
                ax.set_rgrids(np.arange(0.01, max(pdf_sums)+0.01, r_int), angle= np.deg2rad(7.5), weight= 'black')
                ax.set_xlabel(xlabel_dict[key])
                ax.text(-0.07, 1.02, 'b)', transform=ax.transAxes, fontsize=14)
                #ax.set_title('Rose Diagram of Strike PDF"', y=1.10, fontsize=15)
                #plt.savefig('%s_%s_%s_parameter_pdf_rose_diagram.png' % (fig_comment,self.gsim, key), \
                 #       dpi=300, format='png', bbox_inches='tight')
            elif key == 'mag':
                ax =  fig.add_subplot(gs[:,0:2])
                ymax = max(pdf_sums)*1.1
                lbx = [self.min_mag, self.min_mag]
                lby = [0, ymax]
                ubx = [self.max_mag, self.max_mag]
                uby = [0, ymax]
            elif key == 'depth':
                ax =  fig.add_subplot(gs[1,3])
                ymax = max(pdf_sums)*1.1
                lbx = [self.min_depth, self.min_depth]
                lby = [0, ymax]
                ubx = [self.max_depth, self.max_depth]
                uby = [0, ymax]
            elif key == 'dip':
                ax =  fig.add_subplot(gs[1,2])
                ymax = max(pdf_sums)*1.1
                lbx = [self.min_dip, self.min_dip]
                lby = [0, ymax]
                ubx = [self.max_dip, self.max_dip]
                uby = [0, ymax]
            if key == 'mag' or key == 'dip' or key == 'depth' :
                ax.bar(unique_vals, pdf_sums, width, align='center', color='0.5')
                ax.scatter(best_fit_x, best_fit_y, marker = '*', color='0.5', 
                           edgecolor='k', s=200, zorder=10)
                if key != 'latitude' and key != 'longitude':
                    ax.plot(lbx, lby, color='k')
                    ax.plot(ubx, uby, color='k')
                ax.set_ylim(0, ymax)
                ax.set_ylabel('Probability')
                ax.set_xlabel(xlabel_dict[key])
                if key == 'mag':
                    ax.text(0.05, 0.95, 'a)', transform=ax.transAxes, fontsize=14)
                if key == 'dip':
                    ax.text(0.05, 0.92, 'd)', transform=ax.transAxes, fontsize=14)
                if key == 'depth':
                    ax.text(0.05, 0.92, 'e)', transform=ax.transAxes, fontsize=14)
        # Now plot a map of location uncertainty
        # First get 2D pdf
        pdf_sums = []
        all_lons = []
        all_lats = []
        for i, lon in enumerate(self.parameter_space[self.parameter_dict['longitude']]):
            if lon in all_lons:
                continue
            else:
                lat = self.parameter_space[self.parameter_dict['latitude']][i]
#               lat = self.parameter_pdf_values['latitude'][i]
                index = np.intersect1d(np.argwhere(self.parameter_space[self.parameter_dict['longitude']]==lon), \
                                           np.argwhere(self.parameter_space[self.parameter_dict['latitude']]==lat))
                pdf_sum = np.sum(self.uncert_fun.pdf(self.rmse[index]))
                pdf_sums.append(pdf_sum)
                all_lons.append(lon)
                all_lats.append(lat)
        # Normalise pdf sums
        pdf_sums = np.array(pdf_sums)
        pdf_sums = pdf_sums/np.sum(pdf_sums)
        self.parameter_pdf_sums['lon_lat'] = pdf_sums
        all_lons = np.array(all_lons)
        all_lats = np.array(all_lats)

        # Get best fit value
        index = np.argmin(self.rmse)
        best_fit_lon =  self.parameter_space[self.parameter_dict['longitude']][index]
        best_fit_lat =  self.parameter_space[self.parameter_dict['latitude']][index]

        ax =  fig.add_subplot(gs[0,3])
        minlon = min(self.parameter_pdf_values['longitude'])
        maxlon = max(self.parameter_pdf_values['longitude'])
        minlat = min(self.parameter_pdf_values['latitude'])
        maxlat = max(self.parameter_pdf_values['latitude'])
        lat_0 = minlat + (maxlat-minlat)/2.
        lon_0 = minlon + (maxlon-minlon)/2.       
        m = Basemap(projection='tmerc', 
                    lat_0=lat_0, lon_0=lon_0,
                    llcrnrlon=minlon,
                    llcrnrlat=minlat,
                    urcrnrlon=maxlon,
                    urcrnrlat=maxlat, 
                    resolution='i')
        m.drawcoastlines(linewidth=0.5,color='k')
        m.drawcountries(color='0.2')
        m.drawstates(color='0.2')
        if maxlon-minlon < 2:
            gridspace = 0.5
        elif maxlon-minlon < 3:
            gridspace = 1.0
        elif maxlon-minlon < 7:
            gridspace = 2.0
        else:
            gridspace = 5.0
        m.drawparallels(np.arange(-90.,90.,gridspace), labels=[1,0,0,0],
                        fontsize=10, dashes=[2, 2], color='0.5',
                        linewidth=0.5)
        m.drawmeridians(np.arange(0.,360.,gridspace), labels=[0,0,0,1],
                        fontsize=10, dashes=[2, 2], color='0.5',
                        linewidth=0.5)
        max_val = max(pdf_sums)*1.1
        print 'pdf_sums', pdf_sums
        clevs = np.arange(0.0,max_val,(max_val/50))
        cmap = plt.get_cmap('gray_r')
        # Adjust resolution to avoid memory intense interpolations
        res = max((maxlon-minlon)/50., (maxlat-minlat)/50.)
        xy = np.mgrid[minlon:maxlon:res,minlat:maxlat:res]
        xx,yy=np.meshgrid(xy[0,:,0], xy[1][0])
        griddata = interpolate.griddata((all_lons, all_lats), pdf_sums, (xx,yy), method='nearest')
        # now plot filled contours of pdf 
        cs = m.contourf(xx, yy, griddata, clevs, cmap=cmap, vmax=max(pdf_sums), vmin=min(pdf_sums), latlon=True)
        # Mask areas outside of source model
        if limits_filename is not None:
            limits_data = np.genfromtxt(limits_filename, delimiter=',')
            limits_x = limits_data[:,0]
            limits_y = limits_data[:,1]
            limits_x, limits_y = m(limits_x, limits_y) # Convert to map coordinates
            poly = Polygon(np.c_[limits_x, limits_y], closed=True)
            clippath =  poly.get_path()
            ax = plt.gca()
            patch = PathPatch(clippath, transform=ax.transData, facecolor='none', linewidth=0.4)
            print 'Adding patch'
            ax.add_patch(patch)
            for contour in cs.collections:
                contour.set_clip_path(patch)
        # Now add best-fit location on top
        m.scatter(best_fit_lon, best_fit_lat, marker = '*', color='w', 
                   edgecolor='k', s=200, zorder=10, latlon=True)
        #m.text(0.05, 0.95, 'c)', transform=ax.transAxes, fontsize=14)
        plt.annotate('c)', xy=(0.05, 0.9),xycoords='axes fraction', fontsize=14)
        if max_val < 0.001:
            loc_int = 0.0005
        elif max_val < 0.01:
            loc_int = 0.005
        else:
            loc_int = 0.05
        ticks = np.arange(0.0, max_val*1.1, loc_int) 
        cbar = m.colorbar(cs, ticks=ticks)
        cbar.ax.set_ylabel('Probability')
        figname = '%s_%s_all_parameter_pdf.png' % (fig_comment,self.gsim)
        figname = figname.replace('()', '')
        plt.tight_layout()
        plt.savefig(figname, dpi=300, format='png', bbox_inches='tight')

    def uncertainty_slice1D(self, z, x, y, xvalue, yvalue):
        """ Get 1D slices of uncertainty model, e.g. to get range of 
        magnitudes at best-fit location
        :params z:
            Quantity we want uncertainty range on
        :params x:
            Quantity to be used for x-axis
        :params y:
            Quantity to be used for y-axis
        :params xvalue:
            Value of x at which the other dimension will be sliced
        :params yvalue:
            Value of y at which the other dimension will be sliced
        """
        try:
            self.parameter_space
        except AttributeError:
            self.rupture_params_2_array()
        
        xindices = np.where(self.parameter_space[self.parameter_dict[x]] == xvalue)[0]
        yindices = np.where(self.parameter_space[self.parameter_dict[y]] == yvalue)[0]
        i = np.intersect1d(xindices,yindices)
        z_subset = self.parameter_space[self.parameter_dict[z]][i]
        rmse_subset = self.rmse[i]
        indices = np.where(rmse_subset < self.uncert_fun.ppf(0.975))[0]
        zs = z_subset[indices]
        min_z = min(zs)
        max_z = max(zs)
        return min_z, max_z

    def uncertainty_slice2D(self, x, y, z, zvalue, fig_comment=None, limits_filename=None):
        """ Get 2D slices of uncertainty model for plotting
        :params x:
            Quantity to be used for x-axis
        :params y:
            Quantity to be used for y-axis
        :params z:
            Quantity we are slicing across
        :params zvalue:
            Value of z at which the othe dimensions will be sliced
        :params fig_comment:
            String to be appended to figure name
        :params limits_file:
            File containing boundary of parameter space that can be used to mask data, 
            e.g. to not show interpolated data outside the boundaries of the geographic
            extent of the source model. This file should be csv formatted containing
            longitude,latitude. The file should not contain a header.
        """

        # Create an array of rupture parameters
        try:
            self.parameter_space
        except AttributeError:
            self.rupture_params_2_array()

        indices = np.where(self.parameter_space[self.parameter_dict[z]] == zvalue)
        xvalues = self.parameter_space[self.parameter_dict[x], indices][0]
        yvalues = self.parameter_space[self.parameter_dict[y], indices][0]
        rmse_subset = self.rmse[indices]

        # now get lowest rmse for each combination of remaining parameters at each 
        # xy point
        rmse_list = []
        xs = []
        ys = []
        for xi in np.unique(xvalues):
            for yi in np.unique(yvalues):
                j = np.where(xvalues==xi)
                k = np.where(yvalues == yi)
                i = np.intersect1d(j,k) # Get all locations matching both x, y locations
                if len(i) > 0:
                    rmse_list.append(min(rmse_subset[i]))
                    xs.append(xi)
                    ys.append(yi)
        rmse_subset = np.array(rmse_list)
        xs = np.array(xs)
        ys = np.array(ys)
        xx,yy = np.mgrid[min(xvalues):max(xvalues):0.02,min(yvalues):max(yvalues):0.01]
        rmse_grid = interpolate.griddata((xs, ys), rmse_subset, (xx,yy), method='nearest')

        plt.clf()
        # Get data for clipping
        if limits_filename is not None:
            limits_data = np.genfromtxt(limits_filename, delimiter=',')
            limits_x = limits_data[:,0]
            limits_y = limits_data[:,1]
            clippath = Path(np.c_[limits_x, limits_y])
            patch = PathPatch(clippath, facecolor='none')
        try:
            CS1=plt.contour(xx, yy, rmse_grid, levels = [(min(self.rmse) + self.sigma), self.uncert_fun.ppf(0.975)], linewidths=0.5, colors='k')
        except AttributeError:
            pass
        except ValueError:
            print 'only one best fit locations, cannnot plot locations uncertainty'
            pass
#        CS2=plt.contourf(xx, yy, rmse_grid, 8,vmax=np.max(rmse_grid), vmin=np.min(rmse_grid))
        try:
            CS2=plt.contourf(xx, yy, rmse_grid, 8,vmax=max(self.rmse), vmin=min(self.rmse))
            cbar = plt.colorbar(CS2)
            cbar.ax.set_ylabel('RMSE')
        except ValueError:
            print 'only one best fit locations, cannnot plot locations uncertainty'
            pass
#        CS2=plt.contourf(xs, ys, rmse_subset, 8,vmax=max(self.rmse), vmin=min(self.rmse)) 
            # add patch to mask no data area
        try:
            ax = plt.gca()
            ax.add_patch(patch)
        except:
            print 'No limits file'
            pass   

        try:
            figname = '%s_rmse_slice_%s_%.2f_%s_%s_%s.png' % (
                fig_comment, z, zvalue, x, y, self.gsim)
            figname = figname.replace('()', '')
            plt.savefig(figname)
        except ValueError:
            print 'only one best fit locations, cannnot plot locations uncertainty'
            pass
        
                                        
        return xs, ys, rmse_subset

    

