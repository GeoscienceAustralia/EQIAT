"""Contains code to calculate Bayesian posterior distributions for 
parameters calculated using estimate_magnitude.py

Jonathan Griffin
Geoscience Australia 
July 2018
"""

import numpy as np
from mpl_toolkits.basemap import Basemap, maskoceans
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Polygon
from scipy import interpolate

data_file = 'outputs/1847_BooreEtAl2014_parameter_llh.csv'

def update_weights(parameter_space, prior_pdfs):
    """Update weights in a Bayesian sense
    """
    prior_weights = []
    llhs = parameter_space[7]
    print llhs, sum(llhs)
    print max(llhs), min(llhs)
    parameter_space = parameter_space.T
    for combo in parameter_space:
#        print combo
        i0 = np.where(prior_pdfs[0][0]==combo[0])
        i1 = np.where(prior_pdfs[0][1]==combo[1])
        i2 = np.where(prior_pdfs[0][2]==combo[2])
        i3 = np.where(prior_pdfs[0][3]==combo[3])
        i4 = np.where(prior_pdfs[0][4]==combo[4])
        i5 = np.where(prior_pdfs[0][5]==combo[5])
        prior_weight = prior_pdfs[1][0][i0] * prior_pdfs[1][1][i1] * \
            prior_pdfs[1][2][i2] * prior_pdfs[1][3][i3] * \
            prior_pdfs[1][4][i4] * prior_pdfs[1][5][i5]

        prior_weights.append(prior_weight)
#    print prior_weights
    prior_weights = np.array(prior_weights).flatten()
    print 'priors', prior_weights, sum(prior_weights)
    print max(prior_weights), min(prior_weights)
    posterior_probs = llhs*prior_weights/sum(llhs*prior_weights)#prior_weights)#denominator
    print 'updates', posterior_probs, max(posterior_probs), min(posterior_probs)
    print 'sum', sum(posterior_probs)
    return posterior_probs

def parameter_pdf(parameter_space, fig_comment='', limits_filename=None):
    """Calculate a pdf for parameter values based on the uncertainty model                                
    """

    xlabel_dict = {'mag': 'Magnitude ($M_w$)', 'longitude': 'Longitude',
                   'latitude': 'Latitude', 'depth': 'Depth (km)',
                   'strike': 'Strike ($^\circ$)', 'dip':'Dip ($^\circ$)'}
    parameter_dict = {'mag': 0, 'longitude': 1,
                      'latitude': 2, 'depth': 3,
                      'strike': 4, 'dip':5}
    parameter_pdf_sums = {}
    parameter_pdf_values = {}
    plt.clf()
    fig = plt.figure(figsize=(16,8))#, tight_layout=True)                                                 
    gs = plt.GridSpec(2,4)
    for key, value in parameter_dict.iteritems():
        unique_vals = np.unique(parameter_space[value])
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
                    bins = np.arange(min(parameter_space[value]), max(parameter_space[value])+5, 5.)
                elif key == 'depth':
                        bins = np.arange(min(parameter_space[value]), max(parameter_space[value])+20, 20.)
        if bin: # Calculate as histogram                                                                                                           
            hist, bins = np.histogram(parameter_space[value], bins)
            # align to bin centre for plotting                                                                                                     
                #bin_width = bins[1] - bins[0]                                                                                                         
            unique_vals = []
            for i, edge in enumerate(bins):
                try:
                    ind = np.intersect1d(np.where(parameter_space[value] >= edge),
                                         np.where(parameter_space[value] < bins[i+1]))
                except IndexError:
                    ind = np.where(parameter_space[value] >= edge)
                pdf_sum = 0
                for index in ind:
#                    likelihood = np.power((1/(self.sigma*np.sqrt(2*np.pi))), len(self.mmi_obs)) * \
#                        np.exp((-1/2)*((self.sum_squares_list[index]/self.sigma**2)))
                    posterior_prob = parameter_space[7][index]
                    pdf_sum += posterior_prob
                    #pdf_sum = np.sum(self.uncert_fun.pdf(self.rmse[ind]))                                                                             
                pdf_sums.append(pdf_sum)
                unique_vals.append(edge)# + bin_width)      
        else: # Use raw values                                                                                                                     
            for val in unique_vals:
                # Just unique values, as pdf sums are repeated                                                                                     
                ind = np.argwhere(parameter_space[value]==val)
                pdf_sum = 0
                for index in ind:
#                    likelihood = np.power((1/(self.sigma*np.sqrt(2*np.pi))), len(self.mmi_obs)) * \
#                        np.exp((-1/2)*((self.sum_squares_list[index]/self.sigma**2)))
                    posterior_prob = parameter_space[7][index]
                    pdf_sum += posterior_prob
                    #pdf_sum = np.sum(self.uncert_fun.pdf(self.rmse[ind]))                                                                             
                pdf_sums.append(pdf_sum)
        # Normalise pdf sums       
        pdf_sums = np.array(pdf_sums)
        #pdf_sums = pdf_sums/np.sum(pdf_sums)
        print 'pdf_sums', pdf_sums
        parameter_pdf_sums[key] = pdf_sums
        parameter_pdf_values[key] = unique_vals

        # Get the best-fit value for plotting on top                                                                                               
        index = np.argmin(parameter_space[6])
        best_fit_x =  parameter_space[value][index]
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
        elif key == 'mag':
            ax =  fig.add_subplot(gs[:,0:2])
            ymax = max(pdf_sums)*1.1
#            lbx = [self.min_mag, self.min_mag]
#            lby = [0, ymax]
#            ubx = [self.max_mag, self.max_mag]
#            uby = [0, ymax]
        elif key == 'depth':
            ax =  fig.add_subplot(gs[1,3])
            ymax = max(pdf_sums)*1.1
#            lbx = [self.min_depth, self.min_depth]
#            lby = [0, ymax]
#            ubx = [self.max_depth, self.max_depth]
#            uby = [0, ymax]
        elif key == 'dip':
            ax =  fig.add_subplot(gs[1,2])
            ymax = max(pdf_sums)*1.1
#            lbx = [self.min_dip, self.min_dip]
#            lby = [0, ymax]
#            ubx = [self.max_dip, self.max_dip]
#            uby = [0, ymax]
        if key == 'mag' or key == 'dip' or key == 'depth' :
            ax.bar(unique_vals, pdf_sums, width, align='center', color='0.5')
            ax.scatter(best_fit_x, best_fit_y, marker = '*', color='0.5',
                       edgecolor='k', s=200, zorder=10)
            #if key != 'latitude' and key != 'longitude':
#            ax.plot(lbx, lby, color='k')
#            ax.plot(ubx, uby, color='k')
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
    for i, lon in enumerate(parameter_space[parameter_dict['longitude']]):
        if lon in all_lons:
            continue
        else:
            lat = parameter_space[parameter_dict['latitude']][i]
#               lat = self.parameter_pdf_values['latitude'][i]                                                                                         
            index = np.intersect1d(np.argwhere(parameter_space[parameter_dict['longitude']]==lon), \
                                       np.argwhere(parameter_space[parameter_dict['latitude']]==lat))
            pdf_sum = np.sum(parameter_space[7][index])#uncert_fun.pdf(rmse[index]))
            pdf_sums.append(pdf_sum)
            all_lons.append(lon)
            all_lats.append(lat)
    # Normalise pdf sums                                                                                                                           
    pdf_sums = np.array(pdf_sums)
#    pdf_sums = pdf_sums/np.sum(pdf_sums)
    parameter_pdf_sums['lon_lat'] = pdf_sums
    all_lons = np.array(all_lons)
    all_lats = np.array(all_lats)
    # Get best fit value                                                                                                                           
    index = np.argmin(parameter_space[6])
    best_fit_lon =  parameter_space[parameter_dict['longitude']][index]
    best_fit_lat =  parameter_space[parameter_dict['latitude']][index]
    ax =  fig.add_subplot(gs[0,3])
    minlon = min(parameter_pdf_values['longitude'])
    maxlon = max(parameter_pdf_values['longitude'])
    minlat = min(parameter_pdf_values['latitude'])
    maxlat = max(parameter_pdf_values['latitude'])
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
    figname = '%s_all_parameter_pdf.png' % (fig_comment)
    figname = figname.replace('()', '')
    plt.tight_layout()
    plt.savefig(figname, dpi=300, format='png', bbox_inches='tight')

if __name__ == "__main__":
    parameter_space = np.genfromtxt(data_file, delimiter=',', skip_header=1)
    fig_comment = 'figures/' + data_file.split('/')[1][:-4]
    print fig_comment
    parameter_space = parameter_space.T
    # Set up prior pdfs - set-up using basic assumptions and limits of data
    # magnitude - based on Gutenberg-Richter assuming b value = 1, and that
    # CDF from mmin to mmax = 1
    mags = np.unique(parameter_space[0])
    mmax = max(mags)
    mmin = min(mags)
    b=1.
    a = np.log10(1./(np.power(10,-1*b*mmin) - np.power(10, -1*b*mmax)))
    print a
    # Now we need to generate an incremental pdf 
    reversed_mag_priors = []
    reversed_mags = list(reversed(mags))
    for i, mag in enumerate(reversed_mags):
        if i == 0:
            prior = np.power(10, a - b*mag)
        else:
            prior = np.power(10, a - b*mag) - np.power(10, a - b*reversed_mags[i-1])
        reversed_mag_priors.append(prior)
    mag_priors = np.array(list(reversed(reversed_mag_priors)))
    print mags
    print mag_priors, sum(mag_priors)

    # longitude, latitude, strike, depth and dip - uniform across parameter space
    lon_priors = np.ones(len(np.unique(parameter_space[1]))) * \
        (1./len(np.unique(parameter_space[1])))
    lat_priors = np.ones(len(np.unique(parameter_space[2]))) * \
        (1./len(np.unique(parameter_space[2])))
    depth_priors = np.ones(len(np.unique(parameter_space[3]))) * \
        (1./len(np.unique(parameter_space[3])))
    strike_priors = np.ones(len(np.unique(parameter_space[4]))) * \
        (1./len(np.unique(parameter_space[4])))
    dip_priors = np.ones(len(np.unique(parameter_space[5]))) * \
        (1./len(np.unique(parameter_space[5])))

    # Calculate normalising denominator
#    denom = 0
#    for m in mag_priors:
#        for lon in lon_priors:
#            for lat in lat_priors:
#                for dep in depth_priors:
#                    for strike in strike_priors:
#                        for dip in dip_priors:
#                            denom += m*lon*lat*dep*strike*dip
#    print denom
    priors = np.array([[np.unique(parameter_space[0]), np.unique(parameter_space[1]),
                       np.unique(parameter_space[2]), np.unique(parameter_space[3]),
                       np.unique(parameter_space[4]), np.unique(parameter_space[5])],
                      [mag_priors, lon_priors, lat_priors,
                       depth_priors, strike_priors, dip_priors]])
    posterior_probs = update_weights(parameter_space, priors)

    print parameter_space
    parameter_space[7] = posterior_probs
    parameter_pdf(parameter_space, fig_comment = fig_comment)
