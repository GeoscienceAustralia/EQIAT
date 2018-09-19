"""Code to simulate random realisations from a Poisson
process with lamba taken from the mean hazard curve.
This is designed to test the Indonesian earthquake
hazard map with historical intensity observations

Jonathan Griffin
Geoscience Australia
July 2016
"""

import sys
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

def poisson_samples(lam, time, size=1):
    """Generate random poisson samples
    :param lam:
        Int/float: mean (lambda) of the Poisson 
        distribution, i.e. the mean rate taken from
        the earthquake hazard curve for a given 
        intensity level
    :param  time:
        Int/float: time period being sampled in years
    :param size:
        Int/tuple: Number of realisations/output shape.
        If the given shape is, e.g., (m, n, k), 
        then m * n * k samples are drawn.
    :returns:
        Samples: ndarry or scalar of samples depending on
        size parameter."""

    lam_time = lam*time # Rate over time period
    samples = numpy.random.poisson(lam_time, size)

    return samples

def simulate_gm_realisations(incremental_hazard_curve,
                             time, num_realisations):
    """Simulate synthetic ground motion realisations
    from discretised incremental hazard curve
    :param incremental_hazard_curve:
        N*2 array with first column containing intensity measure,
    second containing incremental annual rates
    :param time:
        Time period (years) over which to simulate the
        ground motion realisations
    :param num_realisations:
        Integer for number of realisations to sample
    :returns:
        realisations: ndarry of sample ground motion realisations
    """
    sample_list = []
    for i in range(len(incremental_hazard_curve[0])):
        samples = poisson_samples(incremental_hazard_curve[1][i],
                                  time, num_realisations)
        sample_list.append(samples)
    sample_array = numpy.array(sample_list)
    sample_array = sample_array.T
    #print 'sample_array', sample_array
    return sample_array

def gm_site_amplification(gm, period, site_class):
    """Amplify ground motions based on NEHRP
    site class (Borchedt 1994).
    :param gm:
        Ground motion values (g) for bedrock 
        (vs30 = 760 m/s) 
    :param period:
        Response spectral acceleration period for
        which ground motions are defined
    :param site_class:
        NEHRP site class
    :returns amplified_gm:
        Amplified ground motion values
    """
    from RSA_VS30_to_MMI import Amp_fns
    function = Amp_fns()
    RSA_amp_list = []
    if period < 1.0:
        for i in range(len(gm)):
            RSA_amp = function.short_period_dict[site_class](gm[i])*gm[i]
            RSA_amp_list.append(RSA_amp)
    else:
        for i in range(len(gm)):
            RSA_amp = function.mid_period_dict[site_class](gm[i])*gm[i]
            RSA_amp_list.append(RSA_amp)
    amplified_gm = numpy.array(RSA_amp_list)
    return amplified_gm

def build_mmi_samples(gm, realisations, period):
    """For each ground motion realisation, convert to 
    MMI allowing while sampling uncertainty in the conversion
    :param gm:
        1D Array of ground motion values assuming already
        including site effects.
    :param realisations:
        N*len(gm) array with each row a 1D array of 
        that contains the number of occurences of the
        corresponding ground motion level.
    :param period:
        Response spectral acceleration period for 
        which ground motions are defined
    :returns mmi_samples:
        Array containing number of times each discrete MMI values 
        occurs for each realisation. MMI values
        are binned such that e.g MMI 5 contains all values
        from 4.5 <= MMI < 5.
    """
    from RSA2MMI import rsa2mmi
    from numpy.random import normal
    # Convert gm to MMI
    gm_mmi = rsa2mmi(gm, period=period, include_uncertainty=False)
   # print 'gm_mmi', gm_mmi
    mmi_sample_list = []
    mmi_sample_cumulative_list = []
    for sample in realisations:
        #print sample
        mmi_list = []
        for i in range(len(sample)):
            mmi_samples = normal(gm_mmi[i], 1.0, sample[i])
            mmi_list.append(mmi_samples)
        mmi_array = numpy.array(mmi_list)
        mmi_array = numpy.hstack(mmi_array)
        #print mmi_array
        mmi_array[mmi_array<0]=0
        mmi_array_round = numpy.round(mmi_array, decimals=0)
        mmi_array_round = mmi_array_round.astype(int)
        mmi_occurrences = numpy.bincount(mmi_array_round, minlength = 14)
        # Now calculate exceedence rates
        mmi_occurrences_reversed = mmi_occurrences[::-1]
        mmi_cumsum = numpy.cumsum(mmi_occurrences_reversed)
        mmi_cumsum = mmi_cumsum[::-1]
        mmi_sample_list.append(mmi_occurrences)
        mmi_sample_cumulative_list.append(mmi_cumsum)
   # print 'sample_array'
    #print mmi_sample_list
#    print mmi_samples
    return mmi_sample_list, mmi_sample_cumulative_list

def plot_mmi_samples(mmi_sample_list, mmi_cumulative_list):
    """Plot mmi exceedence curves for each MMI sample
    :param mmi_sample_list:
        list of tuples containing 2D arrays of
        mmi values and number of occurrences
    """
    pyplot.clf()
    x_values = numpy.arange(0,14)
    print x_values
    for i in range(len(mmi_sample_list)):
        pyplot.plot(x_values, mmi_sample_list[i])
    pyplot.savefig('mmi_incremental_occurrences.png')
    # plot cumulative
    pyplot.clf()
    for i in range(len(mmi_cumulative_list)):
        pyplot.semilogy(x_values, mmi_cumulative_list[i])
    pyplot.savefig('mmi_cumulative_occurrences.png')
                       
def plot_hazard_curve(hazard_curve):
    """Plot the hazard curve
    :param hazard_ucrve:
        2D array of hazrad values, incremental_rate
    """
    pyplot.loglog(hazard_curve[0], hazard_curve[1], marker='o')
    pyplot.savefig('incremental_hazard_rates.png')
    pyplot.clf()
    cumulative_rates = []
    hazard_rates_reversed = hazard_curve[1][::-1]
    cumulative_rates = numpy.cumsum(hazard_rates_reversed)
    cumulative_rates = cumulative_rates[::-1]
    #print hazard_curve[0]
    #print cumulative_rates
    pyplot.loglog(hazard_curve[0], cumulative_rates, marker='o')
    pyplot.savefig('cumulative_rates.png')

def calculate_incremental_rates(hazard_curve):
    """Calculate incremental rates from a cumulative
    hazard curve
    :param hazard_curve:
        N*2 array with first column containing intensity measure,
        second column containing rate
    :return incremental_hazard_curve
        N*2 array with first column containing intensity measure,
        second column containing incremental rate.
    """
    rates = hazard_curve[1]
    #print rates
    inc_rates = numpy.zeros(len(rates)) 
    for i in range(len(rates)):
        if i==len(rates)-1:
            inc_rates[i] = rates[i]
        else:
            inc_rates[i] = rates[i] - rates[i+1]
    incremental_hazard_curve = numpy.array([hazard_curve[0], inc_rates])
    return incremental_hazard_curve

def read_mmi_obs(filename):
    """Read the MMI data into a dictionary
    """
    f_in = open(filename, 'r')
    header = f_in.readline()
    mmi_obs_dict = {}
    for line in f_in.readlines():
        row = line.split(',')
        city = row[0]
        try:
            mmi_obs_dict[city][0].append(int(row[1]))
            mmi_obs_dict[city][1].append(int(row[3]))
        except KeyError:
            mmi_obs_dict[city] = [[int(row[1])], [int(row[3])]]
    return mmi_obs_dict

def plot_mmi_hazmap_and_obs(median, percentile1, percentile2, mmi_obs, city,
                            years = 69, ax=None):
    """Plot the hazard map MMI occurrence percentiles against the 
    number of historical observations
    """
    x_values = numpy.arange(0,14)
    if ax is None:
        pyplot.clf()
    # Change 0 values to a small number for plotting purposes
#    median[median==0] = 0.01
    a, = pyplot.semilogy(x_values, median, color='k')
    b, = pyplot.semilogy(x_values, percentile1, color='b', linestyle='--')
    pyplot.semilogy(x_values, percentile2, color='b', linestyle='--')
    c = pyplot.scatter(mmi_obs[city][0], mmi_obs[city][1])
    # Add completeness
    xx = [6.5, 6.5]
    yy = [0.8, 1000]
    pyplot.plot(xx,yy,linestyle = ':', color='0.3')
    pyplot.xlim(2.9,10)
    pyplot.ylim(0.8, 1000)
    ax.set_title(city)
    ax.set_xlabel('MMI')
    ax.set_ylabel('Number of exceedances in %i years' % years)
    return a,b,c
#    pyplot.savefig('mmi_hazmap_and_observations_%s.png' % city)

if __name__ == "__main__":
    filename = 'data/2017_hazard_curves_PGA.csv'
    filename2 = 'data/hazard_curves_PGA.csv'
    #filename2 = 'data/2017_hazard_curves_1.0s.csv'
    #filename = 'data/hazard_curves_1.0s.csv'
    obs_filename = 'data/City_MMI_all_mentions.csv'
    period = 0.0
    site_class = 'C'
    cities = ['Jakarta', 'Bandung', 'Semarang', 'Yogyakarta', 'Surabaya']
#    cities = ['Semarang', 'Yogyakarta', 'Surabaya']
    time = 69 # number of years in time interval - annual rates
    time_others = 69
    time_jkt = 196 # complmeteness for Jakarta
    mmi_obs_dict = read_mmi_obs(obs_filename)
    print mmi_obs_dict                                                                                      
    data = numpy.genfromtxt(filename, delimiter= ',', skip_header=1)
    data2 = numpy.genfromtxt(filename2, delimiter= ',', skip_header=1)
    data_list = [data, data2]
    gm = data[:,0]
    haz_mmi_dict = {}
#    city = 'Jakarta'
    for city in cities:
        print 'Doing city %s' % city
        if city == 'Jakarta':
            time = time_jkt
            site_class = 'D'
            data_index = 1
        elif city == 'Bandung':
            time = time_others
            site_class = 'C'
            data_index = 2
        elif city == 'Semarang':
            time = time_others
            site_class = 'D'
            data_index = 3
        elif city == 'Yogyakarta':
            time = time_others
            site_class = 'C'
            data_index = 4
        elif city == 'Surabaya':
            time = time_others
            site_class = 'D'
            data_index = 5
        # Calculate incremental rate
        hazard_curve_2017 = numpy.array([gm, data_list[0][:,data_index]])
        hazard_curve_2010 = numpy.array([gm, data_list[1][:,data_index]])
        # interpolate hazard curve
        interpolate = False
        if interpolate:
            x_values = numpy.power(10, numpy.arange(-4, 1.7, 0.1))
            hazard_curve_interp = numpy.interp(x_values, gm, data[:,1])
            hazard_curve = numpy.array([x_values, hazard_curve_interp])
        plot_hazard_curve(hazard_curve_2010)
        plot_hazard_curve(hazard_curve_2017)
        haz_mmi_dict_list = []
        for hazard_curve in [hazard_curve_2010, hazard_curve_2017]:
            incremental_hazard_curve = calculate_incremental_rates(hazard_curve)
            realisations = \
                simulate_gm_realisations(incremental_hazard_curve,
                                         time, 100000)
        #Convert to MMI including uncertainty
        #Remove site effects first
        # i.e. for each realisation we want to consider 
        # how we might observe it
            gm_site_effects = gm_site_amplification(hazard_curve[0], period, site_class)
            mmi_sample_database, mmi_sample_cumulative_database = build_mmi_samples(gm_site_effects, realisations, period)

            percentile_97_5 = numpy.percentile(mmi_sample_cumulative_database, 97.5, axis=0,
                                               interpolation='linear') # nearest/linear
            percentile_2_5 = numpy.percentile(mmi_sample_cumulative_database, 2.5, axis=0,
                                              interpolation='linear')
            median = numpy.median(mmi_sample_cumulative_database, axis=0)
            median[median == 0] = 0.0000001
            percentile_2_5[percentile_2_5 == 0] = 0.0000001
            percentile_97_5[percentile_97_5 == 0] = 0.0000001
            print percentile_97_5
            print percentile_2_5
            print 'median', median
            # Don't always need to do this, it takes a long time
            #    plot_mmi_samples(mmi_sample_database, mmi_sample_cumulative_database)
            haz_mmi_dict[city]=[median, percentile_97_5, percentile_2_5]
            haz_mmi_dict_list.append(haz_mmi_dict)
    # Now plot all onto one figure
    pyplot.clf()
    fig = pyplot.figure(figsize=(15,20))
    pyplot.subplots_adjust(hspace=0.4, wspace = 0.2)
    i = 1
    for city in cities:
        ax = pyplot.subplot(3,2,i)
        if city=='Jakarta':
            time = time_jkt
            # Add figure part label
            if period == '0.0':
                ax.text(-0.1, 1.05, 'a)', fontsize=18, transform=ax.transAxes) 
            else:
                ax.text(-0.1, 1.05, 'b)', fontsize=18, transform=ax.transAxes)
        else: 
            time = time_others
        a,b,c = plot_mmi_hazmap_and_obs(haz_mmi_dict_list[0][city][0], haz_mmi_dict_list[0][city][1],
                                haz_mmi_dict_list[0][city][2], mmi_obs_dict, city, years=time,
                                ax=ax)
        a,b,c = plot_mmi_hazmap_and_obs(haz_mmi_dict_list[1][city][0], haz_mmi_dict_list[1][city][1],
                                haz_mmi_dict_list[1][city][2], mmi_obs_dict, city, years=time,
                                ax=ax)
        i+=1
    fig.legend((a,b,c), ('Median number of exceedances \n from the hazard curve',
                         '2.5 and 97.5 percentiles of \n number of exceedances \n from the hazard curve',
                         'Number of observed exceedances'), 'lower right',
               bbox_to_anchor=(0.85, 0.19))
    pyplot.savefig('mmi_2017_hazmap_and_observations_%.1f.png' % period)
        

