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
    """For each gronud motion realisation, convert to 
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
    from matplotlib import pyplot
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
    from matplotlib import pyplot
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

if __name__ == "__main__":
    filename = '../hazard_curves_1.0s.csv'
    period = 1.0
    site_class = 'D'
    time = 69 # number of years in time interval - annual rates
    time_jkt = 196 # complmeteness for Jakarta
    #    f_out = open('results_poisson_simulation_test.csv', 'w')
    data = numpy.genfromtxt(filename, delimiter= ',', skip_header=1)
    f_in = open(filename, 'rU')
    header = f_in.readline()
    f_in.close()

    gm = data[:,0]
#    print gm
    # Calculate incremental rate
    # Interpolation
    #for i in range(1,6):
    #    print i
    # Just test for Jakarta for now
    hazard_curve = numpy.array([gm, data[:,1]])
    # interpolate hazard curve
    interpolate = False
    if interpolate:
        x_values = numpy.power(10, numpy.arange(-4, 1.7, 0.1))
        hazard_curve_interp = numpy.interp(x_values, gm, data[:,1])
    #    print 'hazard_curve_interp', hazard_curve_interp
        hazard_curve = numpy.array([x_values, hazard_curve_interp])
    plot_hazard_curve(hazard_curve)
   # print hazard_curve
    incremental_hazard_curve = calculate_incremental_rates(hazard_curve)
   # print incremental_hazard_curve
    realisations = \
        simulate_gm_realisations(incremental_hazard_curve,
                                 time_jkt, 1000)
#    print realisations
    #Convert to MMI including uncertainty
    #Remove site effects first
    # i.e. for each realisation we want to consider 
    # how we might observe it
    gm_site_effects = gm_site_amplification(hazard_curve[0], period, site_class)
    mmi_sample_database, mmi_sample_cumulative_database = build_mmi_samples(gm_site_effects, realisations, period)
 #   print mmi_sample_cumulative_database
    percentile_97_5 = numpy.percentile(mmi_sample_cumulative_database, 97.5, axis=0)
    percentile_2_5 = numpy.percentile(mmi_sample_cumulative_database, 2.5, axis=0)
    median = numpy.median(mmi_sample_cumulative_database, axis=0)
    print percentile_97_5
    print percentile_2_5
    print median
    plot_mmi_samples(mmi_sample_database, mmi_sample_cumulative_database)
