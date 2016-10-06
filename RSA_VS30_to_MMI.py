"""Converts 1 hz RSA from a bedrock hazard map to soil values
using the NEHRP amplifcation factors
Reference: Borchedt 1994

RSA is then converted to MMI using the formuala of Atkinson and Kaka (2006)

Jonathan Griffin, AIFDR January 2011
"""

import sys, os
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import norm
from RSA2MMI import rsa2mmi9

class Amp_fns(object):
    
    def __init__(self):
        """Define amp factors and interpolate into continuous functions
        """
        # Define amp factors and periods
        ground_motions = [0.0, 0.1, 0.2, 0.3, 0.4]
        classes = ['A','B','C','D','E']
        short_period_amp_factors = [[0.9, 0.9 , 0.9, 1., 1.],
                                    [1., 1., 1., 1., 1.],
                                    [1.3, 1.3, 1.2, 1.1, 1.],
                                    [1.6, 1.6, 1.4, 1.1, 0.9],
                                    [2.0, 2.0, 1.6, 1.2, 0.9]]
        mid_period_amp_factors = [[0.8, 0.8, 0.8, 0.8, 0.8],
                                  [1., 1., 1., 1., 1.],
                                  [1.5, 1.5, 1.5, 1.4, 1.4],
                                  [2.3, 2.3, 2.2, 2.0, 1.8],
                                  [3.5, 3.5, 3.2, 2.8, 2.4]]

        # Interpolate to continuous functions using linear interpolation
        self.short_period_dict = {}
        for i in range(len(short_period_amp_factors)):
            short_period_fn = interp1d(ground_motions, short_period_amp_factors[i],
                                       bounds_error = False,
                                       fill_value = short_period_amp_factors[i][4])
            self.short_period_dict[classes[i]] = short_period_fn
        self.mid_period_dict = {}    
        for i in range(len(mid_period_amp_factors)):
            mid_period_fn = interp1d(ground_motions, mid_period_amp_factors[i],
                                     bounds_error = False,
                                     fill_value = mid_period_amp_factors[i][4])
            self.mid_period_dict[classes[i]] = mid_period_fn


def hazmap2amp(RSA1, NEHRP_class):
    """Function to amplify RSA amplifcation based on BS30 values and the NEHRP
        amplification factors
    """
    function = Amp_fns()
    RSA_amp_list = []
    for i in range(len(RSA1)):
        RSA_amp = function.mid_period_dict[NEHRP_class[i]](RSA1[i])*RSA1[i]
        RSA_amp_list.append(RSA_amp)
    return RSA_amp_list


def read_data(infile):
    """Read data into numpy array
    """
    f_in=open(infile, 'r')
    header = f_in.readline()

    RSA1 = []
    vs30 = []
    for line in f_in.readlines():
        row = line.split(',')
        RSA1.append(float(row[2]))
        vs30.append(float(row[3]))
        
    f_in.close()

    return RSA1, vs30

def vs30_to_NEHRP_class(vs30):
    
    vel1 = 760
    vel2 = 360
    vel3 = 180

    class1 = "B"
    class2 = "C"
    class3 = "D"
    class4 = "E"

    NEHRP_class_list = []

    for vs in vs30:       

        # Find site class using vs30
        if vs >= vel1:
            siteclass = class1
        elif vs < vel1 and vs >= vel2:
            siteclass = class2
        elif vs < vel2 and vs >= vel3:
            siteclass = class3
        elif vs < vel3:
            siteclass = class4

        NEHRP_class_list.append(siteclass)

    return NEHRP_class_list    

def write_data(infile, outfile, NEHRP_class, RSA_amp_list, MMI):
    """Write output data to file
    """
    f_in = open(infile, 'r')
    header = f_in.readline()
    f_out = open(outfile, 'w')
    i = 0

    # Write header
    f_out.write('LONGITUDE,LATITUDE,BEDROCK_RSA1,VS30,SITE_CLASS,SOIL_RSA,MMI\n')
    for line in f_in.readlines():
        row = line.rstrip('\n').rstrip('\r')
        row = row + ',' + NEHRP_class[i] + ',' + str(RSA_amp_list[i]) +  ','\
        + str(MMI[i]) +'\n'
        f_out.write(row)        
        i+=1
    f_in.close()
    f_out.close()

if __name__ == '__main__':
    haz_NEHRP_class_file = sys.argv[1]
    RSA1, vs30 = read_data(haz_NEHRP_class_file)
    NEHRP_class = vs30_to_NEHRP_class(vs30)
    RSA_amp_list = hazmap2amp(RSA1, NEHRP_class)
    MMI = rsa2mmi9(RSA_amp_list, period = 1.0)
    outfile = haz_NEHRP_class_file[:-4] + '_MMI.csv'
    write_data(haz_NEHRP_class_file, outfile, NEHRP_class, RSA_amp_list, MMI)
    
    
    
