"""Test building site collection at interpolated points and
writing to nrml format

Can run in parallel in order to deal with large site collections
Jonathan Griffin
Geoscience Australia
"""

import os, sys
import numpy as np
from time import localtime, strftime, gmtime
import string
import pypar
from get_site_model import get_site_collection
from openquake.hazardlib.geo.point import Point
#from estimate_magnitude import build_site_col

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def build_site_col(sites_data, site_model_file, filename=None):
    """Interpolate vs30 values to sites of interest
    """
#    sites_data = np.genfromtxt(site_file)
    site_points = []
    for i in range(len(sites_data[:,0])):
        site_pt = Point(sites_data[i,0], sites_data[i,1])
        site_points.append(site_pt)
    sitecol = get_site_collection(site_model_file, site_points, None, filename)
    return sitecol


# Set up paralell
proc = pypar.size()                # Number of processors as specified by mpirun                     
myid = pypar.rank()                # Id of of this process (myid in [0, proc-1])                     
node = pypar.get_processor_name()  # Host name on which current process is running                   
print 'I am proc %d of %d on node %s' % (myid, proc, node)
#nruns = 320 # currently hard coded - need to improve this                                            
t0 = pypar.time()


locations_file = 'data/jawa_bali_nt_sulawesi_sites_clean.csv'#jawa_sites.csv'
sites_file = 'data/jawa_bali_nt_sulawesi_site_model.xml'#jawa_site_model.xml'
outfile = 'data/jawa_bali_nt_sulawesi_site_model_full.xml'

#locations_file = 'data/jawa_sites_test.csv'
#sites_file = 'data/jawa_site_model.xml'
#outfile = 'data/jawa_site_model_test.xml'
#if myid == 0:
sites_data = np.genfromtxt(locations_file, delimiter=',')
#print sites_data
# Split sources
list_length = len(sites_data) / (proc)
print list_length
if (len(sites_data) % proc) > 0:
    list_length +=1
pt_list = list(chunks(sites_data, list_length))
#print pt_list

for i in range(0, len(pt_list), 1):
    if i % proc == myid:
        run = "%03d" % i
        # Apply geometrical filtering                                                                        
        print 'Building site model for run %s' % run
        tmp_outfile = outfile[:-4] + '_' + str(run) + '.xml'
        build_site_col(pt_list[i], sites_file, filename=tmp_outfile)
pypar.barrier()
if myid == 0:
    outlines = []
    tmp_pt_source_filename_list = []
    tmp_pt_source_list = []
    # Now combine into one file
    for j in range(0, len(pt_list), 1):
        tmp_pt_filename = outfile[:-4] + '_%03d.xml' % j 
        tmp_pt_source_filename_list.append(tmp_pt_filename)
    for tmp_pt_source_file in tmp_pt_source_filename_list:
        f_in = open(tmp_pt_source_file, 'r')
        for line in f_in.readlines():
            if line.lstrip().startswith('<site '):
                outlines.append(line)
        f_in.close()
    f_out = open(outfile, 'w')
    # write header
    f_out.write('<?xml version="1.0" encoding="utf-8"?>\n')
    f_out.write('<nrml\n')
    f_out.write('xmlns="http://openquake.org/xmlns/nrml/0.4"\n')
    f_out.write('xmlns:gml="http://www.opengis.net/gml"\n')
    f_out.write('  <siteModel>\n')
    #outlines = [oxml + '\n' for oxml in outlines]
    f_out.writelines(outlines)
    f_out.write('  </siteModel>\n')
    f_out.write('</nrml>')
    ss = int(pypar.time() - t0)
    h = ss / 3600
    m = (ss % 3600) / 60
    s = (ss % 3600) % 60
    print "--------------------------------------------------------"
    print 'P0: Total time (%i seconds): %s:%s:%s (hh:mm:ss)' % (ss,
                                                                string.zfill(h, 2),
                                                                string.zfill(m, 2),
                                                                string.zfill(s,2))
    print "--------------------------------------------------------"
pypar.finalize()
