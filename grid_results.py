"""Interpolate csv output results to raster format
Jonathan Griffin
Geoscience Australia June 2017
"""

import os, sys
import numpy as np
from scipy.interpolate import griddata
import ogr, gdal, osr

def grid_results(infile, resolution = 0.01):
    """Read csv file, interpolate and convert to raster
    """
    data = np.genfromtxt(infile, delimiter=',')
    max_lon = max(data[:,0])
    min_lon = min(data[:,0])
    max_lat = max(data[:,1])
    min_lat = min(data[:,1])
    #print max_lon, min_lon, max_lat, min_lat
    xi = np.arange(min_lon, max_lon, resolution)
    yi = np.arange(min_lat, max_lat, resolution)
    XI,YI = np.meshgrid(xi,yi)
    xsize = len(xi)
    ysize = len(yi)

    print 'Interpolating results'
    gridded_results = griddata((data[:,0],data[:,1]),data[:,2],(XI,YI),method='linear')
    #print gridded_results
    outfile = infile.rstrip('().csv') + '_gridded.tif'
    print 'Writing gridded data to %s' % outfile
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(outfile, xsize, ysize, 1, gdal.GDT_Float32)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    gt = [(min_lon - (resolution/2)), resolution, 0, 
          (min_lat - (resolution/2)), 0, resolution]
    ds.SetGeoTransform(gt)
    outband=ds.GetRasterBand(1)
    outband.SetStatistics(np.min(gridded_results), np.max(gridded_results), np.average(gridded_results), np.std(gridded_results))
    outband.WriteArray(gridded_results)

   # outfile = infile.rstrip('.csv') + '_gridded.csv'
   # np.savetxt(outfile, gridded_results)

infile = 'outputs/scenario_gmf_loc_mmi_1699megathrust_AtkinsonBoore2003SInter().csv'
grid_results(infile)
