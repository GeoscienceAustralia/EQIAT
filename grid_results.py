"""Interpolate csv output results to raster format
Jonathan Griffin
Geoscience Australia June 2017
"""

import os, sys
from subprocess import call
from glob import glob
import numpy as np
from scipy.interpolate import griddata
import ogr, gdal, osr

def grid_results(infile, resolution = 0.01, clip_shp = None, 
                 overwrite=True, contour=False):
    """Read csv file, interpolate and convert to raster
    """
    outfile = infile.rstrip('().csv') + '_gridded.tif'
   # if not overwrite:
    if os.path.isfile(outfile):
        if not overwrite:
            print 'Not creating file %s as already exists' % outfile
            print 'To re-create file (e.g if inputs changed) set overwrite=True)'
            return
        else:
            try:
                os.remove(outfile)
                os.remove((outfile.rstrip('.tif') + '_clip.tif'))
            except:
                pass
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
    #outfile = infile.rstrip('().csv') + '_gridded.tif'
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
    # Need to close output dataset before we can do clipping
    ds = None
    # now clip by shapefile
    if clip_shp is not None:
        clipfile = outfile.rstrip('.tif') + '_clip.tif'
        cmd = ['gdalwarp',
               '-cutline',
               clip_shp,
               '-crop_to_cutline',
               '-dstalpha',
               outfile,
               clipfile]
        print cmd
        call(cmd, shell=False)
    if contour is True:
        cmd = 'gdal_contour -i 1 -off 0.5 %s %s.shp' % (outfile,  outfile.rstrip('.tif'))
        print cmd
        call(cmd, shell=True)
        cmd = 'gdal_contour -i 1 -off 0.5 %s %s.shp' % (clipfile, clipfile.rstrip('.tif'))
        print cmd
        call(cmd, shell=True)
   # outfile = infile.rstrip('.csv') + '_gridded.csv'
   # np.savetxt(outfile, gridded_results)

filelist = glob('outputs/*mmi*1780*sq*Chiou*.csv')
#infile = 'outputs/scenario_gmf_loc_mmi_1699megathrust_AtkinsonBoore2003SInter().csv'
clipping_shapefile = '/short/n74/jdg547/eq_hazmap_tests/data/IDN_Admin/IDN_adm0.shp'
for infile in filelist:
    grid_results(infile, resolution=0.01, clip_shp = clipping_shapefile,
                 overwrite=True, contour=True)
