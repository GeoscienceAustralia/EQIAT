from mpl_toolkits.basemap import Basemap, maskoceans
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from osgeo import gdal
import numpy as np
from numpy import linspace
from numpy import meshgrid

m = Basemap(projection='tmerc', 
              lat_0=0, lon_0=110,
              llcrnrlon=94.0,
              llcrnrlat=-10.,
              urcrnrlon=141.,
              urcrnrlat=6.)
m.drawcoastlines(linewidth=0.5,color='k')
m.drawcountries(color='0.2')
m.drawstates(color='0.2')

m.drawparallels(np.arange(-90.,90.,2), labels=[1,0,0,0],
                fontsize=10, dashes=[2, 2], color='0.5',
                linewidth=0.5)
m.drawmeridians(np.arange(0.,360.,2), labels=[0,0,0,1],
                fontsize=10, dashes=[2, 2], color='0.5',
                linewidth=0.5)

infile = 'outputs/scenario_gmf_loc_mmi_1847_ChiouYoungs2008_gridded_clip.tif'

ds = gdal.Open(infile)
data = ds.GetRasterBand(1).ReadAsArray()
#data = np.flipud(data)

# get the edge coordinates and add half the resolution                                               
# to go to center coordinates                                 
gt = ds.GetGeoTransform()
proj = ds.GetProjection()
xres = gt[1]
yres = gt[5]                                       
xmin = gt[0] + xres * 0.5
xmax = gt[0] + (xres * ds.RasterXSize) - xres * 0.5
ymin = gt[3] + (yres * ds.RasterYSize) + yres * 0.5
ymax = gt[3] - yres * 0.5

print xmin,xmax
print ymin,ymax
#ds = None                                                                                           

# create a grid of xy coordinates in the original projection                                         
xy_source = np.mgrid[xmin:xmax+xres:xres, ymax+yres:ymin:yres]
print xy_source
print xy_source.shape

x = linspace(0, m.urcrnrx, data.shape[1])
y = linspace(0, m.urcrnry, data.shape[0])

xy = meshgrid(x, y)
print xy_source[0,:,0]
print xy_source[1][0]
xy=meshgrid(xy_source[0,:,0], xy_source[1][0])
print xy
print xy[0].shape
print xy[1].shape
maskdata = maskoceans(xy[0], xy[1], data)
maskdata = np.ma.masked_where(maskdata==0.0, maskdata)
data = maskdata

m.contourf(xy[0], xy[1], data, latlon=True)
plt.savefig('outputs/test.png')
