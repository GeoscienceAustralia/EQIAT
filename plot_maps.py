from mpl_toolkits.basemap import Basemap, maskoceans
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors, colorbar
from osgeo import gdal
import numpy as np
from numpy import linspace
from numpy import meshgrid
from collections import OrderedDict

mmi_obs_file = 'data/1847HMMI.txt'
infile = 'outputs/scenario_gmf_loc_mmi_1847_ChiouYoungs2008_gridded_clip.tif'
bbox = '104/116/-10/-5'

# function fro roman numerals
def write_roman(num):

    roman = OrderedDict()
    roman[1000] = "M"
    roman[900] = "CM"
    roman[500] = "D"
    roman[400] = "CD"
    roman[100] = "C"
    roman[90] = "XC"
    roman[50] = "L"
    roman[40] = "XL"
    roman[10] = "X"
    roman[9] = "IX"
    roman[5] = "V"
    roman[4] = "IV"
    roman[1] = "I"

    def roman_num(num):
        for r in roman.keys():
            x, y = divmod(num, r)
            yield roman[r] * x
            num -= (r * x)
            if num > 0:
                roman_num(num)
            else:
                break

    return "".join([a for a in roman_num(num)])

# Read observation data
mmi_obs = np.genfromtxt(mmi_obs_file)

# Setup map
figure, ax = plt.subplots(1,figsize=(16,8))
bbox = bbox.split('/')
minlon = float(bbox[0])
maxlon = float(bbox[1])
minlat = float(bbox[2])
maxlat = float(bbox[3])
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

m.drawparallels(np.arange(-90.,90.,2), labels=[1,0,0,0],
                fontsize=10, dashes=[2, 2], color='0.5',
                linewidth=0.5)
m.drawmeridians(np.arange(0.,360.,2), labels=[0,0,0,1],
                fontsize=10, dashes=[2, 2], color='0.5',
                linewidth=0.5)

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
#x = linspace(0, m.urcrnrx, data.shape[1])
#y = linspace(0, m.urcrnry, data.shape[0])
#xy = meshgrid(x, y)
xy=meshgrid(xy_source[0,:,0], xy_source[1][0])
maskdata = maskoceans(xy[0], xy[1], data)
data = np.ma.masked_where(maskdata==0, maskdata)
clevs = np.arange(0.5,9.5,1.0)
cmap = plt.get_cmap('jet')
#norm = BoundaryNorm(clevs, ncolors=cmap.N, clip=True)
#m.imshow(data, cmap=cmap, vmin=0, vmax=9, zorder=0, latlon=True)
m.contourf(xy[0], xy[1], data, clevs, latlon=True)
csm = m.contour(xy[0], xy[1], data, clevs, latlon=True, colors='k')
#plt.clabel(csm, inline=1, fontsize=10, fmt='%0.2f')
#m.pcolormesh(xy[0], xy[1], data, latlon=True)

# Now add historical points on top                                                                            
mmi_labels = []
for obs in mmi_obs[:,2]:
    mmi_labels.append(write_roman(int(obs)))
m.scatter(mmi_obs[:,0], mmi_obs[:,1], c=mmi_obs[:,2], cmap=cmap, 
          vmin=0.5, vmax=8.5, latlon=True)
for label, x, y in zip(mmi_labels, mmi_obs[:,0], mmi_obs[:,1]):
    x,y =  m(x,y)
    plt.annotate(
        label,
        xy=(x, y), xytext=(-2, 2), xycoords='data', fontsize=11,
        textcoords='offset points', ha='right', va='bottom')#,
#        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
#        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
plt.gcf().subplots_adjust(bottom=0.1)
cax = figure.add_axes([0.34,0.05,0.33,0.035]) # setup colorbar axes.
#bounds = np.linspace(0,9,10)
norm = colors.BoundaryNorm(clevs, cmap.N)
#norm = colors.Normalize(vmin=0, vmax=9)
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')
cb.ax.get_xaxis().set_ticks([])
#for j, lab in enumerate(['$0$','$1$','$2$','$>3$']):
clevs = clevs[1:]
labs = []
for j, lab in enumerate(clevs):
    lab = '%i' % lab
    if lab=='8':
        lab = '>8'
    labs.append(lab)
    cb.ax.text((j+0.5)/8.0, 0, lab, ha='center', va='bottom')
#print labs
#cb.ax.set_xticklabels(labs)

figfilename = infile.rstrip('.tif') + '_fig.png'
plt.savefig(figfilename)
