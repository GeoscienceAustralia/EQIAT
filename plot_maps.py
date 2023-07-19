import sys
import re
from mpl_toolkits.basemap import Basemap, maskoceans
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors, colorbar
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from osgeo import gdal
import shapefile
import numpy as np
from numpy import linspace
from numpy import meshgrid
from collections import OrderedDict
from scipy import interpolate
from adjustText import adjust_text # Small package to improve label locations


mmi_obs_file = sys.argv[1]# 'data/1847HMMI.txt'
#infile = 'outputs/scenario_gmf_loc_mmi_1847_ChiouYoungs2008_gridded_clip.tif'
infile = sys.argv[2] #r'outputs/scenario_gmf_loc_mmi_1847_ChiouYoungs2008().csv'
bbox = sys.argv[3] #'104/116/-10/-5'
try:
    shpfile = sys.argv[4]
except:
    shpfile = None

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
        for r in list(roman.keys()):
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
"""
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
"""
clevs = np.arange(0.5,9.5,1.0)
cmap = plt.get_cmap('jet')
#norm = BoundaryNorm(clevs, ncolors=cmap.N, clip=True)
#m.imshow(data, cmap=cmap, vmin=0, vmax=9, zorder=0, latlon=True)

data = np.genfromtxt(infile, delimiter=',')
print(data[:,0])
minlon = np.min(data[:,0])
maxlon = np.max(data[:,0])
minlat = np.min(data[:,1])
maxlat = np.max(data[:,1])
xy = np.mgrid[minlon:maxlon:0.02,minlat:maxlat:0.02]
xx,yy=meshgrid(xy[0,:,0], xy[1][0])
griddata = interpolate.griddata((data[:,0], data[:,1]), data[:,2], (xx,yy), method='nearest')
# Mask Oceans
maskdata = maskoceans(xx, yy, griddata)
m.drawlsmask(resolution='i', grid=2.5,ocean_color='w')
m.contourf(xx, yy, maskdata, clevs, latlon=True)
csm = m.contour(xx, yy, maskdata, clevs, latlon=True, colors='k')

#m.contourf(xy[0], xy[1], data, clevs, latlon=True)
#csm = m.contour(xy[0], xy[1], data, clevs, latlon=True, colors='k')
#plt.clabel(csm, inline=1, fontsize=10, fmt='%0.2f')
#m.pcolormesh(xy[0], xy[1], data, latlon=True)

# Now plot fault shapefile, if it exists
# Note that we need to remove . from shapefile names for them to work with this
if shpfile is not None:
    print(shpfile[:-4])
    m.readshapefile(shpfile[:-4],'rupture', drawbounds=True)
    patches   = []
    for info, shape in zip(m.rupture_info, m.rupture):
        patches.append( Polygon(np.array(shape), True) )  
    ax.add_collection(PatchCollection(patches, facecolor= 'm', edgecolor='k', linewidths=2., zorder=2))

# Now add historical points on top                                                                            
mmi_labels = []
for obs in mmi_obs[:,2]:
    mmi_labels.append(write_roman(int(obs)))
m.scatter(mmi_obs[:,0], mmi_obs[:,1], c=mmi_obs[:,2], cmap=cmap, 
          vmin=0.5, vmax=8.5, s=40, latlon=True)
texts = []
for label, x, y in zip(mmi_labels, mmi_obs[:,0], mmi_obs[:,1]):
    x,y =  m(x,y)
    texts.append(plt.text(x,y,label))
#    texts.append(plt.annotate(label, xy=(x,y),
#                          xycoords='data', fontsize=11))
adjust_text(texts, only_move='xy',
            arrowprops=dict(arrowstyle="->",
                            color='k', lw=0.5))
filename_parts = infile.rstrip('().csv').split('_')
year = filename_parts[4][:4]
print('year', year)
# Some string splitting to get the GMM name
gmm = filename_parts[5]
a=re.sub('([a-z])([1-9])', r'\1 \2', gmm)
b = re.sub('([a-z])([A-Z])', r'\1 \2', a) 
c = re.sub('([1-9])([A-Z])', r'\1 \2', b)
d = re.sub('([A-Z])([A-Z])', r'\1 \2', c).split()
d = [s.replace('Et', 'et') for s in d]
d = [s.replace('Al','al') for s in d]
d = [s.replace('S','Subduction') for s in d]
d = [s.replace('Subductionlab','Slab') for s in d]
d = [s.replace('Inter','Interface') for s in d]
if 'et' not in d:
    d.insert(1, '&')
gmm_name = ' '.join(d)
title = year + '\n' + gmm_name
plt.title(title)
ax = plt.gca()
"""
for child in ax.get_children():
    print type(child)
    if isinstance(child, matplotlib.patches.FancyArrowPatch):
        print 'here'
        sys.exit()
    if isinstance(child, matplotlib.text.Annotation):
        print child.arrowprops
        print child._arrow_relpos
        print child.arrow_patch.get_arrowstyle()
        print child.arrow_patch.get_patch_transform().get_matrix()
        print child.arrow_patch.patchA
        print child.arrow_patch.patchB
#        print child.arrow_patch.get_xy()
#print dir(child)
#        print child.arrow
        print dir(child.arrow_patch)
        print child.arrow_patch._posA_posB
        print child.arrow_patch.get_path()
        print child.arrow_patch.get_path_in_displaycoord()
        #print child.arrow_patch._arrow_relpos
        print child.arrow_patch.shrinkA
        print child.arrow_patch.shrinkB
"""
#    plt.annotate(
#        label,
#        xy=(x, y), xytext=(-2, 2), xycoords='data', fontsize=11,
#        textcoords='offset points', ha='right', va='bottom')#,
#        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
#        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
plt.gcf().subplots_adjust(bottom=0.1)
cax = figure.add_axes([0.34,0.05,0.33,0.035]) # setup colorbar axes.
#bounds = np.linspace(0,9,10)
norm = colors.BoundaryNorm(clevs, cmap.N)
#norm = colors.Normalize(vmin=0, vmax=9)
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')
cb.ax.get_xaxis().set_ticks([])
cb.set_label('MMI', fontsize=16)
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

figfilename = infile.rstrip('().csv') + '_fig.png'
plt.savefig(figfilename)
