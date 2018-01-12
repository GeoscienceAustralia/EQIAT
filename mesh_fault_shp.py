"""Helper script to convert mesh points to a polygon 
"""

import numpy as np
from write_fault_shp import fault2shp 
infile = 'outputs/rupture_mesh_1818_YoungsEtAl1997SInter_Mw6p75_118p882_-8p102_13p62km.txt'

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def mesh2shp(infile):
    """
    """
    output_shp = infile[:-3] + 'shp'
    output_shp = output_shp.replace('mesh', 'scenario')

    file_length = file_len(infile)
    print file_length
    
    # Work out lines that have the required data
    top_lons = 0
    bottom_lons = int(file_length/3 - 1)
    top_lats = int(file_length/3)
    bottom_lats = int((file_length/3)*2 - 1)
    top_depths = int((file_length/3)*2)
    bottom_depths = int(file_length-1)
    print top_lons,bottom_lons,top_lats,bottom_lats,top_depths,bottom_depths
    
    lons = []
    lats = []
    depths = []
    with open(infile) as f:
        for i, l in enumerate(f):
            if i == top_lons:
                lons.append([float(x) for x in l.split(',')])
            if i == bottom_lons:
                tmp_lons = ([float(x) for x in l.split(',')])
                print tmp_lons
                tmp_lons = list(reversed(tmp_lons))
                lons.append(tmp_lons)

            if i == top_lats:
                lats.append([float(x) for x in l.split(',')])
            if i == bottom_lats:
                tmp_lats = ([float(x) for x in l.split(',')])
                print tmp_lats
                tmp_lats = list(reversed(tmp_lats))
                lats.append(tmp_lats)

            if i == top_depths:
                depths.append([float(x) for x in l.split(',')])
            if i == bottom_depths:
                tmp_depths = ([float(x) for x in l.split(',')])
                print tmp_depths
                tmp_depths = list(reversed(tmp_depths))
                depths.append(tmp_depths)
    lons = np.array(lons).flatten()
    lats = np.array(lats).flatten()
    depths = np.array(depths).flatten()
    print lons
    print lats
    print depths

    fault2shp(lons, lats, output_shp, corner_depths = depths, vertice_array=True)

if __name__ == "__main__":
    infiles = [
    'outputs/rupture_mesh_1699megathrust_AtkinsonBoore2003SInter_Mw8p05_107p292_-8p084_42p96km.txt',
    'outputs/rupture_mesh_1699megathrust_YoungsEtAl1997SInter_Mw8p15_107p001_-7p880_39p11km.txt',
    'outputs/rupture_mesh_1699megathrust_ZhaoEtAl2006SInter_Mw8p45_107p105_-8p337_32p28km.txt',
    'outputs/rupture_mesh_1815_AtkinsonBoore2003SInter_Mw7p45_115p084_-7p824_9p08km.txt',
    'outputs/rupture_mesh_1815_YoungsEtAl1997SInter_Mw7p45_114p634_-7p768_9p08km.txt',
    'outputs/rupture_mesh_1815_ZhaoEtAl2006SInter_Mw7p45_114p904_-7p798_9p08km.txt',
    'outputs/rupture_mesh_1818_AtkinsonBoore2003SInter_Mw6p95_118p925_-8p061_11p35km.txt',
    'outputs/rupture_mesh_1818_YoungsEtAl1997SInter_Mw6p75_118p882_-8p102_13p62km.txt',
    'outputs/rupture_mesh_1818_ZhaoEtAl2006SInter_Mw7p05_119p016_-8p058_11p35km.txt',
    'outputs/rupture_mesh_1820_AtkinsonBoore2003SInter_Mw7p65_119p059_-8p017_9p08km.txt',
    'outputs/rupture_mesh_1820_AtkinsonBoore2003SInter_Mw8p35_124p707_-7p920_9p08km.txt',
    'outputs/rupture_mesh_1820_YoungsEtAl1997SInter_Mw7p85_118p288_-8p057_9p08km.txt',
    'outputs/rupture_mesh_1820_YoungsEtAl1997SInter_Mw8p35_122p464_-8p251_9p08km.txt',
    'outputs/rupture_mesh_1820_ZhaoEtAl2006SInter_Mw7p85_118p651_-8p036_9p08km.txt',
    'outputs/rupture_mesh_1820_ZhaoEtAl2006SInter_Mw8p35_122p464_-8p251_9p08km.txt']
    for infile in infiles:
        mesh2shp(infile)
