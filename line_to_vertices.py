""" Function to get vertices from a shapefile"""

import ogr

shpfile = 'outputs/rupture_scenario_1867_CampbellBozorgnia2008()_Mw6.35_110.725_-7.939_20.00km.shp'

def line_to_vertices(shpfile):
    """Extract vertices from a line
    """
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shpfile, 0)
    layer = dataSource.GetLayer(0)
    print(layer.GetFeatureCount())
    feature = layer.GetFeature(0)
    print(feature.GetFieldCount())
    print(feature.items())
    geometry = feature.GetGeometryRef()
    print(geometry.GetGeometryCount())
    line = geometry.GetGeometryRef(0)
    print(line.GetPoint(0))
    print(line.GetPoint(-1))

if __name__=="__main__":
    line_to_vertices(shpfile)
