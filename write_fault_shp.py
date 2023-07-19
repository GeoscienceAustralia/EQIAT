"""Write a fault or rupture geometry to shapefile
"""

#import gdal, osr, ogr
from osgeo import ogr
import os 
from numpy import mean

def fault2shp(corner_lons, corner_lats, output_shp, corner_depths=None, vertice_array=False):
    """Function for writing a fault geometry to a shapefile
    """
    
    # Create a Polygon from the extent tuple
    ring = ogr.Geometry(ogr.wkbLinearRing)
    #for i in range(len(corner_lons)):
    # need to get in right order
    if vertice_array:
        # Assume corner_lons, corner_lats are 2 1D array
        # giving the corrdinates of the polygon boundary
        for i in range(len(corner_lons)):
            ring.AddPoint(corner_lons[i], corner_lats[i])
        ring.AddPoint(corner_lons[0], corner_lats[0]) # close polygon
    else:
        ring.AddPoint(corner_lons[0],corner_lats[0])
        ring.AddPoint(corner_lons[1],corner_lats[1])
        ring.AddPoint(corner_lons[3],corner_lats[3])
        ring.AddPoint(corner_lons[2],corner_lats[2])
        ring.AddPoint(corner_lons[0],corner_lats[0]) # close polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    drv = ogr.GetDriverByName('ESRI Shapefile') 
    # Remove output shapefile if it already exists
    if os.path.exists(output_shp):
        drv.DeleteDataSource(output_shp)
    # Create the output shapefile
    outDataSource = drv.CreateDataSource(output_shp)
    outLayer = outDataSource.CreateLayer("Fault_geom", geom_type=ogr.wkbPolygon)
    # Add an ID field
    idField = ogr.FieldDefn("id", ogr.OFTInteger)
    outLayer.CreateField(idField)
    # Add a depth field
    depthField = ogr.FieldDefn("mean_depth", ogr.OFTReal)
    outLayer.CreateField(depthField)
    # Create the feature and set values
    featureDefn = outLayer.GetLayerDefn()
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(poly)
    feature.SetField("id", 1)
    feature.SetField("mean_depth", mean(corner_depths))
    outLayer.CreateFeature(feature)
    feature = None
    # Save and close
    outDataSource = None
    drv = None

    # Now write upper trace to line shapefile    
    line = ogr.Geometry(ogr.wkbLineString)
    corner_depths = list(corner_depths)
    min_dep = min(corner_depths)
    min_dep_index = corner_depths.index(min(corner_depths))
    print(min_dep_index)
    line.AddPoint(corner_lons[min_dep_index],corner_lats[min_dep_index])
    corner_depths[min_dep_index] = 1e10
    min_dep_2 = min(corner_depths)
    min_dep_index_2 = corner_depths.index(min(corner_depths))
    print(min_dep_index_2)
    line.AddPoint(corner_lons[min_dep_index_2],corner_lats[min_dep_index_2])
    corner_depths[min_dep_index] = min_dep
    print(min_dep, min_dep_2)
    mean_upper_depth = mean([min_dep, min_dep_2])
    print(mean_upper_depth)
    drv = ogr.GetDriverByName('ESRI Shapefile') 
    output_shp = output_shp.rstrip('.shp') + '_upper_edge.shp'
    # Remove output shapefile if it already exists
    if os.path.exists(output_shp):
        drv.DeleteDataSource(output_shp)
    # Create the output shapefile
    outDataSource = drv.CreateDataSource(output_shp)
    outLayer = outDataSource.CreateLayer("Fault_geom", geom_type=ogr.wkbLineString)   
    # Add a depth field
    depthField = ogr.FieldDefn("mean_depth", ogr.OFTReal)
    outLayer.CreateField(depthField)
    # Create the feature and set values
    featureDefn = outLayer.GetLayerDefn()
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(line)
    feature.SetField("mean_depth", mean_upper_depth)
    outLayer.CreateFeature(feature)
    feature = None
    # Save and close
    outDataSource = None
    drv = None

if __name__=="__main__":
    corner_lons =  [ 110.80301631, 110.14964724, 110.86646785, 110.21405031]
    corner_lats = [-9.39397613, -9.31704483, -8.86453769, -8.78772005]
    corner_depths = [ 32.81350292, 32.81350292, 74.01549338, 74.01549338]
    fault2shp(corner_lons, corner_lats, 'test_fault.shp', corner_depths)
    
