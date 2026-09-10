import pandas as pd, numpy as np
import h3
import geopandas as gpd
from pyproj import Geod
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon, Point, mapping, LineString
from shapely.affinity import scale

def generate_h3_grid(lat, # latitude of the center point 
                     lon, # longitude of the center point
                     radius, # radius in KM
                     resolution # 0-15, higher resolution means smaller hexagons
                     ):
    center_h3 = h3.latlng_to_cell(lat, lon, resolution)
    # get the hexagons within the specified radius
    k_distance = radius / (h3.edge_length(h3.origin_to_directed_edges(center_h3)[0], unit='km'))
    print(radius, k_distance)
    hexagons = h3.grid_disk(center_h3, k_distance)

    # convert hexagons to geoJSON geometries
    polygons = []
    cells = []
    for cell in hexagons:
        boundary = h3.cell_to_boundary(cell)
        # flip for Shapely
        boundary = [(lon, lat) for lat, lon in boundary]
        polygons.append(Polygon(boundary))
        cells.append(cell)

    # create a GeoDataFrame
    gdf = gpd.GeoDataFrame({'h3_index': cells, 'geometry': polygons})
    gdf.set_crs(epsg=4326, inplace=True)  # set the coordinate reference system to WGS84
    return gdf

# read in sampling locations
df = pd.read_csv('data/vo_agam_release/gambiae_latlon.txt', header=None, names=['latitude', 'longitude'], dtype='float')
df = df.dropna() # drop any rows with missing values
# centroid of the sampling locations
lat_c = df['latitude'].mean()
lon_c = df['longitude'].mean()
point = Point(lon_c, lat_c)
# calculate the maximum distance from the centroid to any sampling location
geod = Geod(ellps="WGS84")
max_d = df.apply(lambda row: geod.geometry_length(LineString([Point(row['longitude'], row['latitude']), point])), axis=1).max() * 0.001
# generate grid of hexagons around the centroid with a radius equal to the maximum distance
grid = generate_h3_grid(lat_c, lon_c, max_d, resolution=2)
coords = [mapping(grid.geometry[i])['coordinates'][0][:-1] for i in range(len(grid))]
coords = [coord for sublist in coords for coord in sublist]
coords = [tuple(list(coord)[::-1]) for coord in coords] # flip to lat, lon
coords = np.array(list(set(coords)))
# remove any coordinates that are outside the bounding box of the sampling locations
hull = ConvexHull(df[['longitude', 'latitude']].values)
hull = Polygon(df[['longitude', 'latitude']].values[hull.vertices])
hull = scale(hull, xfact=1.1, yfact=1.1, origin='centroid') # scale the hull by 10% to include points just outside the convex hull
coords = np.array([coord for coord in coords if hull.contains(Point(coord[1], coord[0]))])

np.savetxt('data/vo_agam_release/gambiae_grid_coords.csv', coords, delimiter=',', comments='', fmt='%1.6f')