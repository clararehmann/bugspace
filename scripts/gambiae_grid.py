import pandas as pd, numpy as np
import h3
import geopandas as gpd
from pyproj import Geod
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon, Point, mapping, LineString
from shapely.affinity import scale

def get_sampling_radius(locs):
    # centroid of the sampling locations
    lat_c = locs['latitude'].mean()
    lon_c = locs['longitude'].mean()
    point = Point(lon_c, lat_c)
    # calculate the maximum distance from the centroid to any sampling location
    geod = Geod(ellps="WGS84")
    max_d = locs.apply(lambda row: geod.geometry_length(LineString([Point(row['longitude'], row['latitude']), point])), axis=1).max() * 0.001
    return lat_c, lon_c, max_d

def generate_h3_grid(lat, # latitude of the center point 
                     lon, # longitude of the center point
                     radius, # radius in KM
                     resolution, # 0-15, higher resolution means smaller hexagons
                     outpath # where to save shapefile to
                     ):
    center_h3 = h3.latlng_to_cell(lat, lon, resolution)
    # get the hexagons within the specified radius
    k_distance = radius / (h3.edge_length(h3.origin_to_directed_edges(center_h3)[0], unit='km'))
    print(radius, k_distance, h3.edge_length(h3.origin_to_directed_edges(center_h3)[0], unit='km'))
    hexagons = h3.grid_disk(center_h3, k_distance)

    # convert hexagons to geoJSON geometries
    polygons = []
    cells = []
    vertexes = []
    for cell in hexagons:
        # polygon of cell
        boundary = h3.cell_to_boundary(cell)
        # flip for Shapely
        boundary = [(lon, lat) for lat, lon in boundary]
        polygons.append(Polygon(boundary))
        cells.append(cell)

    # create a GeoDataFrame
    gdf = gpd.GeoDataFrame({'h3_index': cells, 'geometry': polygons})
    gdf.set_crs(epsg=4326, inplace=True)  # set the coordinate reference system to WGS84
    # save to shapefile
    gdf.to_file(f'{outpath}_grid_resolution_{resolution}.shp', driver='ESRI Shapefile')
    # save to lat, lon csv
    coords = [mapping(gdf.geometry[i])['coordinates'][0][:-1] for i in range(len(gdf))]
    coords = [coord for sublist in coords for coord in sublist]
    coords = [tuple(list(coord)) for coord in coords] # flip to lat, lon
    coords = np.array(list(set(coords)))
    np.savetxt(f'{outpath}_grid_resolution_{resolution}_coordinates.csv', coords, delimiter=',', comments='', fmt='%1.6f')
    return gdf

def generate_hull(locs, outpath, x_fact=1.1, y_fact=1.1):
    hull = ConvexHull(locs[['longitude', 'latitude']].values)
    hull = Polygon(locs[['longitude', 'latitude']].values[hull.vertices])
    hull = scale(hull, xfact=x_fact, yfact=y_fact, origin='centroid')
    lon, lat = hull.exterior.coords.xy
    hull = np.asarray(tuple(zip(lon,lat)))
    np.savetxt(f'{outpath}_hull.csv', hull, delimiter=',', comments='', fmt='%1.6f')
    return

# read in sampling locations
#df = pd.read_csv('data/test_lonlat.txt', header=None, names=['longitude', 'latitude'], dtype='float')
#df = df.dropna() # drop any rows with missing values
#lat_c, lon_c, max_d = get_sampling_radius(df)
# generate grid of hexagons around the centroid with a radius equal to the maximum distance
#grid = generate_h3_grid(lat_c, lon_c, max_d, resolution=1, outpath='data/test')
#coords = [mapping(grid.geometry[i])['coordinates'][0][:-1] for i in range(len(grid))]
#coords = [coord for sublist in coords for coord in sublist]
#coords = [tuple(list(coord)[::-1]) for coord in coords] # flip to lat, lon
#coords = np.array(list(set(coords)))

# remove any coordinates that are outside the bounding box of the sampling locations
#hull = ConvexHull(df[['longitude', 'latitude']].values)
#hull = Polygon(df[['longitude', 'latitude']].values[hull.vertices])
#hull = scale(hull, xfact=10, yfact=10, origin='centroid') # scale the hull by 10% to include points just outside the convex hull
#coords = np.array([coord for coord in coords if hull.contains(Point(coord[1], coord[0]))])
# remove coordinates not on land
#earth = gpd.read_file('data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp')
#earth = earth.dissolve(by='CONTINENT')
#africa = earth.loc['Africa','geometry']
#coords = np.array([coord for coord in coords if africa.contains(Point(coord[1], coord[0]))])
#coords = pd.DataFrame({'lat':[coord[0] for coord in coords], 'lon':[coord[1] for coord in coords]})
#coords['geometry'] = coords.apply(lambda x: Point((x.lon, x.lat)), axis=1)
#coords = gpd.GeoDataFrame(coords, geometry='geometry')
#coords.set_crs("EPSG:4326", inplace=True)
#coords.to_file('data/vo_agam_release/gambiae_grid.shp', driver='ESRI Shapefile')
#np.savetxt('data/vo_agam_release/gambiae_grid_coords.csv', coords, delimiter=',', comments='', fmt='%1.6f')
# save coordinates of hull
#lon, lat = hull.exterior.coords.xy
#hull = np.asarray(tuple(zip(lon, lat)))
#np.savetxt('data/test_hull_coords.csv', hull, delimiter=',', comments='', fmt='%1.6f')