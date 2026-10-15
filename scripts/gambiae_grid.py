import pandas as pd, numpy as np
from math import ceil
import h3
import geopandas as gpd
from pyproj import Geod, Transformer
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon, Point, LineString
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
                     outpath=None # where to save shapefile to
                     ):
    center_h3 = h3.latlng_to_cell(lat, lon, resolution)
    # get the hexagons within the specified radius
    edge_length = h3.edge_length(h3.origin_to_directed_edges(center_h3)[0], unit='km')
    k_distance = ceil(radius / edge_length)
    print(radius, k_distance, edge_length)
    hexagons = h3.grid_disk(center_h3, k_distance)

    # convert hexagons to geoJSON geometries
    polygons = []
    cells = []
    for cell in hexagons:
        # polygon of cell
        boundary = h3.cell_to_boundary(cell)
        # flip for Shapely
        boundary = [(lon, lat) for lat, lon in boundary]
        polygon = Polygon(boundary)
        polygons.append(polygon)
        cells.append(cell)

    # create a GeoDataFrame
    gdf = gpd.GeoDataFrame({'h3_index': cells, 'geometry': polygons})
    gdf.set_crs(epsg=4326, inplace=True)  # set the coordinate reference system to WGS84
    # save to shapefile
    if outpath:
        local_crs = (
            f'+proj=aeqd +lat_0={lat} +lon_0={lon} '
            '+datum=WGS84 +units=m +no_defs'
        )
        transformer = Transformer.from_crs('EPSG:4326', local_crs, always_xy=True)
        cell_set = set(cells)
        center_coords = {}
        projected_centers = {}
        for cell in cells:
            cell_lat, cell_lon = h3.cell_to_latlng(cell)
            center_coords[cell] = (cell_lon, cell_lat)
            projected_centers[cell] = np.asarray(
                transformer.transform(cell_lon, cell_lat)
            )

        frame_polygons = []
        frame_cells = []
        seen_triangles = set()
        for cell in cells:
            neighbors = [
                neighbor
                for neighbor in h3.grid_disk(cell, 1)
                if neighbor != cell and neighbor in cell_set
            ]
            neighbors.sort(
                key=lambda neighbor: np.arctan2(
                    projected_centers[neighbor][1] - projected_centers[cell][1],
                    projected_centers[neighbor][0] - projected_centers[cell][0],
                )
            )
            if len(neighbors) < 3:
                continue

            for i, first in enumerate(neighbors):
                second = neighbors[(i + 1) % len(neighbors)]
                if not h3.are_neighbor_cells(first, second):
                    continue
                triangle_key = tuple(sorted((cell, first, second)))
                if triangle_key in seen_triangles:
                    continue
                seen_triangles.add(triangle_key)
                frame_polygons.append(
                    Polygon([center_coords[cell], center_coords[first], center_coords[second]])
                )
                frame_cells.append(cell)

        frame_gdf = gpd.GeoDataFrame(
            {'h3_index': frame_cells, 'geometry': frame_polygons},
            crs='EPSG:4326',
        )
        frame_gdf.to_file(f'{outpath}_grid_resolution_{resolution}.shp', driver='ESRI Shapefile')
        coords = np.asarray([center_coords[cell] for cell in cells])
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