import pandas as pd
import h3
import geopandas as gpd
from shapely.geometry import Polygon

def generate_h3_grid(lat, # latitude of the center point 
                     lon, # longitude of the center point
                     radius, # radius in KM
                     resolution # 0-15, higher resolution means smaller hexagons
                     ):
    center_h3 = h3.latlng_to_cell(lat, lon, resolution)

    # get the hexagons within the specified radius
    k_distance = int(radius / (h3.get_hexagon_edge_length_avg(resolution, 'km')))
    hexagons = h3.grid_disk(center_h3, max(k_distance, 1))

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

df = pd.read_csv('data/vo_agam_release/all_samples_metadata_merged.csv')
print(df.head())