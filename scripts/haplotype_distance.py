import numpy as np, pandas as pd
import os, sys
import dask
import dask.array as da
from dask.diagnostics.progress import ProgressBar
# silence some warnings
dask.config.set(**{'array.slicing.split_large_chunks': False})
import allel
import plotly
import plotly.express as px
from gambiae_grid import *
import h3
from shapely.geometry import Polygon, Point, mapping, LineString
import shapely
from pyproj import Geod
import itertools, argparse

import malariagen_data
ag3 = malariagen_data.Ag3()

def get_query(years, sample_sets):
    # generate query string for {year} samples
    years='|'.join([f'year == {year}' for year in years])
    sample_query="taxon=='gambiae' & " + years
    print(sample_query)
    df_samples = ag3.sample_metadata(sample_sets=sample_sets, sample_query=sample_query)
    return df_samples, sample_query

def shared_haplotype(gt1, gt2):
    # where genotypes are the same
    gt_match = np.where(gt1 == gt2)[0]
    # position distance between matches
    gt_diffs = np.diff(gt_match)
    # split the array where the difference is not one
    splits = np.split(gt_match, np.where(gt_diffs != 1)[0] + 1)
    return splits

def haplotype_between(ind1, ind2, gt, pos):
    # compare 4 ways between all chromosomes
    maximum = 0
    for chrom_comb in [(0,0), (0,1), (1,0), (1,1)]:
        print(chrom_comb)
        splits = shared_haplotype(gt[:,ind1,chrom_comb[0]], gt[:,ind2,chrom_comb[1]])
        lengths = [pos[i[-1]] - pos[i[0]] for i in splits]
        length = max(lengths)
        if length > maximum:
            maximum = length
    return maximum

def assign_grid_cells(sample_metadata, resolution=4):
    ## assign each sample to h3 grid cell
    sample_metadata['h3_index'] = [h3.latlng_to_cell(sample_metadata.iloc[i].latitude,
                                                     sample_metadata.iloc[i].longitude,
                                                     resolution) for i in range(len(sample_metadata))]
    #lat_c, lon_c, max_d = get_sampling_radius(sample_metadata)
    #resolution=4
    #grid = generate_h3_grid(lat_c, lon_c, max_d, resolution=resolution)
    # assign cells to samples
    #sample_metadata['h3_index'] = None
    #sample_metadata['point'] = [Point(sample_metadata.iloc[i].longitude, sample_metadata.iloc[i].latitude) for i in range(len(sample_metadata))]
    #for cell in range(len(grid)):
    #    mask = [shapely.contains(grid.iloc[cell].geometry, sample_metadata.iloc[i].point) for i in range(len(sample_metadata))]
    #    sample_metadata.loc[mask, 'h3_index'] = grid.iloc[cell].h3_index
    return sample_metadata


parser=argparse.ArgumentParser()
parser.add_argument('--years', nargs='+', type=int)
parser.add_argument('--region', type=str, default='3L')
parser.add_argument('--resolution', type=int, default=3, help='h3 grid resolution')
parser.add_argument('--out', type=str, help='outpath prefix (appended with _haplotype_distances.csv)')
args = parser.parse_args()

sample_sets = [f'3.{i}' for i in range(17)]
sample_metadata, sample_query = get_query(args.years, sample_sets=sample_sets)

haps = ag3.haplotypes(
    region=args.region,
    sample_query=sample_query
)

# reorder sample metadata to match zarr sample indexing
print('reindexing sample data...')
sample_metadata.set_index('sample_id', inplace=True)
samples = haps['sample_id'].values
sample_metadata = sample_metadata.reindex(samples)
sample_metadata = sample_metadata.reset_index(names=['sample_id'])
print('assigning grid cells...')
sample_metadata = assign_grid_cells(sample_metadata, resolution=args.resolution)
# extract genotypes and positions
genotypes = haps['call_genotype'].values
positions = haps['variant_position'].values
geod = Geod(ellps="WGS84")
# combinations of samples
idxs = np.arange(len(sample_metadata))
comb = np.array(list(itertools.combinations(idxs, 2)))

#lengths = np.empty(len(comb))
#hex_dists = np.empty(len(comb))
#km_dists = np.empty(len(comb))
#sample_1 = np.empty(len(comb), dtype='object')
#sample_2 = np.empty(len(comb), dtype='object')
#month_1 = np.empty(len(comb))
#month_2 = np.empty(len(comb))
#lon_1 = np.empty(len(comb))
#lat_1 = np.empty(len(comb))
#lon_2 = np.empty(len(comb))
#lat_2 = np.empty(len(comb))
#print(genotypes.sizes)
print('getting haplotype distances...')
with open(f'{args.out}_haplotype_distance.csv', 'a') as f:
    f.write('length,hex_distance,km_distance,sample_1,sample_2,month_1,month_2,lon_1,lat_1,lon_2,lat_2\n')
    for c in range(len(comb)):
        try:
            print(c)
            i = comb[c,0]
            j = comb[c,1]
            length = haplotype_between(i, j, genotypes, positions)
            hex_dist = h3.grid_distance(sample_metadata.iloc[i].h3_index, sample_metadata.iloc[j].h3_index)
            km_dist = geod.inv(sample_metadata.iloc[i].longitude,
                                    sample_metadata.iloc[i].latitude,
                                    sample_metadata.iloc[j].longitude,
                                    sample_metadata.iloc[j].latitude)[2]*0.001
            sample_1 = sample_metadata.iloc[i].sample_id
            sample_2 = sample_metadata.iloc[j].sample_id
            month_1 = sample_metadata.iloc[i].month
            month_2 = sample_metadata.iloc[j].month
            lon_1 = sample_metadata.iloc[i].longitude
            lat_1 = sample_metadata.iloc[i].latitude
            lon_2 = sample_metadata.iloc[j].longitude
            lat_2 = sample_metadata.iloc[j].latitude
            f.write(f'{length},{hex_dist},{km_dist},{sample_1},{sample_2},{month_1},{month_2},{lon_1},{lat_1},{lon_2},{lat_2}\n')
        except:
            pass
#df = pd.DataFrame({'length':lengths,
#                   'hex_distance':hex_dists,
#                   'km_distance':km_dists,
#                   'sample_1':sample_1,
#                   'sample_2':sample_2,
#                   'month_1':month_1,
#                   'month_2':month_2,
#                   'lon_1':lon_1,
#                   'lat_1':lat_1,
#                   'lon_2':lon_2,
#                   'lat_2':lat_2})
#df.to_csv(f'{args.out}_haplotype_distances.csv')

#print(sample_metadata.iloc[i1].longitude, sample_metadata.iloc[i1].latitude)
#print(sample_metadata.iloc[i2].longitude, sample_metadata.iloc[i2].latitude)

