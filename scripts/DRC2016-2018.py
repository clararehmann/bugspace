import numpy as np, pandas as pd
import os
import dask
import dask.array as da
from dask.diagnostics.progress import ProgressBar
# silence some warnings
dask.config.set(**{'array.slicing.split_large_chunks': False})
import allel
import plotly.express as px
from ag3_popgen import *
import malariagen_data
ag3 = malariagen_data.Ag3()



def get_year_query(year, sample_sets):
    # generate query string for 2017 samples from Uganda in that month
    df_samples = ag3.sample_metadata(sample_sets=sample_sets, sample_query=f"country=='The Democratic Republic of Congo' & year=={year} & taxon=='gambiae'")
    samples = df['sample_id'].tolist()
    sample_query = '|'.join([f"sample=='{s}'" for s in samples])
    return sample_query

outpath='out/DRC/
os.makedirs(outpath, exist_ok=True)

years=[2016,2017,2018]
sample_sets = [f'3.{i}' for i in range(17)]
# get metadata
year=years[0]
sample_query = get_year_query(year, sample_sets=sample_sets)
sample_metadata = ag3.sample_metadata(sample_query=sample_query)
print(f'Year {year}: {len(sample_metadata)} samples')
# run PCA on samples
pca_plot = run_plot_pca(sample_query)
pca_plot.write_image(f'{outpath}{year}_PCA.png')
# get popgen stats
## decide if these should be all saved to one dataframe over the course of the year...
stats_df = run_popgen_stats(sample_query)
print(stats_df) 
genotypes = get_zarr_genotypes(sample_query)
# get between-location FST
sample_metadata, fst_values = pairwise_fst(sample_metadata, genotypes)
print(f'mean Fst: {np.mean(fst_values.values())}\nmedian Fst: {np.median(fst_values.values())}\nminimum Fst: {min(fst_values.values())}\nmaximum Fst: {max(fst_values.values()))}'