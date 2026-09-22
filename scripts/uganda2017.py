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
from ag3_popgen import *
import malariagen_data
ag3 = malariagen_data.Ag3()

"""
Investigating how population structure and diversity changes over the course of a year
Does seasonality influence this? And if so, it should be incorporated into our spatial model
"""

def get_month_query(month, sample_sets):
    # generate query string for 2017 samples from Uganda in that month
    df_samples = ag3.sample_metadata(sample_sets=sample_sets, sample_query=f"country=='Uganda' & year==2017 & taxon=='gambiae' & month=={month}")
    samples = df_samples['sample_id'].tolist()
    sample_query = '|'.join([f"sample_id=='{s}'" for s in samples])
    return df_samples, sample_query
outpath='out/Uganda2017/'
os.makedirs(outpath, exist_ok=True)

months=[3,4,5,6,7,8,9,11]
sample_sets = [f'3.{i}' for i in range(17)]

# get metadata
month=sys.argv[1]
sample_metadata, sample_query = get_month_query(month, sample_sets=sample_sets)
print(f'Month {month}: {len(sample_metadata)} samples')
sample_query=f"country=='Uganda' & year==2017 & taxon=='gambiae' & month=={month}"
# run PCA on samples
print('Running PCA...')
pca_plot = run_plot_pca(sample_query)
plotly.offline.plot(pca_plot, filename=f'{outpath}{month}_PCA.html')
# get popgen stats
## decide if these should be all saved to one dataframe over the course of the year...
print('Calculating popgen stats...')
cohort_sizes = [len(sample_metadata[sample_metadata.cohort_admin2_year==c]) for c in np.unique(sample_metadata.cohort_admin2_year)]
stats_df = run_popgen_stats(sample_query, min(cohort_sizes))
stats_df.to_csv(f'{outpath}{month}_stats.csv')
print('Calculating between-cohort FST')
#genotypes = get_zarr_genotypes(sample_query)
# get between-location FST
#sample_metadata, fst_values = pairwise_fst(sample_metadata, genotypes)
#print(np.mean(fst_values.values()), np.median(fst_values.values()), min(fst_values.values()), max(fst_values.values()))
ccombs, fsts, stderrs = between_cohort_fst(sample_metadata)
ch1 = [c[0] for c in ccombs]
ch2 = [c[1] for c in ccombs]
fst_df = pd.DataFrame({'cohort_1':ch1, 'cohort_2':ch2, 'average_fst':fsts, 'standard_error':stderrs})
fst_df.to_csv(f'{outpath}{month}_fst.csv')