import numpy as np, pandas as pd
import os
import dask
import dask.array as da
from dask.diagnostics.progress import ProgressBar
# silence some warnings
dask.config.set(**{'array.slicing.split_large_chunks': False})
import allel
import plotly.express as px
import plotly
from ag3_popgen import *
import malariagen_data
ag3 = malariagen_data.Ag3()



def get_year_query(year, sample_sets):
    # generate query string for {year} samples from DRC
    sample_query=f"country=='The Democratic Republic of Congo' & year=={year} & taxon=='gambiae'"
    df_samples = ag3.sample_metadata(sample_sets=sample_sets, sample_query=sample_query)
    return df_samples, sample_query

outpath='out/DRC/
os.makedirs(outpath, exist_ok=True)

years=[2016,2017,2018]
sample_sets = [f'3.{i}' for i in range(17)]
# get metadata
year=sys.argv[1]
sample_metadata, sample_query = get_year_query(year, sample_sets=sample_sets)
print(f'Year {year}: {len(sample_metadata)} samples')
# run PCA on samples
print('Running PCA...')
pca_plot = run_plot_pca(sample_query)
plotly.offline.plot(pca_plot, filename=f'{outpath}{year}_PCA.png')
# get popgen stats
print('Calculating popgen stats...')
cohort_sizes = [len(sample_metadata[sample_metadata.cohort_admin2_year==c]) for c in np.unique(sample_metadata.cohort_admin2_year)]
stats_df = run_popgen_stats(sample_query, min(cohort_sizes))
stats_df.to_csv(f'{outpath}{year}_stats.csv')
# get between-location FST
ccombs, fsts, stderrs = between_cohort_fst(sample_metadata)
ch1 = [c[0] for c in ccombs]
ch2 = [c[1] for c in ccombs]
fst_df = pd.DataFrame({'cohort_1':ch1, 'cohort_2':ch2, 'average_fst':fsts, 'standard_error':stderrs})
fst_df.to_csv(f'{outpath}{year}_fst.csv')