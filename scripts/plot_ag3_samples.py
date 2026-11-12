import malariagen_data
import plotly.express as px
import pandas as pd

ag3=malariagen_data.Ag3()
sample_sets = [f'3.{i}' for i in range(17)]
sample_metadata = ag3.sample_metadata(sample_sets=sample_sets, sample_query="taxon=='gambiae'")

# exclude UK samples
sample_metadata = sample_metadata[sample_metadata['latitude'] < 26]
sample_counts = sample_metadata.groupby(['latitude','longitude','year']).size().reset_index(name='count')
print(sample_counts.head())
fig = px.scatter_map(sample_counts,
                     lat='latitude',
                     lon='longitude',
                     size='count',
                     color='year'
                     )
fig.write_html('out/ag3_gambiae_samples.html', include_plotlyjs=True, full_html=True)