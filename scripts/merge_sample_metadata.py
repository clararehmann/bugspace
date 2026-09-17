import sys
import pandas as pd

def merge_data(species_calls_df, sample_metadata_df, sample_accession_df, sample_cohorts_df):
    """
    Merge species calls DataFrame with sample metadata DataFrame and accession DataFrame on 'sample_id'.
    """
    common=list(set(species_calls_df.columns).intersection(list(sample_metadata_df.columns)))
    print(common)
    merged_df=pd.merge(species_calls_df, sample_metadata_df, on=common, how='left')
    common=list(set(merged_df.columns).intersection(list(sample_accession_df.columns)))
    merged_df=pd.merge(merged_df, sample_accession_df, on='sample_id', how='left')
    common=list(set(merged_df.columns).intersection(list(sample_cohorts_df.columns)))
    merged_df=pd.merge(merged_df, sample_cohorts_df, on='sample_id', how='left')
    return merged_df

def main():
    args=sys.argv[1:]
    """
    Merge sample metadata with species calls.
    first argument: path to species calls file
    second argument: path to sample metadata file
    third argument (optional): path to output file, if not provided then will save to 'merged_sample_data.csv'
    """
    species_calls_path=args[0]
    sample_metadata_path=args[1]
    sample_accession_path=args[2]
    sample_cohorts_path=args[3]
    output_path=args[4] if len(args) > 4 else 'merged_sample_data.csv'

    species_calls_df=pd.read_csv(species_calls_path)
    sample_metadata_df=pd.read_csv(sample_metadata_path)
    sample_accession_df=pd.read_csv(sample_accession_path)
    sample_cohorts_df=pd.read_csv(sample_cohorts_path)
    merged_df=merge_data(species_calls_df, sample_metadata_df, sample_accession_df, sample_cohorts_df)
    merged_df.to_csv(output_path, index=False)

if __name__=="__main__":
    main()
