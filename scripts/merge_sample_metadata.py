import sys
import pandas as pd

def read_species_calls(species_calls_path):
    """
    Read species calls file and return a pandas DataFrame.
    """
    species_calls_df=pd.read_csv(species_calls_path)
    return species_calls_df

def read_sample_metadata(sample_metadata_path):
    """
    Read sample metadata file and return a pandas DataFrame.
    """
    sample_metadata_df=pd.read_csv(sample_metadata_path)
    return sample_metadata_df

def merge_data(species_calls_df, sample_metadata_df):
    """
    Merge species calls DataFrame with sample metadata DataFrame on 'sample_id'.
    """
    merged_df=pd.merge(species_calls_df, sample_metadata_df, on='sample_id', how='left')
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
    output_path=args[2] if len(args) > 2 else 'merged_sample_data.csv'

    species_calls_df=read_species_calls(species_calls_path)
    sample_metadata_df=read_sample_metadata(sample_metadata_path)
    merged_df=merge_data(species_calls_df, sample_metadata_df)
    merged_df.to_csv(output_path, index=False)

if __name__=="__main__":
    main()
