#!/bin/bash
# script to streamline downloading and combining 

# download agam releases from google cloud
releases=("3" "3.1" "3.10" "3.11" "3.12" "3.13" "3.14" "3.15" "3.16" "3.17" "3.2" "3.3" "3.4" "3.5" "3.6" "3.7" "3.8" "3.9")
for v in "${releases[@]}"; do
    echo "Downloading agam release $v"
    mkdir -p data/vo_agam_release/v$v/
    # download metadata
    gcloud storage cp -r gs://vo_agam_release_master_us_central1/v$v/metadata data/vo_agam_release/v$v/
    # merge cohort metadata files
    awk 'NR == 1 || FNR > 1' data/vo_agam_release/v$v/metadata/general/*/samples.meta.csv > data/vo_agam_release/v$v/metadata/all_samples_metadata.csv
    # merge species metadata files
    awk 'NR == 1 || FNR > 1' data/vo_agam_release/v$v/metadata/species_calls_aim_20220528/*/samples.species_aim.csv > data/vo_agam_release/v$v/metadata/all_samples_species_calls.csv
    # merge the two files
    Rscript scripts/merge_sample_metadata.R data/vo_agam_release/v$v/metadata/all_samples_species_calls.csv data/vo_agam_release/v$v/metadata/all_samples_metadata.csv data/vo_agam_release/v$v/metadata/all_samples_merged.csv
done
# combine all releases into one file
awk 'NR == 1 || FNR > 1' data/vo_agam_release/v*/metadata/all_samples_merged.csv > data/vo_agam_release/all_samples_metadata_merged.csv