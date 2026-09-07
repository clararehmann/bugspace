#!/bin/bash

# wget agam vcfs from sanger website

URL="https://vo_agam_output.cog.sanger.ac.uk/"

# gambiae
YEARS=("2015" "2016" "2017" "2018" "2019" "2020" "2021" "2022" "2023")
for y in "${YEARS[@]}"; do
    echo "Downloading gambiae vcfs for year $y"
    mkdir -p data/vo_agam_release/vcfs/gambiae/$y/
    grep $y data/vo_agam_release/all_samples_metadata_merged.csv | grep "gambiae" | cut -d, -f1 | while read sample; do
        sample=$(echo $sample | tr -d '"')
        # wget -O data/vo_agam_release/vcfs/gambiae/$y/$sample.vcf.gz $URL$sample
        echo $URL$sample
    done
done

# coluzzii
YEARS=("2019" "2020" "2021" "2022")
for y in "${YEARS[@]}"; do
    echo "Downloading coluzzii vcfs for year $y"
    mkdir -p data/vo_agam_release/vcfs/coluzzii/$y/
    grep $y data/vo_agam_release/all_samples_metadata_merged.csv | grep "coluzzii" | cut -d, -f1 | while read sample; do
        sample=$(echo $sample | tr -d '"')
        # wget -O data/vo_agam_release/vcfs/coluzzii/$y/$sample.vcf.gz $URL$sample
        echo $URL$sample
    done
done

# arabiensis
YEARS=("2017" "2018" "2019" "2020" "2021" "2022" "2023")
for y in "${YEARS[@]}"; do
    echo "Downloading arabiensis vcfs for year $y"
    mkdir -p data/vo_agam_release/vcfs/arabiensis/$y/
    grep $y data/vo_agam_release/all_samples_metadata_merged.csv | grep "arabiensis" | cut -d, -f1 | while read sample; do
        sample=$(echo $sample | tr -d '"')
        # wget -O data/vo_agam_release/vcfs/arabiensis/$y/$sample.vcf.gz $URL$sample
        echo $URL$sample
    done
done