#!/bin/bash

# wget agam vcfs from sanger website

# gambiae
YEARS=("2017" "2018" "2019" "2020" "2021" "2022")
for y in "${YEARS[@]}"; do
    echo "Downloading gambiae vcfs for year $y"
    sbatch scripts/agam_gambiae_download.batch $y
done

# coluzzii
#YEARS=("2019" "2020" "2021" "2022")
#for y in "${YEARS[@]}"; do
#    sbatch scripts/agam_coluzzii_download.batch $y
#done

# scrap arabiensis for now
# arabiensis
#YEARS=("2017" "2018" "2019" "2020" "2021" "2022" "2023")
#for y in "${YEARS[@]}"; do
#    sbatch scripts/agam_arabiensis_download.batch $y
#done