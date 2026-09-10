#!/bin/bash
rm data/vo_agam_release/gambiae_samples.txt
rm data/vo_agam_release/gambiae_latlon.txt

YEARS=("2017" "2018" "2019" "2020" "2021" "2022")
for y in "${YEARS[@]}"; do
    grep 'gambiae' data/vo_agam_release/all_samples_metadata_merged.csv | grep $y | cut -d',' -f 1,14,15 | while IFS=, read sample lat lon; do
        if [ -n $lat ] && [ -n $lon ]; then
            if [[ $(echo "$lat < 40" | bc -l) == "1" ]]; then 
            echo $sample >> data/vo_agam_release/gambiae_samples.txt
            echo $lat,$lon >> data/vo_agam_release/gambiae_latlon.txt
            fi
        fi
    done
done