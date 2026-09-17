#!/bin/bash
rm data/vo_agam_release/gambiae_samples.txt
rm data/vo_agam_release/gambiae_latlon.txt

YEARS=("2017" "2018" "2019" "2020" "2021" "2022")
for y in "${YEARS[@]}"; do
    grep 'gambiae' data/vo_agam_release/all_samples_metadata_merged.csv | grep $y | cut -d',' -f 1,12,14,15 | while IFS=, read sample year lat lon; do
        if [ -n $lat ] && [ -n $lon ]; then
            if [[ $(echo "$lat < 40" | bc -l) == "1" ]]; then 
                if [[ "$y" == "$year" ]]; then
                    echo $sample >> data/vo_agam_release/gambiae_samples.txt
                    echo $lat,$lon >> data/vo_agam_release/gambiae_latlon.txt
                fi
            fi
        fi
    done
done

# make test with just two samples
samples=$(pixi exec -c bioconda bcftools -l /hps/nobackup/jlees/rehmann/vo_agam_release/vcfs/gambiae/gambiae_merged_1451_1500.vcf.gz)
for s in "${samples[@]}"; do
    grep $s data/vo_agam_release/all_samples_metadata_merged.csv | cut -d',' -f 1,14,15 | while IFS=, read sample lat lon; do
        echo $lat,$lon >> data/test_latlon.txt
    done
done