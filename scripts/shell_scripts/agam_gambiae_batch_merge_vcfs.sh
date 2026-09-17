#!/bin/bash
VCFS=$(find /hps/nobackup/jlees/rehmann/vo_agam_release/vcfs/gambiae/20* -name "*.vcf.gz" -size +10)

STEP=50
END=$STEP
LENGTH=$(echo "$VCFS" | wc -l)
#LENGTH=10
for START in `seq 1 $STEP $LENGTH`; do
    SUBSET=$(tail -n +$START <<< "$VCFS" | head -n $STEP)
    SUBSET=$(echo "$SUBSET" | tr '\n' ' ')
    END=$(($START + $STEP - 1))
    if [ $END -gt $LENGTH ]; then
        END=$LENGTH
    fi
    sbatch scripts/agam_gambiae_batch_merge_vcfs.batch "$SUBSET" $START $END
done
