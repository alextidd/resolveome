#!/bin/bash

# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling/Bioskryb_Inc_c8_eval.lic

# test
(
  cd out/BaseJumper/
  grep -e "biosampleName\|49900" samplesheet.csv > samplesheet_49900.csv
  nextflow run $wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling \
    --input_csv samplesheet_49900.csv \
    --publish_dir ./ \
    -c ~/.nextflow/config \
    -c $wd/config/bj-somatic-variantcalling.config \
    --is_bam \
    -w $wd/work/BaseJumper/ \
    -profile singularity \
    --architecture "x86_64" \
    --process.containerOptions '--bind /lustre,/nfs,/data,/software' \
    -resume
)

# run
(
  cd out/BaseJumper/
  grep -e "biosampleName\|49686" samplesheet.csv > samplesheet_49686.csv
  nextflow run $wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling \
    --input_csv samplesheet_49686.csv \
    --publish_dir ./ \
    -c ~/.nextflow/config \
    -c $wd/config/bj-somatic-variantcalling.config \
    --is_bam \
    -w $wd/work/BaseJumper/
)
