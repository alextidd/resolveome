#!/bin/bash

wd=$(pwd)
(
  cd out/BaseJumper/
  grep -e "biosampleName\|49686" samplesheet.csv > samplesheet_49686.csv
  nextflow run $wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling \
    --input_csv samplesheet_49686.csv \
    --publish_dir ./ \
    -c ~/.nextflow/config \
    -c $wd/config/bj-somatic-variantcalling.config \
    --is_bam \
    -w $wd/
)
