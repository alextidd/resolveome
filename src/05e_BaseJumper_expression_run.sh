#!/bin/bash

# dirs
wd=$(pwd)

# modules
module load singularity

# run bj-expression
(
  cd out/BaseJumper/bj-expression/
  nextflow run $wd/../nextflow/external/BaseJumper/bj-expression \
    --input_csv samplesheet.csv \
    --publish_dir ./ \
    -c ~/.nextflow/config \
    -c $wd/config/bj-expression.config \
    -w $wd/work/BaseJumper/bj-expression/ \
    --architecture "x86_64" \
    --dnascope_model_selection bioskryb129 \
    -profile singularity \
    --process.containerOptions '--bind /lustre,/nfs,/data,/software' \
    -resume
)