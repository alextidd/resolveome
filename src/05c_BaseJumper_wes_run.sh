#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q week -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J bj-wes -o log/%J_bj-wes.out -e log/%J_bj-wes.err 'bash src/05c_BaseJumper_wes_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# run bj-wgs
(
  cd out/BaseJumper/bj-wes
  nextflow run $wd/../nextflow/external/BaseJumper/bj-wes \
    --input_csv samplesheet.csv \
    --publish_dir ./ \
    -c ~/.nextflow/config \
    -c $wd/config/bj-wes.config \
    -w $wd/work/BaseJumper/bj-wes/ \
    --architecture "x86_64" \
    --dnascope_model_selection bioskryb129 \
    -profile singularity \
    --process.containerOptions '--bind /lustre,/nfs,/data,/software' \
    -resume
)