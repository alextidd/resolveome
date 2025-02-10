#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M20000 -R 'span[hosts=1] select[mem>20000] rusage[mem=20000]' -J bj-wgs -o log/%J_bj-wgs.out -e log/%J_bj-wgs.err 'bash src/05b_BaseJumper_wgs_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# run bj-wgs
(
  cd out/BaseJumper/bj-wgs
  nextflow run $wd/../nextflow/external/BaseJumper/bj-wgs \
    --input_csv samplesheet.csv \
    --publish_dir ./ \
    --sentieon_license /lustre/scratch125/casm/team268im/at31/nextflow/external/BaseJumper/bj-wgs/sentieon_eval.lic \
    -c ~/.nextflow/config \
    -c $wd/config/bj-wgs.config \
    -w $wd/work/BaseJumper/bj-wgs/ \
    -profile singularity \
    --architecture "x86_64" \
    --dnascope_model_selection bioskryb129 \
    --process.containerOptions '--bind /lustre,/nfs,/data,/software' \
    -resume
)