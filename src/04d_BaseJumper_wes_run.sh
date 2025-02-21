#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 04d_BaseJumper_wes_run -o log/%J_04d_BaseJumper_wes_run.out -e log/%J_04d_BaseJumper_wes_run.err 'bash src/04d_BaseJumper_wes_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# run bj-wgs
(
  cd out/BaseJumper/bj-wes/
  nextflow run $wd/../nextflow/external/BaseJumper/bj-wes \
    --input_csv samplesheet.csv \
    --publish_dir ./ \
    -c ~/.nextflow/config \
    -c $wd/config/bj-wes.config \
    -w $wd/work/BaseJumper/bj-wes/ \
    -c $wd/config/basejumper.config \
    --architecture "x86_64" \
    --dnascope_model_selection bioskryb129 \
    -profile singularity \
    -resume \
    -N at31@sanger.ac.uk
)