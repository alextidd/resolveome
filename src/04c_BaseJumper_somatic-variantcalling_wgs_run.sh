#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J 04c_BaseJumper_somatic-variantcalling_wgs_run -o log/%J_04c_BaseJumper_somatic-variantcalling_wgs_run.out -e log/%J_04c_BaseJumper_somatic-variantcalling_wgs_run.err 'bash src/04c_BaseJumper_somatic-variantcalling_wgs_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# # sentieon license
# export SENTIEON_LICENSE=$wd/config/Bioskryb_Inc_c14_eval.lic

# run
(
  cd out/BaseJumper/bj-somatic-variantcalling/wgs/
  
  nextflow run $wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling \
    --input_csv samplesheet.csv \
    --publish_dir PD63118 \
    --variant_workflow_type somatic_heuristic_filter \
    -c $wd/config/bj-somatic-variantcalling.config \
    -c $wd/config/basejumper.config \
    -w $wd/work/BaseJumper/bj-somatic-variantcalling/wgs/ \
    -profile singularity \
    --architecture "x86_64" \
    --dnascope_model_selection bioskryb129 \
    -resume \
    -N at31@sanger.ac.uk
)