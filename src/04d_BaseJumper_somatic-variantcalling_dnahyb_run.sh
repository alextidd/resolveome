#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 04d_BaseJumper_somatic-variantcalling_dnahyb_run -o log/%J_04d_BaseJumper_somatic-variantcalling_dnahyb_run.out -e log/%J_04d_BaseJumper_somatic-variantcalling_dnahyb_run.err 'bash src/04d_BaseJumper_somatic-variantcalling_dnahyb_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling/sentieon_eval.lic 
export LSB_EXCLUSIVE=Y

# run
(
  cd out/BaseJumper/bj-somatic-variantcalling/all_cells/dnahyb/
  
  nextflow run $wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling \
    --input_csv samplesheet.csv \
    --publish_dir PD63118 \
    --variant_workflow_type somatic_heuristic_filter \
    --dnascope_model_selection bioskryb129 \
    --mode exome \
    --exome_panel TWIST_immune \
    -c $wd/config/bj-somatic-variantcalling.config \
    -c $wd/config/basejumper.config \
    -w $wd/work/BaseJumper/bj-somatic-variantcalling/all_cells/dnahyb/ \
    -profile singularity \
    --architecture "x86_64" \
    -resume \
    -N at31@sanger.ac.uk
)