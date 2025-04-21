#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 03e_BaseJumper_somatic-variantcalling_dnahyb_run -o log/%J_03e_BaseJumper_somatic-variantcalling_dnahyb_run.out -e log/%J_03e_BaseJumper_somatic-variantcalling_dnahyb_run.err 'bash src/03e_BaseJumper_somatic-variantcalling_dnahyb_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling/sentieon_eval.lic
export LSB_EXCLUSIVE=Y

# run
(
  cd out/BaseJumper/bj-somatic-variantcalling/dnahyb/
  
  nextflow run $wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling \
    --input_csv samplesheet.csv \
    --publish_dir PD63118 \
    --timestamp run \
    --variant_workflow_type somatic_heuristic_filter \
    --dnascope_model_selection bioskryb129 \
    --skip_variant_annotation false \
    --skip_sigprofile false \
    -c $wd/config/bj-somatic-variantcalling_dnahyb.config \
    -c $wd/config/bj-somatic-variantcalling.config \
    -c $wd/config/basejumper.config \
    -w $wd/work/BaseJumper/bj-somatic-variantcalling/dnahyb/ \
    -profile singularity \
    --architecture "x86_64" \
    -resume \
    -with-tower \
    -N at31@sanger.ac.uk
)

# bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J 03e_BaseJumper_somatic-variantcalling_dnahyb_symlink -o "log/%J_03e_BaseJumper_somatic-variantcalling_dnahyb_symlink.out" "replace_symlinks out/BaseJumper/bj-somatic-variantcalling/dnahyb/PD63118_run/"