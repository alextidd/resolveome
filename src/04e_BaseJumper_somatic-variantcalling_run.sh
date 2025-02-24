#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 04e_BaseJumper_somatic-variantcalling_run -o log/%J_04e_BaseJumper_somatic-variantcalling_run.out -e log/%J_04e_BaseJumper_somatic-variantcalling_run.err 'bash src/04e_BaseJumper_somatic-variantcalling_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling/sentieon_eval.lic 

# run bj-somatic-variantcalling
(
  cd out/BaseJumper/bj-wes-somatic-variantcalling/
  
  # run
  # Viren: switch --variant_workflow_type to somatic_heuristic_filter - doesn't 
  # require a matched normal (for bj-somatic-variantcalling)
  nextflow run $wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling \
    --input_csv samplesheet.csv \
    --publish_dir ./ \
    --variant_workflow_type somatic_heuristic_filter \
    --is_bam \
    -c ~/.nextflow/config \
    -c $wd/config/bj-somatic-variantcalling.config \
    -c $wd/config/basejumper.config \
    -w $wd/work/BaseJumper/bj-somatic-variantcalling/ \
    -profile singularity \
    --architecture "x86_64" \
    -resume
)