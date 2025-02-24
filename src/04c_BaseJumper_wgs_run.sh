#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 04c_BaseJumper_wgs_run -o log/%J_04c_BaseJumper_wgs_run.out -e log/%J_04c_BaseJumper_wgs_run.err 'bash src/04c_BaseJumper_wgs_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# run bj-wgs
(
  cd out/BaseJumper/bj-wgs/
  nextflow run $wd/../nextflow/external/BaseJumper/bj-wgs \
    --input_csv samplesheet.csv \
    --publish_dir ./ \
    -c ~/.nextflow/config \
    -c $wd/config/bj-wgs.config \
    -c $wd/config/basejumper.config \
    -w $wd/work/BaseJumper/bj-wgs/ \
    --architecture "x86_64" \
    --dnascope_model_selection bioskryb129 \
    -profile singularity \
    #-resume \
    -N at31@sanger.ac.uk
)


# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling/sentieon_eval.lic 

# run bj-wgs + bj-somatic-variantcalling
(
  cd out/BaseJumper/bj-somatic-variantcalling/
  
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