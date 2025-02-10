#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J bj-somatic-variantcalling -o log/%J_bj-somatic-variantcalling.out -e log/%J_bj-somatic-variantcalling.err 'bash src/05d_BaseJumper_somatic-variantcalling_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling/Bioskryb_Inc_c8_eval.lic

# run bj-somatic-variantcalling
(
  cd out/BaseJumper/bj-somatic-variantcalling
  nextflow run $wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling \
    --input_csv samplesheet.csv \
    --publish_dir ./ \
    -c ~/.nextflow/config \
    -c $wd/config/bj-somatic-variantcalling.config \
    --is_bam \
    -w $wd/work/BaseJumper/ \
    -profile singularity \
    --architecture "x86_64" \
    --process.containerOptions '--bind /lustre,/nfs,/data,/software' \
    -resume
)
