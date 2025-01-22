#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J nf-sarek -o log/%J_nf-sarek.out -e log/%J_nf-sarek.err 'bash src/02b_nf-sarek_run.sh'

# dir
wd=$(pwd)

# run
(
  cd out/nf-sarek/
  nextflow run nf-core/sarek \
    -r 3.5.0 \
    --input samplesheet.csv \
    --outdir ./ \
    --genome GATK.GRCh38 \
    -w $wd/work/nf-sarek/ \
    --step variant_calling \
    --tools haplotypecaller \
    -profile singularity \
    -resume \
    -N at31@sanger.ac.uk
)