#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J get_bams -o log/get_bams_%J.out -e log/get_bams_%J.err 'bash src/04a_BaseJumper_get_bams.sh'

# run nf-get_bam
nextflow run ../nextflow/nf-get_bam \
  --samplesheet out/nf-resolveome/PD63118/samplesheet.tsv \
  --out_dir out/BaseJumper/ \
  --cram_to_bam \
  -resume