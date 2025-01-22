#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J get_bams -o log/get_bams_%J.out -e log/get_bams_%J.err 'bash src/01_get_bams.sh'

# module
module load IRODS/1.0

# run nf-get_bam for DNA
nextflow run ../nextflow/nf-get_bam \
  --samplesheet data/resolveome/DNA/samplesheet_irods.csv \
  --out_dir data/resolveome/DNA/ \
  --cram_to_bam \
  -resume

# run nf-get_bam for RNA
nextflow run ../nextflow/nf-get_bam \
  --samplesheet data/resolveome/RNA/samplesheet_irods.csv \
  --out_dir data/resolveome/RNA/ \
  --cram_to_bam \
  -resume