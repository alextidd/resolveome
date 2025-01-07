#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J caveman -o log/caveman_%J.out -e log/caveman_%J.err 'bash src/02b_nf-caveman_run.sh'

# run
module load singularity
nextflow run ../nextflow/nf-caveman \
  --samplesheet out/nf-caveman/samplesheet.csv \
  --outdir out/nf-caveman/ \
  --tum_cn_default 10 \
  --norm_cn_default 2 \
  --seqType genome \
  -c ../nextflow/nf-caveman/config/sanger_hg37.config \
  -w work/caveman/ \
  -resume \
  -N at31@sanger.ac.uk