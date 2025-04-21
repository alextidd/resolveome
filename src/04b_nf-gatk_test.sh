#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 04b_nf-gatk_test -o log/%J_04b_nf-gatk_test.out -e log/%J_04b_nf-gatk_test.err 'bash src/04b_nf-gatk_test.sh'

# dir
wd=$(pwd)
export TMPDIR=/lustre/scratch125/casm/team268im/at31/tmp/bash/

# modules
module load singularity/3.11.4

# run
(
  cd out/test/nf-gatk/
  nextflow run $wd/../nextflow/nf-gatk/ \
    --samplesheet samplesheet.csv \
    --outdir ./ \
    --location local \
    --fasta /lustre/scratch125/casm/team268im/fa8/117/PTA_49686/PTATO/resources/hg38/Homo_sapiens.GRCh38.dna_sm.toplevel.fa \
    --realign \
    -w $wd/work/nf-gatk_test/ \
    -N at31@sanger.ac.uk \
    -resume
)