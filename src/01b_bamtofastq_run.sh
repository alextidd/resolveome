#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J bamtofastq -o log/%J_bamtofastq.out -e log/%J_bamtofastq.err 'bash src/01b_bamtofastq_run.sh'

# dirs
wd=$(pwd)
mkdir -p out/bamtofastq/

# modules
module load singularity/3.11.4

# # generate samplesheet
# cat data/resolveome/DNA/samplesheet_local.csv |
# awk -F',' 'BEGIN {OFS=","; print "mapped,index,sample_id,file_type"} NR>1 {print $5, $5 ".bai", $7, "bam"}' \
# > out/bamtofastq/samplesheet.csv

# run on DNA
(
  cd out/bamtofastq/
  nextflow run nf-core/bamtofastq \
    -profile singularity,sanger \
    --input samplesheet.csv \
    --outdir . \
    -w $wd/work/bamtofastq/ \
    -resume
)

# run on RNA