#!/bin/bash

# modules
module load HGI/softpack/users/eh19/test-subread/1

# run featurecounts
featureCounts \
  -a ../reference/gencode/GRCh38/gencode.v40.annotation.gtf.gz  \
  -t exon \
  -g gene_id \
  --primary \
  -p \
  -C \
  -o out/featureCounts/PD63118_49901.featureCounts.txt \
  data/resolveome/RNA/PD63118_49901/plex*/bam/plex*.bam