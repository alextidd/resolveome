#!/bin/bash

# modules
module load bwa-0.7.17
module load samtools-1.19.2/python-3.11.6

# dirs
wd=$(pwd)
fastq_dir=$wd/out/bamtofastq/reads/

# generate samplesheet
echo "donor_id,id,fastq_1,fastq2" > out/nf-gatk/samplesheet.csv
while read -r donor_id id ; do
  echo "$donor_id,$id,$fastq_dir/${id}_1.merged.fastq.gz,${id}_2.merged.fastq.gz" \
  >> out/bwa_mem/samplesheet.csv
done < <(cat data/resolveome/DNA/samplesheet_local.csv | sed 1d | cut -d, -f6,7)

# remap
bwa mem \
  /lustre/scratch125/casm/team268im/fa8/117/PTA_49686/PTATO/resources/hg38/Homo_sapiens.GRCh38.dna_sm.toplevel.fa \
  out/bamtofastq/reads/49882_plex1_1.merged.fastq.gz \
  out/bamtofastq/reads/49882_plex1_2.merged.fastq.gz |
samtools sort -@ 8 -o realigned.bam
