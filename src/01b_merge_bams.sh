#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 01b_merge_bams -o log/%J_01b_merge_bams.out -e log/%J_01b_merge_bams.err 'bash src/01b_merge_bams.sh'

# dirs
wd=$(pwd)

# create rna samplesheet
cat data/resolveome/samplesheet_local.csv | head -1 \
> out/merge_bams/samplesheet.csv
cat data/resolveome/samplesheet.csv |
awk -F, -v OFS="," \
  'NR>1 && $5 == "rna" && ($6 == 49901 || $6 == 50072) {
    $2 = $1"_"$5"_merged"
    print $0
  }' >> out/merge_bams/samplesheet.csv

# run
(
  cd out/merge_bams/
  nextflow run $wd/../nextflow/nf-get_bam \
    --samplesheet samplesheet.csv \
    --location local \
    --out_dir ./ \
    --merge_bams \
    -resume \
    -w $wd/work/merge_bams/ \
    -N at31@sanger.ac.uk
)
