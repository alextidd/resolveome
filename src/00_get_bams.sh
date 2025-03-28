#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 00_get_bams -o log/%J_00_get_bams.out -e log/%J_00_get_bams.err 'bash src/00_get_bams.sh'

# module
module load IRODS/1.0

# dirs
wd=$(pwd)

# # run
# (
#   cd data/resolveome/
#   nextflow run $wd/../nextflow/nf-get_bam \
#     --samplesheet samplesheet_irods.csv \
#     --location irods \
#     --out_dir ./ \
#     --cram_to_bam \
#     --merge_bams \
#     -resume \
#     -w $wd/work/get_bams/ \
#     -N at31@sanger.ac.uk
# )

# run
(
  cd data/resolveome/
  cat samplesheet_irods.csv | grep "cell_id\|plate10_" \
  > samplesheet_irods.tmp.csv
  nextflow run $wd/../nextflow/nf-get_bam \
    --samplesheet samplesheet_irods.tmp.csv \
    --location irods \
    --out_dir ./ \
    --cram_to_bam \
    -resume \
    -w $wd/work/get_bams/ \
    -N at31@sanger.ac.uk
)