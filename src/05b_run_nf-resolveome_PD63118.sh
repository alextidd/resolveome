#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J resolveome_PD63118 -o log/resolveome_PD63118_%J.out -e log/resolveome_PD63118_%J.err 'bash src/05b_run_nf-resolveome_PD63118.sh'

wd=$(pwd)
(
  cd out/nf-resolveome/PD63118
  nextflow run $wd/../nextflow/nf-resolveome \
    --samplesheet samplesheet.tsv \
    --bait_set $wd/data/immune_panel/Probes_merged_ok_combined_Sanger_Immune-v1_TE-91661256_hg19.bed \
    --genes TNFRSF14,CD274 \
    --fasta /lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa \
    --location irods \
    --out_dir ./ \
    -w $wd/work/resolveome/ \
    -resume \
    -N at31@sanger.ac.uk \
    -with-tower
)