#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 03b_nf-resolveome_dna_run -o log/%J_03b_nf-resolveome_dna_run.out -e log/%J_03b_nf-resolveome_dna_run.err 'bash src/03b_nf-resolveome_dna_run.sh'

# dirs
wd=$(pwd)

# run
(
  cd out/nf-resolveome/dna/
  nextflow run $wd/../nextflow/nf-resolveome \
    --samplesheet samplesheet.csv \
    --bait_set $wd/data/immune_panel/Probes_merged_ok_combined_Sanger_Immune-v1_TE-91661256_hg19.bed \
    --genes TNFRSF14,CD274,LTB,TNFAIP3,TET2,DUSP2,CCR6,DNMT3A,RFTN1,CBL,RASA2,CXCR3,ACTG1,KLHL6 \
    --fasta /lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa \
    --location local \
    --bamtofastq false \
    --out_dir ./ \
    -w $wd/work/nf-resolveome/dna/ \
    -resume \
    -N at31@sanger.ac.uk \
    -with-tower
)
