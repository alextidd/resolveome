#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 02c_nf-resolveome_dnahyb_run -o log/%J_02c_nf-resolveome_dnahyb_run.out -e log/%J_02c_nf-resolveome_dnahyb_run.err 'bash src/02c_nf-resolveome_dnahyb_run.sh'

# dirs
wd=$(pwd)

# run
(
  cd out/nf-resolveome/dnahyb/
  nextflow run $wd/../nextflow/nf-resolveome \
    --samplesheet samplesheet.csv \
    --bait_set_hyb $wd/data/immune_panel/Probes_merged_ok_combined_Sanger_Immune-v1_TE-91661256_hg19.bed \
    --bait_set_vdj $wd/out/vdj_coverage/regions/ig_tcr_genes.bed \
    --genes $wd/data/driver_genes/PD63118_driver_genes.txt \
    --fasta /lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCh37d5/genome.fa \
    --location local \
    --out_dir ./ \
    -w $wd/work/nf-resolveome/dnahyb/ \
    -resume \
    -N at31@sanger.ac.uk \
    -with-tower
)

