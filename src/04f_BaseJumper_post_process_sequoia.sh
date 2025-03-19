#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M15000 -R 'span[hosts=1] select[mem>15000] rusage[mem=15000]' -J 04f_BaseJumper_post_process_sequoia -o log/%J_04f_BaseJumper_post_process_sequoia.out -e log/%J_04f_BaseJumper_post_process_sequoia.err 'bash src/04f_BaseJumper_post_process_sequoia.sh'

# modules
module load ISG/rocker/rver/4.4.0

# dirs
seq_dir=out/BaseJumper/bj-somatic-variantcalling/dna/PD63118_250306_232001/SOMATIC_VARIANT_WORKFLOW_Heuristic_Filter_SEQUOIA/
out_dir=out/BaseJumper/bj-somatic-variantcalling/dna/PD63118_250306_232001/post_process_sequoia/
mkdir $out_dir

# run sequoia post-processing
Rscript bin/post_process_sequoia.R $seq_dir $out_dir