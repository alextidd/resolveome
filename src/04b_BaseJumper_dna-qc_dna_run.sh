#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J 04b_BaseJumper_dna-qc_dna_run -o log/%J_04b_BaseJumper_dna-qc_dna_run.out -e log/%J_04b_BaseJumper_dna-qc_dna_run.err 'bash src/04b_BaseJumper_dna-qc_dna_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/BaseJumper/bj-dna-qc/sentieon_eval.lic
export LSB_EXCLUSIVE=Y

# run
(
  cd out/BaseJumper/bj-dna-qc/dna/ ;
  
  nextflow run $wd/../nextflow/external/BaseJumper/bj-dna-qc \
    --input_csv samplesheet.csv \
    --publish_dir PD63118 \
    --timestamp run \
    --skip_ginkgo false \
    -c $wd/config/basejumper.config \
    -w $wd/work/BaseJumper/bj-dna-qc/dna/ \
    -profile singularity \
    -N at31@sanger.ac.uk \
    -resume 
)

# bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J 04b_BaseJumper_dna-qc_dna_symlink -o "log/%J_04b_BaseJumper_dna-qc_dna_symlink.out" ". ~/.bashrc ; replace_symlinks out/BaseJumper/bj-dna-qc/dna/PD63118_run/"