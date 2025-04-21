#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 03c_BaseJumper_expression_rna_run -o log/%J_03c_BaseJumper_expression_rna_run.out -e log/%J_03c_BaseJumper_expression_rna_run.err 'bash src/03c_BaseJumper_expression_rna_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/BaseJumper/bj-expression/sentieon_eval.lic 
export LSB_EXCLUSIVE=Y

# run bj-expression
(
  cd out/BaseJumper/bj-expression/rna/
  nextflow run $wd/../nextflow/external/BaseJumper/bj-expression \
    --input_csv samplesheet.csv \
    --publish_dir PD63118 \
    --genome GRCh38 \
    --skip_subsampling \
    --timestamp run \
    --tmp_dir $TMPDIR \
    -w $wd/work/BaseJumper/bj-expression/filter_cells/ \
    -c $wd/config/basejumper.config \
    -c $wd/config/bj-expression.config \
    -profile singularity \
    -resume
)

# bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J 03c_BaseJumper_expression_rna_symlink -o "log/%J_03c_BaseJumper_expression_rna_symlink.out" ". ~/.bashrc ; replace_symlinks out/BaseJumper/bj-expression/rna/PD63118_run/"