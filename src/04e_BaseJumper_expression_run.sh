#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J 04e_BaseJumper_expression_run -o log/%J_04e_BaseJumper_expression_run.out -e log/%J_04e_BaseJumper_expression_run.err 'bash src/04e_BaseJumper_expression_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/BaseJumper/bj-expression/sentieon_eval.lic 

# run bj-expression
(
  cd out/BaseJumper/bj-expression/
  nextflow run $wd/../nextflow/external/BaseJumper/bj-expression \
    --input_csv samplesheet.csv \
    --publish_dir PD63118 \
    --genome GRCh38 \
    --skip_subsampling \
    --tmp_dir $TMPDIR \
    -w $wd/work/BaseJumper/bj-expression/ \
    -c $wd/config/bj-expression.config \
    -profile singularity \
    -resume
)