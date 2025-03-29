#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 04i_BaseJumper_expression_filter_cells_run -o log/%J_04i_BaseJumper_expression_filter_cells_run.out -e log/%J_04i_BaseJumper_expression_filter_cells_run.err 'bash src/04i_BaseJumper_expression_filter_cells_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/BaseJumper/bj-expression/sentieon_eval.lic 
export LSB_EXCLUSIVE=Y

# run bj-expression
(
  cd out/BaseJumper/bj-expression/filter_cells/
  nextflow run $wd/../nextflow/external/BaseJumper/bj-expression \
    --input_csv samplesheet.csv \
    --publish_dir PD63118 \
    --genome GRCh38 \
    --skip_subsampling \
    --tmp_dir $TMPDIR \
    -w $wd/work/BaseJumper/bj-expression/filter_cells/ \
    -c $wd/config/basejumper.config \
    -profile singularity \
    -resume
)