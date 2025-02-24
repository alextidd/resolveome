#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 04f_BaseJumper_expression_run -o log/%J_04f_BaseJumper_expression_run.out -e log/%J_04f_BaseJumper_expression_run.err 'bash src/04f_BaseJumper_expression_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# run bj-expression
(
  cd out/BaseJumper/bj-expression/
  nextflow run $wd/../nextflow/external/BaseJumper/bj-expression \
    --input_csv samplesheet.csv \
    --publish_dir ./ \
    --genome GRCh38 \
    --tmp_dir $TMPDIR \
    -w $wd/work/BaseJumper/bj-expression/ \
    -c $wd/config/bj-expression.config \
    -profile singularity \
    --process.containerOptions '--bind /lustre,/nfs,/data,/software'
)