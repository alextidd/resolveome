#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M20000 -R 'span[hosts=1] select[mem>20000] rusage[mem=20000]' -J 04h_BaseJumper_somatic-variantcalling_dnahyb_filter_cells_run -o log/%J_04h_BaseJumper_somatic-variantcalling_dnahyb_filter_cells_run.out -e log/%J_04h_BaseJumper_somatic-variantcalling_dnahyb_filter_cells_run.err 'bash src/04h_BaseJumper_somatic-variantcalling_dnahyb_filter_cells_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling/sentieon_eval.lic
export LSB_EXCLUSIVE=Y

# run
(
  cd out/BaseJumper/bj-somatic-variantcalling/filter_cells/dnahyb/
  
  nextflow run $wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling \
    --input_csv samplesheet.csv \
    --publish_dir PD63118 \
    --variant_workflow_type somatic_heuristic_filter \
    --dnascope_model_selection bioskryb129 \
    --skip_variant_annotation false \
    --skip_sigprofile false \
    -c $wd/config/bj-somatic-variantcalling_dnahyb.config \
    -c $wd/config/bj-somatic-variantcalling.config \
    -c $wd/config/basejumper.config \
    -w $wd/work/BaseJumper/bj-somatic-variantcalling/filter_cells/dnahyb/ \
    -profile singularity \
    --architecture "x86_64" \
    -resume \
    -with-tower \
    -N at31@sanger.ac.uk
)

# bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J symlink_replace -o "log/%J_cp_basejumper_symlinks.out" "find out/BaseJumper/bj-somatic-variantcalling/filter_cells/dnahyb/PD63118_250321_170558 -type l | while read link; do target=\$(readlink -f \"\$link\"); echo \"Replacing symlink: \$link -> \$target\"; if [ -e \"\$target\" ]; then cp -a \"\$target\" \"\$link.tmp\" && mv \"\$link.tmp\" \"\$link\"; else echo \"Warning: Target does not exist for \$link\"; fi; done"