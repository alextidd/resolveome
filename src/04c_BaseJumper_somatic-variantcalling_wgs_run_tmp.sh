#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J 04c_BaseJumper_somatic-variantcalling_wgs_run_tmp -o log/%J_04c_BaseJumper_somatic-variantcalling_wgs_run_tmp.out -e log/%J_04c_BaseJumper_somatic-variantcalling_wgs_run_tmp.err 'bash src/04c_BaseJumper_somatic-variantcalling_wgs_run_tmp.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling/sentieon_eval.lic
export LSB_EXCLUSIVE=Y

# run
(
  cd out/BaseJumper/bj-somatic-variantcalling/wgs/

  cat samplesheet.csv |
  head -1 samplesheet.csv > samplesheet_tmp.csv
  grep plate3_wellA12_dna_run49882 samplesheet.csv >> samplesheet_tmp.csv
  
  nextflow run $wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling \
    --input_csv samplesheet_tmp.csv \
    --publish_dir PD63118 \
    --variant_workflow_type somatic_heuristic_filter \
    -c $wd/config/bj-somatic-variantcalling.config \
    -c $wd/config/basejumper.config \
    -w $wd/work/BaseJumper/bj-somatic-variantcalling/wgs/ \
    -profile singularity \
    --architecture "x86_64" \
    --dnascope_model_selection bioskryb129 \
    -resume \
    -N at31@sanger.ac.uk
)

# bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J symlink_replace -o "log/%J_cp_basejumper_symlinks.out" "find out/BaseJumper/bj-somatic-variantcalling/wgs/PD63118_250306_232001 -type l | while read link; do target=\$(readlink -f \"\$link\"); echo \"Replacing symlink: \$link -> \$target\"; if [ -e \"\$target\" ]; then cp -a \"\$target\" \"\$link.tmp\" && mv \"\$link.tmp\" \"\$link\"; else echo \"Warning: Target does not exist for \$link\"; fi; done"