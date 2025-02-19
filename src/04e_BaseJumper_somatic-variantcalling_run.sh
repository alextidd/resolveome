#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 04e_BaseJumper_somatic-variantcalling_run -o log/%J_04e_BaseJumper_somatic-variantcalling_run.out -e log/%J_04e_BaseJumper_somatic-variantcalling_run.err 'bash src/04e_BaseJumper_somatic-variantcalling_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling/Bioskryb_Inc_c8_eval.lic

# run bj-somatic-variantcalling
(
  cd out/BaseJumper/bj-somatic-variantcalling/

  # create samplesheet
  echo "biosampleName,read1,read2,groups,isbulk,bam" > samplesheet.csv
  while read -r bam ; do
    id=$(basename ${bam/.bam/})
    echo "$id,,,PD63118,false,$bam" >> samplesheet.csv
  done < <( ls $wd/out/BaseJumper/bj-wgs/_250217_104202/secondary_analyses/alignment/*.bam )
  
  # run
  # Viren: switch --variant_workflow_type to somatic_heuristic_filter - doesn't 
  # require a matched normal (for bj-somatic-variantcalling)
  nextflow run $wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling \
    --input_csv samplesheet.csv \
    --publish_dir ./ \
    -c ~/.nextflow/config \
    -c $wd/config/bj-somatic-variantcalling.config \
    -c $wd/config/basejumper.config \
    --is_bam \
    -w $wd/work/BaseJumper/bj-somatic-variantcalling/ \
    -profile singularity \
    --architecture "x86_64" \
    --process.containerOptions '--bind /lustre,/nfs,/data,/software' \
    -resume
)