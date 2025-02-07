#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M20000 -R 'span[hosts=1] select[mem>20000] rusage[mem=20000]' -J ptato -o log/ptato_%J.out -e log/ptato_%J.err 'bash src/03b_PTATO_run.sh'

# dirs
wd=$(pwd)

(
  cd out/ptato/filtered_run/

  # run
  module load singularity
  module load ISG/rocker/rver/4.4.0 
  export R_LIBS_USER=$HOME/R-tmp-4.4
  nextflow run $wd/../tools/PTATO/ptato.nf \
    --run.svs false \
    --run.cnvs false \
    --smurf.time '20d' \
    --smurf.cpus 24 \
    --walker.time '7d' \
    --shapeit.reference /lustre/scratch125/casm/team268im/at31/tools/PTATO/resources/hg38/shapeit/Phasing_reference_no_chr/ \
    --shapeit.maps /lustre/scratch125/casm/team268im/at31/tools/PTATO/resources/hg38/shapeit/shapeit_maps_no_chr/ \
    --optional.short_variants.somatic_vcfs_dir $wd/out/ptato/intermediate/short_variants/somatic_vcfs/ \
    --optional.short_variants.phased_vcfs_dir $wd/out/ptato/intermediate/short_variants/shapeit/ \
    --optional.short_variants.ab_tables_dir $wd/out/ptato/intermediate/short_variants/ab/ \
    --optional.short_variants.context_beds_dir $wd/out/ptato/intermediate/short_variants/context/ \
    --optional.short_variants.features_beds_dir $wd/out/ptato/intermediate/short_variants/features/ \
    -c ~/.nextflow/config \
    -c $wd/../tools/PTATO/configs/run-template.config \
    -c $wd/config/ptato.config \
    --out_dir ./ \
    -w $wd/work/ptato/ \
    -resume \
    -with-tower \
    -N at31@sanger.ac.uk
)