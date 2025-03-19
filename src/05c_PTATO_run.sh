#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 05c_PTATO_run -o log/%J_05c_PTATO_run.out -e log/%J_05c_PTATO_run.err 'bash src/05c_PTATO_run.sh'

# dirs
wd=$(pwd)

(
  cd out/ptato/filtered_run/

  # run
  module load singularity
  module load ISG/rocker/rver/4.4.0 
  export R_LIBS_USER=$HOME/R-tmp-4.4
  nextflow run $wd/../nextflow/external/PTATO/ptato.nf \
    --out_dir ./ \
    -c $wd/config/ptato_run.config \
    -c $wd/config/ptato.config \
    -w $wd/work/ptato/ \
    --smurf.time '20d' \
    --smurf.cpus 24 \
    --walker.time '7d' \
    --shapeit.reference /lustre/scratch125/casm/team268im/at31/nextflow/external/PTATO/resources/hg38/shapeit/Phasing_reference_no_chr/ \
    --shapeit.maps /lustre/scratch125/casm/team268im/at31/nextflow/external/PTATO/resources/hg38/shapeit/shapeit_maps_no_chr/ \
    --optional.short_variants.somatic_vcfs_dir $wd/out/ptato/intermediate/short_variants/somatic_vcfs/ \
    --optional.short_variants.phased_vcfs_dir $wd/out/ptato/intermediate/short_variants/shapeit/ \
    --optional.short_variants.ab_tables_dir $wd/out/ptato/intermediate/short_variants/ab/ \
    --optional.short_variants.context_beds_dir $wd/out/ptato/intermediate/short_variants/context/ \
    --optional.short_variants.features_beds_dir $wd/out/ptato/intermediate/short_variants/features/ \
    -resume \
    -with-tower \
    -N at31@sanger.ac.uk
)


