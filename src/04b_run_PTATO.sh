#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J ptato -o log/ptato_%J.out -e log/ptato_%J.err 'bash src/04b_run_PTATO.sh'

# dirs
wd=$(pwd)
cd out/ptato/

# run
module load singularity
module load ISG/rocker/rver/4.4.0 
export R_LIBS_USER=$HOME/R-tmp-4.4
nextflow run $wd/../tools/PTATO/ptato.nf \
  --run.svs false \
  --run.cnvs false \
  --smurf.time '20d' \
  --smurf.cpus 24 \
  --shapeit.reference /lustre/scratch125/casm/team268im/at31/tools/PTATO/resources/hg38/shapeit/Phasing_reference_no_chr/ \
  --shapeit.maps /lustre/scratch125/casm/team268im/at31/tools/PTATO/resources/hg38/shapeit/shapeit_maps_no_chr/ \
  -c ~/.nextflow/config \
  -c $wd/../tools/PTATO/configs/run-template.config \
  -c $wd/../tools/PTATO/configs/nextflow.config \
  -c $wd/../tools/PTATO/configs/process.config \
  -c $wd/../tools/PTATO/configs/resources.config \
  -c $wd/config/ptato.config \
  --out_dir ./ \
  -w $wd/work/ptato/ \
  -resume \
  -with-tower \
  -N at31@sanger.ac.uk