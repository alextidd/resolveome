#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 06c_vdj_coverage_filter_cells_run -o log/%J_06c_vdj_coverage_filter_cells_run.out -e log/%J_06c_vdj_coverage_filter_cells_run.err 'bash src/06c_vdj_coverage_filter_cells_run.sh'

# run mosdepth
wd=$(pwd)
(
  cd out/vdj_coverage/filter_cells/
  nextflow run $wd/../nextflow/nf-mosdepth \
    --samplesheet samplesheet.csv \
    --bait_set ../regions/ig_tcr_genes.bed \
    --fasta /lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCh37d5/genome.fa \
    --location local \
    --out_dir ./ \
    -w $wd/work/mosdepth/filter_cells/ \
    -resume \
    -N at31@sanger.ac.uk \
    -with-tower
)

# intersect per-base coverage bed with genes and regions
module load bedtools2-2.29.0/python-3.10.10
module load tabix/1.18

for bed in out/vdj_coverage/filter_cells/PD63118/*/mosdepth/*.per-base.bed.gz ; do
  echo $bed
  bedtools intersect \
    -a <(zcat $bed) \
    -b out/vdj_coverage/regions/ig_tcr_genes.bed \
    -wa -wb \
  | bgzip > ${bed/.bed.gz/.ig_tcr_genes.bed.gz}
  bedtools intersect \
    -a <(zcat $bed) \
    -b out/vdj_coverage/regions/ig_tcr_regions.bed \
    -wa -wb \
  | bgzip > ${bed/.bed.gz/.ig_tcr_regions.bed.gz}
  echo "Done!"
done