#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M20000 -R 'span[hosts=1] select[mem>20000] rusage[mem=20000]' -J vdj_coverage -o log/%J_vdj_coverage.out -e log/%J_vdj_coverage.err 'bash src/08b_vdj_coverage_run.sh'

# run mosdepth
wd=$(pwd)
(
  cd out/vdj_coverage/
  nextflow run $wd/../nextflow/nf-mosdepth \
    --samplesheet samplesheet.csv \
    --bait_set regions/ig_tcr_genes.bed \
    --fasta /lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa \
    --location local \
    --out_dir ./ \
    -w $wd/work/mosdepth/ \
    -resume \
    -N at31@sanger.ac.uk \
    -with-tower
)

# intersect per-base coverage bed with genes and regions
module load bedtools2-2.29.0/python-3.10.10
module load tabix/1.18

for bed in out/vdj_coverage/PD63118/*/mosdepth/*.per-base.bed.gz ; do
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