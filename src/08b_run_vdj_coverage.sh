#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M20000 -R 'span[hosts=1] select[mem>20000] rusage[mem=20000]' -J vdj_mosdepth -o log/%J_vdj_mosdepth.out -e log/%J_vdj_mosdepth.err 'bash src/08b_run_MOSDEPTH_on_VDJ.sh'

# run mosdepth
wd=$(pwd)
(
  cd out/vdj_reconstruction/
  nextflow run $wd/../nextflow/nf-mosdepth \
    --samplesheet samplesheet.csv \
    --bait_set ig_tcr_genes.bed \
    --fasta /lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa \
    --location local \
    --out_dir ./ \
    -w $wd/work/mosdepth/ \
    -resume \
    -N at31@sanger.ac.uk \
    -with-tower
)

# intersect per-base coverage bed with intervals
module load bedtools2-2.29.0/python-3.10.10
module load tabix/1.18
for bed in out/vdj_reconstruction/PD63118_49882/*/mosdepth/*.per-base.bed.gz ; do
  echo $bed
  bedtools intersect \
    -a <(zcat $bed) \
    -b out/vdj_reconstruction/ig_tcr_genes.bed \
    -wa -wb \
  | bgzip > ${bed/.bed.gz/.regions.bed.gz}
  echo "Done!"
done