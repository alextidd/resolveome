#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M50000 -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]' -J get_het_snps -o log/get_het_snps_%J.out -e log/get_het_snps_%J.err 'bash src/03b_nf-resolveome_get_het_snps.sh'

# modules
module load bcftools-1.9/python-3.11.6

# dirs
outdir=out/nf-resolveome/PD63118/snps/
mkdir $outdir

# run
bcftools view -a -s PD63118  \
  /lustre/scratch125/casm/team268im/fa8/117/PTA_49686/vcfs_dir/donor1/final.all.sorted.vcf.gz |
bcftools view -g "het" - |
gzip > $outdir/final.all.sorted.vcf.gz
