#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J nf-gatk -o log/%J_nf-gatk.out -e log/%J_nf-gatk.err 'bash src/02b_nf-gatk_run.sh'

# dir
wd=$(pwd)
export TMPDIR=/lustre/scratch125/casm/team268im/at31/tmp/bash/

# modules
module load singularity/3.11.4

# run
(
  cd out/nf-gatk/
  nextflow run $wd/../nextflow/nf-gatk/ \
    --samplesheet samplesheet.csv \
    --outdir ./ \
    --fasta /lustre/scratch125/casm/team268im/fa8/117/PTA_49686/PTATO/resources/hg38/Homo_sapiens.GRCh38.dna_sm.toplevel.fa \
    --realign \
    -w $wd/work/nf-gatk/ \
    -N at31@sanger.ac.uk \
    -resume
)

# # modules
# module load pcap-core
# module add gatk/4.5.0.0

# # files and dirs
# wd=$(pwd)
# fasta=/lustre/scratch125/casm/team268im/fa8/117/PTA_49686/PTATO/resources/hg38/Homo_sapiens.GRCh38.dna_sm.toplevel.fa
# outdir=$wd/out/resolveome/bams/grch38/

# # remap to GRCh38
# while IFS="," read -r run plex_n lane_n seq_type bam donor_id id ; do
#   echo "$donor_id $id"
#   bsub \
#     -q long  -e $wd/log/%J_remap_$id.err -o $wd/log/%J_remap_$id.out -J remap_$id \
#     -M 64000 -R "select[mem>64000] rusage[mem=64000] span[hosts=1]" -n 16 \
#     bwa_mem.pl -b '-T 30 -Y' -o $outdir/$id/ -reference $fasta -sample $id -fragment 10000000 -threads 16 $bam
# done < <(sed 1d data/resolveome/DNA/samplesheet_local.csv)

# # link matched normal
# mkdir $outdir/PD63118/
# ln -s \
#   /lustre/scratch125/casm/team268im/fa8/117/PTA_49686/PD63118_MATCHED_NORMAL/PD63118b_lo0001.sample.dupmarked.bam \
#   $outdir/PD63118/PD63118.bam

# # run GATK
# while IFS="," read -r run plex_n lane_n seq_type bam donor_id id ; do
#   echo "$donor_id $id"
#   bsub -q long \
#     -M 10000 -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" \
#     -e log/%J_gatk_haplotypecaller_$id.err \
#     -o log/%J_gatk_haplotypecaller_$id.out \
#   gatk --java-options -Xmx8500m HaplotypeCaller -R $fasta -I $outdir/$id.grch38.bam -O r$outdir/$id.grch38.gvcf -ERC GVCF
# done

#  bsub10000 -q long -e gatk.matchednormal.err -o gatk.matchednormal.err gatk --java-options -Xmx8500m HaplotypeCaller \
#   -R $fasta \
#   -I PD63118_MATCHED_NORMAL/PD63118.grch38.bam/PD63118.bam \
#   -O PD63118_MATCHED_NORMAL/PD63118.grch38.bam/PD63118.gvcf \
#   -ERC GVCF

# # Then, combine all your gvcf files:

# # By chr:
# for chr in 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT
#   do echo $chr
#   cd /lustre/scratch125/casm/team268im/fa8/117/PTA_49686
#   module add gatk/4.5.0.0
#   JAVA_TOOL_OPTIONS=-Xmx78000m
#   bsub80000 -q long -e gatk.combine.basement.$chr.err -o gatk.combine.basement.$chr.err gatk --java-options -Xmx78000m CombineGVCFs \
#    -R $fasta \
#    --variant PD63118_MATCHED_NORMAL/PD63118.grch38.bam/PD63118.gvcf \
#    --variant remapped_bams/10.grch38.bam/PLEX10.gvcf \
#    --variant remapped_bams/11.grch38.bam/PLEX11.gvcf \
#    --variant remapped_bams/12.grch38.bam/PLEX12.gvcf \
#    --variant remapped_bams/13.grch38.bam/PLEX13.gvcf \
#    --variant remapped_bams/14.grch38.bam/PLEX14.gvcf \
#    --variant remapped_bams/15.grch38.bam/PLEX15.gvcf \
#    --variant remapped_bams/16.grch38.bam/PLEX16.gvcf \
#    --variant remapped_bams/17.grch38.bam/PLEX17.gvcf \
#    --variant remapped_bams/18.grch38.bam/PLEX18.gvcf \
#    --variant remapped_bams/19.grch38.bam/PLEX19.gvcf \
#    --variant remapped_bams/1.grch38.bam/PLEX1.gvcf \
#    --variant remapped_bams/2.grch38.bam/PLEX2.gvcf \
#    --variant remapped_bams/3.grch38.bam/PLEX3.gvcf \
#    --variant remapped_bams/4.grch38.bam/PLEX4.gvcf \
#    --variant remapped_bams/5.grch38.bam/PLEX5.gvcf \
#    --variant remapped_bams/6.grch38.bam/PLEX6.gvcf \
#    --variant remapped_bams/7.grch38.bam/PLEX7.gvcf \
#    --variant remapped_bams/8.grch38.bam/PLEX8.gvcf \
#    --variant remapped_bams/9.grch38.bam/PLEX9.gvcf \
#    -L $chr \
#    -O combined.g_vBasement.$chr.vcf.gz
# done


# # Finally, perform joint-genotyping on the combined gvcf:
# for chr in 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT
#   do echo $chr
#   cd /lustre/scratch125/casm/team268im/fa8/117/PTA_49686
#   module add gatk/4.5.0.0
#   bsub80000 -q normal -e gatk.call.$chr.err -o gatk.call.$chr.err gatk --java-options -Xmx78g GenotypeGVCFs \
#    -R $fasta \
#    -V combined.g_vBasement.$chr.vcf.gz \
#    -O final.$chr.vcf.gz
# done

# # Now concatenate final vcfs with vcf-concat: vcf-concat is useless!!!!!
# # vcf-concat [OPTIONS] A.vcf.gz B.vcf.gz C.vcf.gz > out.vcf
# # vcf-concat final*vcf.gz > final.all.vcf
# bcftools concat -o final.all.vcf final*vcf.gz 
# bcftools sort -o final.all.sorted.vcf -O v final.all.vcf
