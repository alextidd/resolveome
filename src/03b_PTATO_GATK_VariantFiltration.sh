#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M20000 -R 'span[hosts=1] select[mem>20000] rusage[mem=20000]' -J GATK_varfilt -o log/%J_GATK_varfilt.out -e log/%J_GATK_varfilt.err 'bash src/03b_PTATO_GATK_VariantFiltration.sh'

# dirs
donor_id=PD63118
out_dir=out/GATK/VariantFiltration/$donor_id/

nextflow run ../nextflow/nf-gatk \
  --samplesheet data/resolveome/DNA/samplesheet_local.csv \
  --fasta /lustre/scratch125/casm/team268im/fa8/117/PTA_49686/PTATO/resources/hg38/Homo_sapiens.GRCh38.dna_sm.toplevel.fa \
  -resume


# modules
module load bcftools-1.9/python-3.11.6
module load gatk-4.5.0.0/python-3.12.0

# function
filter_variants() {
  local input_vcf="$1"
  local output_vcf="$2"
  
  gatk VariantFiltration \
    --reference /lustre/scratch125/casm/team268im/fa8/117/PTA_49686/PTATO/resources/hg38/Homo_sapiens.GRCh38.dna_sm.toplevel.fa \
    --variant "$input_vcf" \
    --output "$output_vcf" \
    --filter-name 'SNP_LowQualityDepth' \
    --filter-expression 'QD < 2.0' \
    --filter-name 'SNP_MappingQuality' \
    --filter-expression 'MQ < 40.0' \
    --filter-name 'SNP_StrandBias' \
    --filter-expression 'FS > 60.0' \
    --filter-name 'SNP_MQRankSumLow' \
    --filter-expression 'MQRankSum < -12.5' \
    --filter-name 'SNP_ReadPosRankSumLow' \
    --filter-expression 'ReadPosRankSum < -8.0' \
    --filter-name 'SNP_LowCoverage' \
    --filter-expression 'DP < 5' \
    --filter-name 'SNP_VeryLowQual' \
    --filter-expression 'QUAL <30' \
    --filter-name 'SNP_LowQual' \
    --filter-expression 'QUAL >= 30.0 && QUAL < 50.0' \
    --filter-name 'SNP_SOR' \
    --filter-expression 'SOR > 4.0' \
    --filter-name 'SNP_HardToValidate' \
    --filter-expression 'MQ0 >= 4 && ((MQ0/(1.0 * DP)) > 0.1)' \
	  --filter-name 'SNP_HaplotypeScoreHigh' \
    --filter-expression 'HaplotypeScore > 13.0' \
    --cluster-size 3 \
    --cluster-window-size 10
}

# run
filter_variants \
  /lustre/scratch125/casm/team268im/fa8/117/PTA_49686/vcfs_dir/donor1/final.all.sorted.vcf.gz \
  $out_dir/final.all.sorted.vcf.gz

# filter to pass only
bcftools view \
  -f PASS \
  $out_dir/final.all.sorted.vcf.gz \
  -Oz \
  -o $out_dir/final.all.sorted.filtered.vcf.gz
bcftools index $out_dir/final.all.sorted.filtered.vcf.gz
tabix -p vcf $out_dir/final.all.sorted.filtered.vcf.gz

# save filtering summary
echo -e "filter\tn" > $out_dir/filtering_summary.tsv
zcat $out_dir/final.all.sorted.vcf.gz |
grep -v '#' | cut -f7 | sort | uniq -c | awk '{print $2 "\t" $1}' \
>> $out_dir/filtering_summary.tsv

# stage for PTATO pipeline
mkdir -p out/ptato/vcfs/$donor_id/
cp $out_dir/final.all.sorted.filtered.vcf.gz* out/ptato/vcfs/$donor_id/