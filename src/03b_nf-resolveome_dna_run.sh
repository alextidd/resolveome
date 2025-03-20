#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 03b_nf-resolveome_dna_run -o log/%J_03b_nf-resolveome_dna_run.out -e log/%J_03b_nf-resolveome_dna_run.err 'bash src/03b_nf-resolveome_dna_run.sh'

# dirs
wd=$(pwd)

# run
(
  cd out/nf-resolveome/dna/
  nextflow run $wd/../nextflow/nf-resolveome \
    --samplesheet samplesheet.csv \
    --bait_set $wd/data/immune_panel/Probes_merged_ok_combined_Sanger_Immune-v1_TE-91661256_hg19.bed \
    --genes TNFRSF14,CD274,LTB,TNFAIP3,TET2,DUSP2,CCR6,DNMT3A,RFTN1,CBL,RASA2,CXCR3,ACTG1,KLHL6 \
    --fasta /lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa \
    --location local \
    --bamtofastq false \
    --out_dir ./ \
    -w $wd/work/nf-resolveome/dna/ \
    -resume \
    -N at31@sanger.ac.uk \
    -with-tower
)

# create NR and NV matrices
mkdir -p out/nf-resolveome/dna/PD63118/genotyping/
files=(out/nf-resolveome/dna/PD63118/*/genotyping/*_genotyped_mutations.tsv)

# initiate NV and NR matrices with mut ids as rownames
awk -F"\t" 'NR == 1 {print ""} ; NR > 1 && $21 == "nanoseq_mutations" {print $1"-"$2"-"$3"-"$4}' ${files[0]} \
> out/nf-resolveome/dna/PD63118/genotyping/NV_genotyped_mutations.tsv
cp out/nf-resolveome/dna/PD63118/genotyping/NV_genotyped_mutations.tsv \
  out/nf-resolveome/dna/PD63118/genotyping/NR_genotyped_mutations.tsv

# add NV and NR per sample
for file in ${files[@]} ; do
  id=$(echo $file | cut -d/ -f5)
  echo $id
  awk -F"\t" -v id=$id 'NR == 1 {print id} ; NR > 1 && $21 == "nanoseq_mutations" {print $5}' $file \
  | paste out/nf-resolveome/dna/PD63118/genotyping/NR_genotyped_mutations.tsv - \
  > temp && mv temp out/nf-resolveome/dna/PD63118/genotyping/NR_genotyped_mutations.tsv
  awk -F"\t" -v id=$id 'NR == 1 {print id} ; NR > 1 && $21 == "nanoseq_mutations" {print $7}' $file \
  | paste out/nf-resolveome/dna/PD63118/genotyping/NV_genotyped_mutations.tsv - \
  > temp && mv temp out/nf-resolveome/dna/PD63118/genotyping/NV_genotyped_mutations.tsv
done