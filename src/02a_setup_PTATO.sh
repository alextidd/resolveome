#!/bin/bash
cd /lustre/scratch125/casm/team268im/at31/nextflow/external/PTATO

# clone the repo
git clone git@github.com:ToolsVanBox/PTATO.git

# download the sif
module load singularity
singularity pull ptato_1.2.0.sif docker://vanboxtelbioinformatics/ptato:1.2.0

# open tars
tar -xzvf resources/hg38/gripss/gridss_pon_breakpoint.tar.gz -C resources/hg38/gripss/
tar -xzvf resources/hg38/cobalt/COBALT_PTA_Normalized_Full.tar.gz -C resources/hg38/cobalt/
tar -xzvf resources/hg38/smurf/Mutational_blacklists/Fetal_15x_raw_variants_hg38.tar.gz -C resources/hg38/smurf/Mutational_blacklists/
tar -xzvf resources/hg38/smurf/Mutational_blacklists/MSC_healthyBM_raw_variants_hg38.tar.gz -C resources/hg38/smurf/Mutational_blacklists/

# get reference genome
ln -s \
  /lustre/scratch125/casm/team268im/fa8/117/PTA_49686/PTATO/resources/hg38/Homo_sapiens.GRCh38.dna_sm.toplevel.* \
  resources/hg38/

# download shapeit phasing reference
mkdir -p resources/hg38/shapeit/Phasing_reference/
wget -qO- https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/ | \
grep -oP '(?<=href=")[^"]*\.(vcf\.gz|vcf\.gz\.tbi)' | \
sed "s|^|https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/|" \
> resources/hg38/shapeit/Phasing_reference/README.txt
wget \
  -i resources/hg38/shapeit/Phasing_reference/README.txt \
  -P resources/hg38/shapeit/Phasing_reference/

# download shapeit maps
mkdir resources/hg38/shapeit/shapeit_maps/
wget https://github.com/odelaneau/shapeit4/blob/master/maps/genetic_maps.b38.tar.gz
tar -xzvf genetic_maps.b38.tar.gz -C resources/hg38/shapeit/shapeit_maps/
rm genetic_maps.b38.tar.gz