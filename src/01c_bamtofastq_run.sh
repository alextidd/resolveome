
#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 01c_bamtofastq_run -o log/%J_01c_bamtofastq_run.out -e log/%J_01c_bamtofastq_run.err 'bash src/01c_bamtofastq_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity/3.11.4

# run
(
  cd out/bamtofastq/
  nextflow run nf-core/bamtofastq \
    -profile singularity,sanger \
    --input samplesheet.csv \
    --outdir . \
    -w $wd/work/bamtofastq/ \
    -resume
)

# run multiqc
mkdir out/multiqc/
multiqc out/bamtofastq/ \
  -f \
  -o out/multiqc/dna \
  --ignore-samples "*rna*" --ignore-samples "*dnayb*"

multiqc out/bamtofastq/ \
  -f \
  -o out/multiqc/dnahyb \
  --ignore-samples "_rna_" \
  --ignore-samples "_dna_"
multiqc out/bamtofastq/ \
  -o out/multiqc/rna \
  --ignore-samples "_dna_" \
  --ignore-samples "_dnahyb_"