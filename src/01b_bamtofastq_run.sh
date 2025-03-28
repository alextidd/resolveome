
#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 01b_bamtofastq_run -o log/%J_01b_bamtofastq_run.out -e log/%J_01b_bamtofastq_run.err 'bash src/01b_bamtofastq_run.sh'

# dirs
wd=$(pwd)

# modules
module load singularity/3.11.4

# run
(
  cd out/bamtofastq/
  cat samplesheet.csv | grep "^sample_id\|^plate10_" > samplesheet.tmp.csv
  nextflow run nf-core/bamtofastq \
    -profile singularity,sanger \
    --input samplesheet.tmp.csv \
    --outdir . \
    -w $wd/work/bamtofastq/ \
    -resume
)