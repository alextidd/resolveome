
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

NFCORE_BAMTOFASTQ:BAMTOFASTQ:PRE_CONVERSION_QC:SAMTOOLS_IDXSTATS (plate1_wellC2_dna_run49686)

for id in plate1_wellC4_dna_run49686 plate1_wellE3_dna_run49686 plate3_wellA2_dna_run49882 plate3_wellB2_dna_run49882 plate3_wellC2_dna_run49882 ; do
  (
    echo $id
    cd $id/bam
    samtools index $id.bam
  )
done