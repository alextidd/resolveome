#!/bin/bash

# dirs
mkdir out/nf-caveman/

# create samplesheet
ls /lustre/scratch125/casm/team268im/fa8/117/PTA_49686/lane4-5/plex{1..19}/*cram | \
grep -v 'phix' | \
awk '
BEGIN {
    FS = "/"; OFS = ",";
    print "donor_id,sample_id,tumour_bam,normal_bam";
}
{
    split($NF, parts, "#");                 # Split filename to extract plex number
    plex_num = parts[2];                    # Get the plex number from the filename
    donor_id = "PD63118";                   # Static donor ID
    sample_id = donor_id "_plex" plex_num;  # Generate sample_id
    sub(/\.cram$/, "", sample_id);          # Remove .cram suffix
    tumour_bam = $0;                        # Full path of the tumour BAM
    normal_bam = "/lustre/scratch125/casm/team268im/fa8/119/live/3464/PD63118b_lo0001/PD63118b_lo0001.sample.dupmarked.bam";
    print donor_id, sample_id, tumour_bam, normal_bam;
}' > out/nf-caveman/samplesheet.csv

# convert to bams
