#!/bin/bash
# 1 = dunovo_temp_directory
# 2 = forward_reads
# 3 = reverse_reads
# 4 = dunovo_prefix
# 5 = dunovo_invariant
# 6 = barcode_length
# 7 = dunovo_dist
# 8 = dunovo_output

# Making families
source activate bargoat_dunovo
samtools sort
rm -rf $1
mkdir -p "$1"/refdir
paste $2 $3 | \
    paste - - - - | \
    awk -f "$4"/make-barcodes.awk -v INVARIANT="$5" TAG_LEN="$6" | \
    sort > "$1"/uncorrected_families.tsv;
# Aligning families
"$4"/baralign.sh "$1"/uncorrected_families.tsv "$1"/refdir "$1"/temp.bam;
echo ALL GOOD
# Correcting barcodes
samtools view -f 256 "$1"/temp.bam | \
    "$4"/correct.py -d $7 "$1"/uncorrected_families.tsv "$1"/refdir/barcodes.fa | \
    sort > $8
# Cleaning intermediate files...
rm -rf $1
source deactivate
