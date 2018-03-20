base=$1
t=$2
ref=$3

umi_tools extract -I "$base".R1.fastq --bc-pattern=NNNNNNNN --bc-pattern2=NNNNNNNN --read2-in "$base".R2.fastq --stdout "$base".R1.extracted.fastq  --read2-out "$base".R2.extracted.fastq

awk '{if (NR %2 == 0) {print substr($0,2)} else {print}}' "$base".R1.extracted.fastq > "$base".R1.extracted.trim.fastq
awk '{if (NR %2 == 0) {print substr($0,2)} else {print}}' "$base".R2.extracted.fastq > "$base".R2.extracted.trim.fastq

bwa mem -t $2 $3 "$base".R1.extracted.trim.fastq "$base".R2.extracted.trim.fastq | \
awk 'BEGIN{OFS = "\t"; n = 0} (substr($0,1,1)=="@") || (and($2,3852) && $7=="=") {if ($1 in rids) {$1=rids[$1]""substr($1,length($1)-16)} else {rids[$1] = n++} print}' > "$base".extracted.trim.mapping.renamed.sam


samtools sort -@ $2 -O bam -m 2G -T samtemp -o "$base".extracted.trim.mapping.sorted.bam "$base".extracted.trim.mapping.renamed.sam
samtools index "$base".extracted.trim.mapping.sorted.bam
umi_tools group -I "$base".extracted.trim.mapping.sorted.bam --group-out "$base".extracted.trim.mapping.sorted.cluster.tsv --paired --method cluster --edit-distance-threshold 2


awk 'and($2,64)  {print "@"n++"\n"substr($1, length($1) - 15, 8)$10"\n+\n"$11}' "$base".extracted.trim.mapping.renamed.sam >  "$base".extracted.trim.mapping.R1.fastq
awk 'and($2,128) {print "@"n++"\n"substr($1, length($1) - 15, 8)$10"\n+\n"$11}' "$base".extracted.trim.mapping.renamed.sam >  "$base".extracted.trim.mapping.R2.fastq
