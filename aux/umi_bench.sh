base=$1
t=$2
ref=$3
samtools_path=$4

umi_tools extract -I "$base".R1.fastq --bc-pattern=NNNNNNNN --bc-pattern2=NNNNNNNN --read2-in "$base".R2.fastq --stdout "$base".R1.extracted.fastq  --read2-out "$base".R2.extracted.fastq

awk '{if (NR %2 == 0) {print substr($0,2)} else {print}}' "$base".R1.extracted.fastq > "$base".R1.extracted.trim.fastq
awk '{if (NR %2 == 0) {print substr($0,2)} else {print}}' "$base".R2.extracted.fastq > "$base".R2.extracted.trim.fastq

bwa mem -t $2 $3 "$base".R1.extracted.trim.fastq "$base".R2.extracted.trim.fastq > "$base".extracted.trim.sam

"$samtools_path"samtools view -F 3852 -h "$base".extracted.trim.sam | \
   awk 'substr($0,1,1)=="@" || $5 > 0' | \
   awk '
BEGIN{OFS = "\t"; n = 0}
{
   if(substr($0,1,1)=="@") {
       print;
   } else if ($7=="=") {
       key = $1;
       $1 = substr($1,length($1)-16);
       if (key in first_mates) {
           print n""first_mates[key];
           print n""$0;
           n++;
       } else {
           first_mates[key] = $0;
       }
   }
}' > \
   "$base".extracted.trim.mapping.renamed.sam


"$samtools_path"samtools sort -@ $2 -O bam -m 2G -T samtemp -o "$base".extracted.trim.mapping.sorted.bam "$base".extracted.trim.mapping.renamed.sam
"$samtools_path"samtools index "$base".extracted.trim.mapping.sorted.bam
umi_tools group -I "$base".extracted.trim.mapping.sorted.bam --group-out "$base".extracted.trim.mapping.sorted.cluster.tsv --paired --method cluster --edit-distance-threshold 2 --soft-clip-threshold 9
awk 'BEGIN{prev = -1} NR > 1 {if ($9 != prev) {prev = $9; print "# "$9} rid = substr($1, 1,length($1)-17); print $2"."$3,rid,"@"rid,"s1","q1","@"rid,"s2","q2" }' "$base".extracted.trim.mapping.sorted.cluster.tsv > "$base".extracted.trim.mapping.sorted.cluster

"$samtools_path"samtools view -f 64  "$base".extracted.trim.mapping.renamed.sam | awk '{print "@"substr($1,1,length($1)-17)"\n"substr($1, length($1) - 15  , 8)$10"\n+\nKKKKKKKK"$11}' > "$base".extracted.trim.mapping.R1.fastq
"$samtools_path"samtools view -f 128 "$base".extracted.trim.mapping.renamed.sam | awk '{print "@"substr($1,1,length($1)-17)"\n"substr($1, length($1) - 15+8, 8)$10"\n+\nKKKKKKKK"$11}' > "$base".extracted.trim.mapping.R2.fastq
