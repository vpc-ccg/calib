fq_1=$1
fq_2=$2
out_path=$3
working_dir=$4
bwa=$5
ref=$6
samtools=$7

mkdir -p $working_dir
umi_tools extract \
    -I $fq_1 --read2-in $fq_2 \
    --bc-pattern=NNNNNNNN --bc-pattern2=NNNNNNNN \
    --stdout "$working_dir/barcode-extracted.1.fastq"  --read2-out "$working_dir/barcode-extracted.2.fastq"

bwa mem -t 32 $ref "$working_dir/barcode-extracted.1.fastq" "$working_dir/barcode-extracted.2.fastq" | \
    awk 'substr($0,1,1) == "@" || $2 < 256' | \
    $samtools sort - -O BAM -@ 32 -o "$working_dir/barcode-extracted.bam"
$samtools index "$working_dir/barcode-extracted.bam"
umi_tools group  \
    --paired \
    --method cluster \
    --edit-distance-threshold 2 \
    --soft-clip-threshold 9 \
    -I "$working_dir/barcode-extracted.bam" \
    --group-out "$out_path"
