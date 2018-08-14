fq_1=$1
fq_2=$2
panel=$3
out_dir=$4
barcode_length=8
gtime="gtime"
calib="./calib"
rainbow="aux/other_tools/rainbow/rainbow"
run_rainbow="aux/run_rainbow.sh"
convert_rainbow="aux/convert_rainbow_to_cluster.sh"
calib_cons="consensus/calib_cons"
bam_readcount="consensus/bam-readcount_v0.8.0/bin/bam-readcount"
bwa="consensus/bwa_v0.7.17/bwa"
ref="simulating/genomes/hg38.fa"
samtools="samtools"
sinvict="consensus/sinvict/sinvict"

# mkdir -p $out_dir
# # Calib
# $gtime -o $out_dir"/calib.gtime" -v $calib -f $fq_1 -r $fq_2 -o $out_dir"/calib." -l $barcode_length --no-sort -c 4
# Rainbow
# $gtime -o $out_dir"/rainbow.gtime" -v $run_rainbow $rainbow $fq_1 $fq_2 2 "true" "$out_dir/rainbow.out"
# $convert_rainbow "$out_dir/rainbow.out" > "$out_dir/rainbow.cluster"
# rm "$out_dir/rainbow.out"
# # calib_cons
# $calib_cons "$out_dir/calib.cluster"   $fq_1 "$out_dir/calib.1"   $fq_2 "$out_dir/calib.2"
# $calib_cons "$out_dir/rainbow.cluster" $fq_1 "$out_dir/rainbow.1" $fq_2 "$out_dir/rainbow.2"
# # BWA MEM & samtools
# $bwa mem $ref "$out_dir/calib.1.fastq"   "$out_dir/calib.2.fastq"   -t 32 | $samtools sort - -O BAM -@ 16 -m 2G -o "$out_dir/calib.bam"
# $bwa mem $ref "$out_dir/rainbow.1.fastq" "$out_dir/rainbow.2.fastq" -t 32 | $samtools sort - -O BAM -@ 16 -m 2G -o "$out_dir/rainbow.bam"
# $samtools index "$out_dir/calib.bam"
# $samtools index "$out_dir/rainbow.bam"
# bam-readcount
mkdir -p "$out_dir/calib.bam-readcount"
$bam_readcount -l $panel -f $ref -w 1 "$out_dir/calib.bam"   > "$out_dir/calib.bam-readcount/out.tsv"
mkdir -p "$out_dir/rainbow.bam-readcount"
$bam_readcount -l $panel -f $ref -w 1 "$out_dir/rainbow.bam" > "$out_dir/rainbow.bam-readcount/out.tsv"
# SiNVICT
mkdir -p "$out_dir/calib.sinvict"
$sinvict -m 3 -t "$out_dir/calib.bam-readcount"   -o "$out_dir/calib.sinvict"
mkdir -p "$out_dir/rainbow.sinvict"
$sinvict -m 3 -t "$out_dir/rainbow.bam-readcount" -o "$out_dir/rainbow.sinvict"
