#!/bin/bash
#SBATCH --job-name=real_data
#SBATCH -c 32
#SBATCH --mem 512G
#SBATCH -t 23:59:59
#SBATCH --output=real_data.out
#SBATCH --error=real_data.err
#SBATCH --export=all
#SBATCH -p debug,express,normal,big-mem,long
fq_1="real-run-9/1.fastq"
fq_2="real-run-9/2.fastq"
panel="real-run-9/panel.bed"
out_dir="real-run-9_auto"
barcode_length=8
gtime="gtime"
calib="./calib"
rainbow="aux/other_tools/rainbow/rainbow"
run_rainbow="aux/run_rainbow.sh"
run_umi_tools="aux/run_umi-tools.sh"
convert_rainbow="aux/convert_rainbow_to_cluster.sh"
convert_umitools="aux/convert_umitools_to_cluster.sh"
restore_cluster_missing_reads="aux/restore_cluster_missing_reads.sh"
calib_cons="consensus/calib_cons"
bam_readcount="consensus/bam-readcount_v0.8.0/bin/bam-readcount"
bwa="consensus/bwa_v0.7.17/bwa"
ref="simulating/genomes/hg38.fa"
samtools="samtools"
sinvict="consensus/sinvict/sinvict"

mkdir -p $out_dir
# Calib
$gtime -o $out_dir"/calib.gtime" -v $calib -f $fq_1 -r $fq_2 -o $out_dir"/calib." -l $barcode_length --no-sort
# Rainbow
$gtime -o $out_dir"/rainbow.gtime" -v $run_rainbow $rainbow $fq_1 $fq_2 2 "true" "$out_dir/rainbow.out"
$convert_rainbow "$out_dir/rainbow.out" > "$out_dir/rainbow.cluster.temp"
$restore_cluster_missing_reads "$out_dir/rainbow.cluster.temp" $fq_1 > "$out_dir/rainbow.cluster"
rm "$out_dir/rainbow.out" "$out_dir/rainbow.cluster.temp"
# # umi-tools
echo $gtime -o $out_dir"/umi-tools.gtime" -v $run_umi_tools $fq_1 $fq_2 $out_dir $out_dir"/umi_tools_work" $bwa $ref $samtools
$convert_umitools "$out_dir/umi-tools.out" > "$out_dir/umi-tools.cluster.temp"
$restore_cluster_missing_reads "$out_dir/umi-tools.cluster.temp" $fq_1 > "$out_dir/umi-tools.cluster"
rm "$out_dir/umi-tools.out" "$out_dir/umi-tools.cluster.temp"
# calib_cons
$calib_cons "$out_dir/calib.cluster"   $fq_1 "$out_dir/calib.1"   $fq_2 "$out_dir/calib.2"
$calib_cons "$out_dir/rainbow.cluster" $fq_1 "$out_dir/rainbow.1" $fq_2 "$out_dir/rainbow.2"
$calib_cons "$out_dir/umi-tools.cluster" $fq_1 "$out_dir/umi-tools.1" $fq_2 "$out_dir/umi-tools.2"
# BWA MEM & samtools
$bwa mem $ref "$out_dir/calib.1.fastq"     "$out_dir/calib.2.fastq"     -t 32 | $samtools sort - -O BAM -@ 16 -m 2G -o "$out_dir/calib.bam"
$bwa mem $ref "$out_dir/rainbow.1.fastq"   "$out_dir/rainbow.2.fastq"   -t 32 | $samtools sort - -O BAM -@ 16 -m 2G -o "$out_dir/rainbow.bam"
$bwa mem $ref "$out_dir/umi-tools.1.fastq" "$out_dir/umi-tools.2.fastq" -t 32 | $samtools sort - -O BAM -@ 16 -m 2G -o "$out_dir/umi-tools.bam"
$samtools index "$out_dir/calib.bam"
$samtools index "$out_dir/rainbow.bam"
$samtools index "$out_dir/umi-tools.bam"
#bam-readcount
mkdir -p "$out_dir/calib.bam-readcount"
$bam_readcount -l $panel -f $ref -w 1 "$out_dir/calib.bam"     > "$out_dir/calib.bam-readcount/out.tsv"
mkdir -p "$out_dir/rainbow.bam-readcount"
$bam_readcount -l $panel -f $ref -w 1 "$out_dir/rainbow.bam"   > "$out_dir/rainbow.bam-readcount/out.tsv"
mkdir -p "$out_dir/umi-tools.bam-readcount"
$bam_readcount -l $panel -f $ref -w 1 "$out_dir/umi-tools.bam" > "$out_dir/umi-tools.bam-readcount/out.tsv"
# SiNVICT
mkdir -p "$out_dir/calib.sinvict"
$sinvict -m 3 -t "$out_dir/calib.bam-readcount"     -o "$out_dir/calib.sinvict"
mkdir -p "$out_dir/rainbow.sinvict"
$sinvict -m 3 -t "$out_dir/rainbow.bam-readcount"   -o "$out_dir/rainbow.sinvict"
mkdir -p "$out_dir/umi-tools.sinvict"
$sinvict -m 3 -t "$out_dir/umi-tools.bam-readcount" -o "$out_dir/umi-tools.sinvict"
