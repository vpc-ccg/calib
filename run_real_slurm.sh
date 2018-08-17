fq_1=$1
fq_2=$2
panel=$3
out_dir=${4%/}
barcode_length=8
# Pipeline components
gtime="gtime"
calib_cons="consensus/calib_cons"
bam_readcount="consensus/bam-readcount_v0.8.0/bin/bam-readcount"
bwa="consensus/bwa_v0.7.17/bwa"
ref="simulating/genomes/hg38.fa"
samtools="samtools"
sinvict="consensus/sinvict/sinvict"
# Clustering tools and scripts
calib="./calib"
rainbow="aux/other_tools/rainbow/rainbow"
run_rainbow="aux/run_rainbow.sh"
convert_rainbow="aux/convert_rainbow_to_cluster.sh"
run_umi_tools="aux/run_umi-tools.sh"
convert_umitools="aux/convert_umitools_to_cluster.sh"
restore_cluster_missing_reads="aux/restore_cluster_missing_reads.sh"

mkdir -p $out_dir
for tool in calib rainbow umi-tools raw
do
    # Slurm header
    echo "Preparing things for $tool"
    echo "#!/bin/bash"                                  > "$out_dir/$tool.pbs"
    echo "#SBATCH --job-name=$tool"                     >> "$out_dir/$tool.pbs"
    echo "#SBATCH -c 32"                                >> "$out_dir/$tool.pbs"
    echo "#SBATCH --mem 512G"                           >> "$out_dir/$tool.pbs"
    echo "#SBATCH -t 23:59:59"                          >> "$out_dir/$tool.pbs"
    echo "#SBATCH --output=$out_dir/$tool.pbs.out"      >> "$out_dir/$tool.pbs"
    echo "#SBATCH --error=$out_dir/$tool.pbs.err"       >> "$out_dir/$tool.pbs"
    echo "#SBATCH --export=all"                         >> "$out_dir/$tool.pbs"
    echo "#SBATCH -p debug,express,normal,big-mem,long" >> "$out_dir/$tool.pbs"
    # Tool specific commands
    case $tool in
    "calib")
    echo "$gtime -o $out_dir/calib.gtime -v $calib -f $fq_1 -r $fq_2 -o $out_dir/calib. -l $barcode_length --no-sort" >> "$out_dir/$tool.pbs"
    ;;
    "rainbow")
    echo "$gtime -o $out_dir/rainbow.gtime -v $run_rainbow $rainbow $fq_1 $fq_2 2 true $out_dir/rainbow.out" >> "$out_dir/$tool.pbs"
    echo "$convert_rainbow $out_dir/rainbow.out > $out_dir/rainbow.cluster.temp"                             >> "$out_dir/$tool.pbs"
    echo "$restore_cluster_missing_reads $out_dir/rainbow.cluster.temp $fq_1 > $out_dir/rainbow.cluster"     >> "$out_dir/$tool.pbs"
    ;;
    "umi-tools")
    echo "$gtime -o $out_dir/umi-tools.gtime -v $run_umi_tools $fq_1 $fq_2 $out_dir $out_dir/umi-tools_work $bwa $ref $samtools" >> "$out_dir/$tool.pbs"
    echo "$convert_umitools $out_dir/umi-tools.out > $out_dir/umi-tools.cluster.temp"                                            >> "$out_dir/$tool.pbs"
    echo "$restore_cluster_missing_reads $out_dir/umi-tools.cluster.temp $fq_1 > $out_dir/umi-tools.cluster"                     >> "$out_dir/$tool.pbs"
    ;;
    esac
    # Consensus, mapping, and SNV calling
    if [ $tool == "raw" ]
    then
        echo "$bwa mem $ref $fq_1 $fq_2 -t 32 | $samtools sort - -O BAM -@ 32 -m 2G -o $out_dir/$tool.bam" >> "$out_dir/$tool.pbs"
        echo "$samtools index $out_dir/$tool.bam" >> "$out_dir/$tool.pbs"
        echo "mkdir -p $out_dir/$tool.bam-readcount" >> "$out_dir/$tool.pbs"
        echo "$bam_readcount -l $panel -f $ref -w 1 $out_dir/$tool.bam > $out_dir/$tool.bam-readcount/out.tsv" >> "$out_dir/$tool.pbs"
        echo "$sinvict -m 55 -t $out_dir/umi-tools.bam-readcount -o $out_dir/umi-tools.sinvict" >> "$out_dir/$tool.pbs"
    else
        echo "$calib_cons $out_dir/$tool.cluster $fq_1 $out_dir/$tool.1 $fq_2 $out_dir/$tool.2" >> "$out_dir/$tool.pbs"
        echo "$bwa mem $ref $out_dir/$tool.1.fastq $out_dir/$tool.2.fastq -t 32 | $samtools sort - -O BAM -@ 32 -m 2G -o $out_dir/$tool.bam" >> "$out_dir/$tool.pbs"
        echo "$samtools index $out_dir/$tool.bam" >> "$out_dir/$tool.pbs"
        echo "mkdir -p $out_dir/$tool.bam-readcount" >> "$out_dir/$tool.pbs"
        echo "$bam_readcount -l $panel -f $ref -w 1 $out_dir/$tool.bam > $out_dir/$tool.bam-readcount/out.tsv" >> "$out_dir/$tool.pbs"
        echo "$sinvict -m 3 -t $out_dir/umi-tools.bam-readcount -o $out_dir/umi-tools.sinvict" >> "$out_dir/$tool.pbs"
    fi
    sbatch "$out_dir/$tool.pbs"
done
