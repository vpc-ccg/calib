fq_1=$1
fq_2=$2
panel=$3
out_dir=${4%/}
barcode_length=8
# Pipeline components
gtime="aux/other_tools/time-1.9/time"
calib_cons_directory="calib_cons/"
calib_cons="consensus/calib_cons"
bam_readcount_directory="consensus/bam-readcount_v0.8.0/"
bam_readcount="consensus/bam-readcount_v0.8.0/bin/bam-readcount"
bwa_directory="aux/other_tools/bwa/"
bwa="aux/other_tools/bwa/bwa"
ref="simulating/genomes/hg19.fa"
samtools="samtools"
sinvict_directory="aux/other_tools/sinvict/"
sinvict="aux/other_tools/sinvict/sinvict"
# Clustering tools and scripts
calib="./calib"
starcode_umi="aux/other_tools/starcode/starcode-umi"
starcode="aux/other_tools/starcode/starcode"
convert_starcode="aux/convert_starcode_to_cluster.sh"
rainbow="aux/other_tools/rainbow/rainbow"
run_rainbow="aux/run_rainbow.sh"
convert_rainbow="aux/convert_rainbow_to_cluster.sh"
run_umitools="aux/run_umitools.sh"
convert_umitools="aux/convert_umitools_to_cluster.sh"
restore_cluster_missing_reads="aux/restore_cluster_missing_reads.sh"
# Plotting
run_plotting="slurm_scripts/real_tests_plotting.py"
verified_snv="experiments/real_tests/SNV.on_panel.verified.hg19.tsv"
reported_snv="experiments/real_tests/SNV.on_panel.reported.hg19.tsv"

./benchmark starcode rainbow umitools bwa calib $gtime
make -C $calib_cons_directory
make -C $bwa_directory
make -C $bam_readcount_directory
make -C $sinvict_directory

if [ ! -f $ref ]; then
    echo "Reference not found. Downloading..."
    wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz' -O chromFa.tar.gz
    tar xvzfO chromFa.tar.gz > $ref
    rm chromFa.tar.gz
fi
if [ ! -f "$ref.bwt" ]; then
    echo "Reference BWT index is not there. Indexing..."
    $bwa index $ref
fi

last_job_id=""
deps=""
mkdir -p $out_dir
for tool in umitools calib starcode rainbow raw
do
    # Slurm header
    echo "Preparing things for $tool"
    echo "#!/bin/bash"                                   > "$out_dir/$tool.pbs"
    echo "#SBATCH --job-name=$tool"                      >> "$out_dir/$tool.pbs"
    echo "#SBATCH -c 32"                                 >> "$out_dir/$tool.pbs"
    echo "#SBATCH --mem 512G"                            >> "$out_dir/$tool.pbs"
    echo "#SBATCH -t 23:59:59"                           >> "$out_dir/$tool.pbs"
    echo "#SBATCH --output=$out_dir/$tool.pbs.out"       >> "$out_dir/$tool.pbs"
    echo "#SBATCH --error=$out_dir/$tool.pbs.err"        >> "$out_dir/$tool.pbs"
    echo "#SBATCH --export=all"                          >> "$out_dir/$tool.pbs"
    echo "#SBATCH -p debug,express,normal,big-mem,long"  >> "$out_dir/$tool.pbs"
    echo "touch $out_dir/$tool.pbs.no_success"           >> "$out_dir/$tool.pbs"
    # Tool specific commands
    case $tool in
    "calib")
    echo "$gtime -o $out_dir/calib.gtime -v \\"  >> "$out_dir/$tool.pbs"
    echo "    $calib \\"                         >> "$out_dir/$tool.pbs"
    echo "        -f $fq_1 \\"                   >> "$out_dir/$tool.pbs"
    echo "        -r $fq_2 \\"                   >> "$out_dir/$tool.pbs"
    echo "        -o $out_dir/calib. \\"         >> "$out_dir/$tool.pbs"
    echo "        -l $barcode_length \\"         >> "$out_dir/$tool.pbs"
    echo "        --no-sort"                     >> "$out_dir/$tool.pbs"
    ;;
    "starcode")
    echo "$gtime -o $out_dir/starcode.gtime -v \\"   >> "$out_dir/$tool.pbs"
    echo "    $starcode_umi \\"                      >> "$out_dir/$tool.pbs"
    echo "        --starcode-path $starcode \\"      >> "$out_dir/$tool.pbs"
    echo "        --umi-len $barcode_length \\"      >> "$out_dir/$tool.pbs"
    echo "        --seq-trim 50 \\"                  >> "$out_dir/$tool.pbs"
    echo "        --umi-d 2 \\"                      >> "$out_dir/$tool.pbs"
    echo "        --seq-d 5 \\"                      >> "$out_dir/$tool.pbs"
    echo "        --umi-cluster-ratio 1 \\"          >> "$out_dir/$tool.pbs"
    echo "        --seq-cluster-ratio 1 \\"          >> "$out_dir/$tool.pbs"
    echo "        --seq-id \\"                       >> "$out_dir/$tool.pbs"
    echo "        $fq_1 \\"                          >> "$out_dir/$tool.pbs"
    echo "        $fq_2 \\"                          >> "$out_dir/$tool.pbs"
    echo "            > $out_dir/starcode.out"       >> "$out_dir/$tool.pbs"
    echo "$convert_starcode \\"                      >> "$out_dir/$tool.pbs"
    echo "    $out_dir/starcode.out \\"              >> "$out_dir/$tool.pbs"
    echo "        > $out_dir/starcode.cluster.temp"  >> "$out_dir/$tool.pbs"
    echo "$restore_cluster_missing_reads \\"         >> "$out_dir/$tool.pbs"
    echo "    $out_dir/starcode.cluster.temp \\"     >> "$out_dir/$tool.pbs"
    echo "    $fq_1 \\"                              >> "$out_dir/$tool.pbs"
    echo "        > $out_dir/starcode.cluster"       >> "$out_dir/$tool.pbs"
    ;;
    "rainbow")
    echo "$gtime -o $out_dir/rainbow.gtime -v \\"   >> "$out_dir/$tool.pbs"
    echo "    $run_rainbow \\"                      >> "$out_dir/$tool.pbs"
    echo "        $rainbow \\"                      >> "$out_dir/$tool.pbs"
    echo "        $fq_1 \\"                         >> "$out_dir/$tool.pbs"
    echo "        $fq_2 \\"                         >> "$out_dir/$tool.pbs"
    echo "        2 \\"                             >> "$out_dir/$tool.pbs"
    echo "        true \\"                          >> "$out_dir/$tool.pbs"
    echo "        $out_dir/rainbow.out"             >> "$out_dir/$tool.pbs"
    echo "$convert_rainbow \\"                      >> "$out_dir/$tool.pbs"
    echo "    $out_dir/rainbow.out \\"              >> "$out_dir/$tool.pbs"
    echo "        > $out_dir/rainbow.cluster.temp"  >> "$out_dir/$tool.pbs"
    echo "$restore_cluster_missing_reads \\"        >> "$out_dir/$tool.pbs"
    echo "    $out_dir/rainbow.cluster.temp \\"     >> "$out_dir/$tool.pbs"
    echo "    $fq_1 \\"                             >> "$out_dir/$tool.pbs"
    echo "        > $out_dir/rainbow.cluster"       >> "$out_dir/$tool.pbs"
    ;;
    "umitools")
    echo "$gtime -o $out_dir/umitools.gtime -v \\"   >> "$out_dir/$tool.pbs"
    echo "    $run_umitools \\"                      >> "$out_dir/$tool.pbs"
    echo "        $fq_1 \\"                          >> "$out_dir/$tool.pbs"
    echo "        $fq_2 \\"                          >> "$out_dir/$tool.pbs"
    echo "        $out_dir/umitools.out \\"          >> "$out_dir/$tool.pbs"
    echo "        $out_dir/umitools.work \\"         >> "$out_dir/$tool.pbs"
    echo "        $bwa \\"                           >> "$out_dir/$tool.pbs"
    echo "        $ref \\"                           >> "$out_dir/$tool.pbs"
    echo "        $samtools"                         >> "$out_dir/$tool.pbs"
    echo "$convert_umitools \\"                      >> "$out_dir/$tool.pbs"
    echo "    $out_dir/umitools.out \\"              >> "$out_dir/$tool.pbs"
    echo "        > $out_dir/umitools.cluster.temp"  >> "$out_dir/$tool.pbs"
    echo "$restore_cluster_missing_reads \\"         >> "$out_dir/$tool.pbs"
    echo "    $out_dir/umitools.cluster.temp \\"     >> "$out_dir/$tool.pbs"
    echo "    $fq_1 \\"                              >> "$out_dir/$tool.pbs"
    echo "        > $out_dir/umitools.cluster"       >> "$out_dir/$tool.pbs"
    ;;
    esac
    # Consensus, mapping, and SNV calling
    if [ $tool == "raw" ]
    then
        echo "cp $fq_1 $out_dir/$tool.1.fastq"  >> "$out_dir/$tool.pbs"
        echo "cp $fq_2 $out_dir/$tool.2.fastq"  >> "$out_dir/$tool.pbs"
        sinvict_depth=55
    else
        echo "$calib_cons \\"                    >> "$out_dir/$tool.pbs"
        echo "    -c $out_dir/$tool.cluster \\"  >> "$out_dir/$tool.pbs"
        echo "    -q \\"                         >> "$out_dir/$tool.pbs"
        echo "        $fq_1 \\"                  >> "$out_dir/$tool.pbs"
        echo "        $fq_2 \\"                  >> "$out_dir/$tool.pbs"
        echo "    -o \\"                         >> "$out_dir/$tool.pbs"
        echo "        $out_dir/$tool.1 \\"       >> "$out_dir/$tool.pbs"
        echo "        $out_dir/$tool.2 \\"       >> "$out_dir/$tool.pbs"
        echo "    -t 8"                          >> "$out_dir/$tool.pbs"
        sinvict_depth=3
    fi
    echo "$bwa mem \\"                         >> "$out_dir/$tool.pbs"
    echo "    -t 32 \\"                        >> "$out_dir/$tool.pbs"
    echo "    $ref \\"                         >> "$out_dir/$tool.pbs"
    echo "    $out_dir/$tool.1.fastq \\"       >> "$out_dir/$tool.pbs"
    echo "    $out_dir/$tool.2.fastq \\"       >> "$out_dir/$tool.pbs"
    echo "        | $samtools sort \\"         >> "$out_dir/$tool.pbs"
    echo "             - \\"                   >> "$out_dir/$tool.pbs"
    echo "             -O BAM \\"              >> "$out_dir/$tool.pbs"
    echo "             -@ 32 \\"               >> "$out_dir/$tool.pbs"
    echo "             -m 2G \\"               >> "$out_dir/$tool.pbs"
    echo "             -o $out_dir/$tool.bam"  >> "$out_dir/$tool.pbs"
    echo "$samtools index \\"      >> "$out_dir/$tool.pbs"
    echo "    $out_dir/$tool.bam"  >> "$out_dir/$tool.pbs"
    echo "mkdir -p $out_dir/$tool.bam-readcount"           >> "$out_dir/$tool.pbs"
    echo "$bam_readcount \\"                               >> "$out_dir/$tool.pbs"
    echo "    -l $panel \\"                                >> "$out_dir/$tool.pbs"
    echo "    -f $ref \\"                                  >> "$out_dir/$tool.pbs"
    echo "    -w 1 \\"                                     >> "$out_dir/$tool.pbs"
    echo "    $out_dir/$tool.bam \\"                       >> "$out_dir/$tool.pbs"
    echo "        > $out_dir/$tool.bam-readcount/out.tsv"  >> "$out_dir/$tool.pbs"
    for seq_err in 0.01 0.001 0.0001
    do
        echo "mkdir -p $out_dir/$tool.sinvict/$seq_err"  >> "$out_dir/$tool.pbs"
        echo "$sinvict \\"                               >> "$out_dir/$tool.pbs"
        echo "    -5 1 \\"                               >> "$out_dir/$tool.pbs"
        echo "    -e $seq_err \\"                        >> "$out_dir/$tool.pbs"
        echo "    -m $sinvict_depth \\"                  >> "$out_dir/$tool.pbs"
        echo "    -t $out_dir/$tool.bam-readcount \\"    >> "$out_dir/$tool.pbs"
        echo "    -o $out_dir/$tool.sinvict/$seq_err"    >> "$out_dir/$tool.pbs"
    done
    echo -e "rm $out_dir/$tool.pbs.no_success"  >> "$out_dir/$tool.pbs"
    last_job_id=$(sbatch "$out_dir/$tool.pbs")
    deps=$deps":"${last_job_id##* }
done

echo "Preparing things for plotting"
echo "#!/bin/bash"                                         >  "$out_dir/plotting.pbs"
echo "#SBATCH --job-name=real_tests_plotting"              >> "$out_dir/plotting.pbs"
echo "#SBATCH -c 1"                                        >> "$out_dir/plotting.pbs"
echo "#SBATCH --mem 16G"                                   >> "$out_dir/plotting.pbs"
echo "#SBATCH -t 02:59:59"                                 >> "$out_dir/plotting.pbs"
echo "#SBATCH --output=$out_dir/plotting.pbs.out"          >> "$out_dir/plotting.pbs"
echo "#SBATCH --error=$out_dir/plotting.pbs.err"           >> "$out_dir/plotting.pbs"
echo "#SBATCH --export=all"                                >> "$out_dir/plotting.pbs"
echo "#SBATCH -p debug,express,normal,big-mem,long"        >> "$out_dir/plotting.pbs"
echo "#SBATCH --dependency=afterany$deps"                  >> "$out_dir/plotting.pbs"
echo "$run_plotting $out_dir $verified_snv $reported_snv"  >> "$out_dir/plotting.pbs"
sbatch "$out_dir/plotting.pbs"
