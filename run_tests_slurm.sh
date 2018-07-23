#!/bin/bash

# slurm params
# 1 = filename
# 2 = job_name
# 3 = mem
# 4 = time
# 5 = command
function slurm {
    echo "#!/bin/bash"                                  >  $1
    echo "#SBATCH --job-name=$2"                        >> $1
    echo "#SBATCH -n 1"                                 >> $1
    echo "#SBATCH --mem $3"                             >> $1
    echo "#SBATCH -t $4"                                >> $1
    echo "#SBATCH --output=$1.out"                      >> $1
    echo "#SBATCH --error=$1.err"                       >> $1
    echo "#SBATCH --export=all"                         >> $1
    echo "#SBATCH -p debug,express,normal,big-mem,long" >> $1
    echo "$5"                                           >> $1
    sbatch $1
}
for dataset in $@
do
    case $dataset in
    "small") echo "Test dataset one"
    num_barcodes=100
    num_molecules=100000
    ;;
    "medium")  echo "Test dataset two"
    num_barcodes=5000
    num_molecules=1000000
    ;;
    "large")  echo "Test dataset three"
    num_barcodes=25000
    num_molecules=1000000
    ;;
    "giant")  echo "Test dataset four"
    num_barcodes=25000
    num_molecules=2000000
    ;;
    "tiny")  echo "Test dataset five"
    num_barcodes=100
    num_molecules=10
    ;;
    esac
    # Initializing the expirement by creating the log files and dataset
    echo Initializing the expirement by creating the log files and dataset
    ./benchmark simulate reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes
    ./benchmark make_log_files reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes

    slurm_path=slurm_pbs/"$dataset"
    mkdir -p "$slurm_path"

   # calib_log
    for error_tolerance in 1 2
    do
        for kmer_size in 4 8
        do
            for param_set in 3,1 3,2 4,1 4,2 4,3 5,2 5,3 5,4;
            do
                IFS=",";
                set -- $param_set;
                job_name="calib.$error_tolerance.$kmer_size.$1.$2"
                filename="$slurm_path/$job_name.pbs"
                mem="102400"
                tim="05:59:59"
                command="./benchmark calib_log reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes barcode_error_tolerance=$error_tolerance kmer_size=$kmer_size minimizers_num=$1 minimizers_threshold=$2"
                slurm "$filename" "$job_name" "$mem" "$tim" "$command"
            done
        done
    done

    # starcode_log
    for starcode_dist in 1 2 3 4 5 6 7 8
    do
      for starcode_ratio in 1 2 3 4 5
      do
          job_name="starcode.$starcode_dist.$starcode_ratio"
          filename="$slurm_path/$job_name.pbs"
          mem="614400"
          tim="11:59:59"
          command="./benchmark starcode_log reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes starcode_dist=$starcode_dist starcode_ratio=$starcode_ratio"
          slurm "$filename" "$job_name" "$mem" "$tim" "$command"
      done
    done

    # rainbow_log
    for rainbow_mismatch in 1 2 3 4 5 6 7 8 9
    do
        for div in "true" "false"
        do
            job_name="rainbow.$rainbow_mismatch.$div"
            filename="$slurm_path/$job_name.pbs"
            mem="51200"
            tim="05:59:59"
            command="./benchmark rainbow_log reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes rainbow_mismatch=$rainbow_mismatch rainbow_div=$div"
            slurm "$filename" "$job_name" "$mem" "$tim" "$command"
        done
    done

    # cd-hit-est
    for cdhitest_dist in 0.85 0.95 0.96 0.97 0.98
    do
        job_name="cdhitest.$cdhitest_dist"
        filename="$slurm_path/$job_name.pbs"
        mem="204800"
        tim="23:59:59"
        command="./benchmark cdhitest_log reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes cdhitest_dist=$cdhitest_dist"
        slurm "$filename" "$job_name" "$mem" "$tim" "$command"
    done

done
exit
