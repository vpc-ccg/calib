#!/bin/bash

last_job_id=""
function slurm {
    filename=$1
    job_name=$2
    mem=$3
    tim=$4
    command=$5
    thread_count=$6
    dependencies=$7
    echo "#!/bin/bash"                                        >  $filename
    echo "#SBATCH --job-name=$job_name"                       >> $filename
    echo "#SBATCH -c $thread_count"                           >> $filename
    echo "#SBATCH --mem $mem"                                 >> $filename
    echo "#SBATCH -t $tim"                                    >> $filename
    echo "#SBATCH --output=$filename.out"                     >> $filename
    echo "#SBATCH --error=$filename.err"                      >> $filename
    echo "#SBATCH --export=all"                               >> $filename
    echo "#SBATCH -p debug,express,normal,big-mem,long"       >> $filename
    if [ $dependencies != "NONE" ]
    then
        echo "#SBATCH --dependency=afterany$dependencies"     >> $filename
    fi
    echo -e "$command"                                        >> $filename
    last_job_id=$(sbatch $filename)
}

random_seed=42
barcode_length=8
molecule_size_mu=300
read_length=150
sequencing_machine=HS25

./benchmark annotation reference reference_name=hg38 gnu_time cdhitest starcode rainbow
make
./benchmark panel reference_name=hg38 gene_list_name=COSMIC_cancer_genes random_seed=$random_seed

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
    num_molecules=25
    ;;
    esac
    slurm_path=slurm_pbs/"$dataset"
    mkdir -p "$slurm_path"


    # Making barcodes
    barcodes_deps=""
    job_name="barcodes"
    filename="$slurm_path/$job_name.pbs"
    mem="10240"
    tim="00:29:59"
    command="./benchmark barcodes "
    command=$command"num_barcodes=$num_barcodes "
    command=$command"barcode_length=$barcode_length "
    command=$command"random_seed=$random_seed "
    depends="NONE"
    slurm "$filename" "$job_name" "$mem" "$tim" "$command" "1" "$depends"
    barcodes_deps=$barcodes_deps":"${last_job_id##* }

    # Making molecules
    molecules_deps=""
    job_name="molecules"
    filename="$slurm_path/$job_name.pbs"
    mem="51200"
    tim="02:59:59"
    command="./benchmark molecules "
    command=$command"random_seed=$random_seed "
    command=$command"reference_name=hg38 "
    command=$command"gene_list_name=COSMIC_cancer_genes "
    command=$command"num_molecules=$num_molecules "
    command=$command"molecule_size_mu=$molecule_size_mu "
    command=$command"read_length=$read_length "
    depends="$barcodes_deps"
    slurm "$filename" "$job_name" "$mem" "$tim" "$command" "1" "$depends"
    molecules_deps=$molecules_deps":"${last_job_id##* }

    # Making simulation and log files
    simulate_deps=""
    job_name="simulate"
    filename="$slurm_path/$job_name.pbs"
    mem="51200"
    tim="05:59:59"
    command="./benchmark simulate make_log_files "
    command=$command"random_seed=$random_seed "
    command=$command"reference_name=hg38 "
    command=$command"gene_list_name=COSMIC_cancer_genes "
    command=$command"num_molecules=$num_molecules "
    command=$command"num_barcodes=$num_barcodes "
    command=$command"barcode_length=$barcode_length "
    command=$command"molecule_size_mu=$molecule_size_mu "
    command=$command"read_length=$read_length "
    command=$command"sequencing_machine=$sequencing_machine "
    depends="$barcodes_deps""$molecules_deps"
    slurm "$filename" "$job_name" "$mem" "$tim" "$command" "1" "$depends"
    simulate_deps=$simulate_deps":"${last_job_id##* }

    # calib_log
    job_name="calib"
    filename="$slurm_path/$job_name.pbs"
    mem="51200"
    tim="05:59:59"
    command="./benchmark calib_log "
    command=$command"log_comment=$dataset""_calib "
    command=$command"random_seed=$random_seed "
    command=$command"reference_name=hg38 "
    command=$command"gene_list_name=COSMIC_cancer_genes "
    command=$command"num_molecules=$num_molecules "
    command=$command"num_barcodes=$num_barcodes "
    command=$command"barcode_length=$barcode_length "
    command=$command"molecule_size_mu=$molecule_size_mu "
    command=$command"read_length=$read_length "
    command=$command"sequencing_machine=$sequencing_machine "
    depends="$simulate_deps"
    slurm "$filename" "$job_name" "$mem" "$tim" "$command" "1" "$depends"

    # starcode_log
    for starcode_seq_trim in 0 50 100
    do
        for starcode_umi_dist in 1 2 3
        do
            for starcode_umi_ratio in 1 3 5
            do
                for starcode_seq_dist in 3 5 8
                do
                    for starcode_seq_ratio in 1 3 5
                    do
                        job_name="starcode_"$starcode_umi_dist"_"$starcode_umi_ratio"_"$starcode_seq_dist"_"$starcode_seq_ratio"_"$starcode_seq_ratio""
                        filename="$slurm_path/$job_name.pbs"
                        mem="614400"
                        tim="05:59:59"
                        command="./benchmark starcode_log "
                        command=$command"log_comment=$dataset""_starcode "
                        command=$command"random_seed=$random_seed "
                        command=$command"reference_name=hg38 "
                        command=$command"gene_list_name=COSMIC_cancer_genes "
                        command=$command"num_molecules=$num_molecules "
                        command=$command"num_barcodes=$num_barcodes "
                        command=$command"barcode_length=$barcode_length "
                        command=$command"molecule_size_mu=$molecule_size_mu "
                        command=$command"read_length=$read_length "
                        command=$command"sequencing_machine=$sequencing_machine "
                        command=$command"starcode_seq_trim=$starcode_seq_trim "
                        command=$command"starcode_umi_dist=$starcode_umi_dist "
                        command=$command"starcode_umi_ratio=$starcode_umi_ratio "
                        command=$command"starcode_seq_dist=$starcode_seq_dist "
                        command=$command"starcode_seq_ratio=$starcode_seq_ratio "
                        depends="$simulate_deps"
                        slurm "$filename" "$job_name" "$mem" "$tim" "$command" "1" "$depends"
                    done
                done
            done
        done
    done

    # rainbow_log
    for rainbow_mismatch in 1 2 3 4 5 6 7 8 9
    do
        for rainbow_div in "true" "false"
        do
            # rainbow_log
            job_name="rainbow_"$rainbow_mismatch"_"$rainbow_div""
            filename="$slurm_path/$job_name.pbs"
            mem="51200"
            tim="05:59:59"
            command="./benchmark rainbow_log "
            command=$command"log_comment=$dataset""_rainbow "
            command=$command"random_seed=$random_seed "
            command=$command"reference_name=hg38 "
            command=$command"gene_list_name=COSMIC_cancer_genes "
            command=$command"num_molecules=$num_molecules "
            command=$command"num_barcodes=$num_barcodes "
            command=$command"barcode_length=$barcode_length "
            command=$command"molecule_size_mu=$molecule_size_mu "
            command=$command"read_length=$read_length "
            command=$command"sequencing_machine=$sequencing_machine "
            command=$command"rainbow_mismatch=$rainbow_mismatch "
            command=$command"rainbow_div=$rainbow_div "
            depends="$simulate_deps"
            slurm "$filename" "$job_name" "$mem" "$tim" "$command" "1" "$depends"
        done
    done

    # cd-hit-est
    for cdhitest_dist in 0.85 0.95 0.96 0.97 0.98
    do
        job_name="cdhitest_"$cdhitest_dist""
        filename="$slurm_path/$job_name.pbs"
        mem="204800"
        tim="05:59:59"
        command="./benchmark cdhitest_log "
        command=$command"log_comment=$dataset""_cdhit "
        command=$command"random_seed=$random_seed "
        command=$command"reference_name=hg38 "
        command=$command"gene_list_name=COSMIC_cancer_genes "
        command=$command"num_molecules=$num_molecules "
        command=$command"num_barcodes=$num_barcodes "
        command=$command"barcode_length=$barcode_length "
        command=$command"molecule_size_mu=$molecule_size_mu "
        command=$command"read_length=$read_length "
        command=$command"sequencing_machine=$sequencing_machine "
        command=$command"cdhitest_dist=$cdhitest_dist "
        depends="$simulate_deps"
        slurm "$filename" "$job_name" "$mem" "$tim" "$command" "1" "$depends"
    done
done
exit
