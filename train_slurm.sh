#!/bin/bash

# slurm params
# 1 = filename
# 2 = job_name
# 3 = mem
# 4 = time
# 5 = command
# 6 = thread_count
function slurm {
    echo "#!/bin/bash"                                  >  $1
    echo "#SBATCH --job-name=$2"                        >> $1
    echo "#SBATCH -n $6"                                 >> $1
    echo "#SBATCH --mem $3"                             >> $1
    echo "#SBATCH -t $4"                                >> $1
    echo "#SBATCH --output=$1.out"                      >> $1
    echo "#SBATCH --error=$1.err"                       >> $1
    echo "#SBATCH --export=all"                         >> $1
    echo "#SBATCH -p debug,express,normal,big-mem,long" >> $1
    echo "$5"                                           >> $1
    sbatch $1
}

num_barcodes=5000
num_molecules=1000
thread_count=8
make
for barcode_length in {4,8,12};
do
    # Initializing the expirement by creating the log files and dataset
    for length_set in 75,HS20 150,HS25 250,MS;
    do
        IFS=",";
        set -- $length_set;
        read_length=$1;
        sequencing_machine=$2;
        molecule_size_mu=300
        ./benchmark simulate \
            reference_name=hg38 \
            bed=Panel.hg38 \
            num_molecules=$num_molecules \
            num_barcodes=$num_barcodes \
            barcode_length=$barcode_length \
            molecule_size_mu=$molecule_size_mu \
            read_length=$read_length \
            sequencing_machine=$sequencing_machine
        ./benchmark make_calib_log_file \
            reference_name=hg38 \
            bed=Panel.hg38 \
            num_molecules=$num_molecules \
            num_barcodes=$num_barcodes \
            barcode_length=$barcode_length \
            molecule_size_mu=$molecule_size_mu \
            read_length=$read_length \
            sequencing_machine=$sequencing_machine
    done
done
for barcode_length in {4,8,12};
do
    # Initializing the expirement by creating the log files and dataset
    for length_set in 75,HS20 150,HS25 250,MS;
    do
        IFS=",";
        set -- $length_set;
        read_length=$1;
        sequencing_machine=$2;
        molecule_size_mu=300

        slurm_path=slurm_pbs/
        mkdir -p "$slurm_path"
        # calib_log
        for error_tolerance in 1 2
        do
            for kmer_size in 4 8
            do
                for param_set in 3,1 3,2 4,1 4,2 4,3 5,2 5,3 5,4 6,2 6,3 6,4 6,5 7,2 7,3 7,4 7,5 7,6;
                do
                    IFS=",";
                    set -- $param_set;
                    minimizers_num=$1;
                    minimizers_threshold=$2;
                    job_name="calib.bl_$barcode_length.rl_$read_length.e_$error_tolerance.k_$kmer_size.m_$minimizers_num.t_$minimizers_threshold"
                    filename="$slurm_path/$job_name.pbs"
                    mem="51200"
                    tim="02:59:59"
                    command="./benchmark calib_log reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes barcode_length=$barcode_length molecule_size_mu=$molecule_size_mu read_length=$read_length sequencing_machine=$sequencing_machine barcode_error_tolerance=$error_tolerance kmer_size=$kmer_size minimizers_num=$minimizers_num minimizers_threshold=$minimizers_threshold thread_count=$thread_count"
                    slurm "$filename" "$job_name" "$mem" "$tim" "$command" "$thread_count"
                done
            done
        done
    done
done


exit
