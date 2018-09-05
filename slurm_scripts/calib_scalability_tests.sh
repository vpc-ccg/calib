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
    echo -e "touch $filename.no_success"                      >> $filename
    echo -e "$command"                                        >> $filename
    echo -e "rm $filename.no_success"                         >> $filename
    last_job_id=$(sbatch $filename)
}

random_seed=1963
num_barcodes=5000
barcode_length=8
output_html="slurm_scripts/calib_scalability_tests_out.html"
slurm_path="slurm_scripts/slurm_files/calib_scalability_tests"
mkdir -p "$slurm_path"

./benchmark annotation reference reference_name=hg38 gnu_time calib
./benchmark barcodes num_barcodes=$num_barcodes random_seed=$random_seed barcode_length=$barcode_length
./benchmark panel gene_list_name=COSMIC_cancer_genes reference_name=hg38 random_seed=$random_seed

molecule_size_mu=300
read_length=150
sequencing_machine=HS25
calib_deps=""
for num_molecules in 100000 1000000 2000000
do
    echo "Molecule count $num_molecules"
    simulate_deps=""
    job_name="simulate.bl_$barcode_length.nm_$num_molecules.rl_$read_length"
    filename="$slurm_path/$job_name.pbs"
    mem="102400"
    tim="05:59:59"
    command="./benchmark simulate calib_log_file "
    command=$command"random_seed=$random_seed "
    command=$command"reference_name=hg38 "
    command=$command"gene_list_name=COSMIC_cancer_genes "
    command=$command"num_molecules=$num_molecules "
    command=$command"num_barcodes=$num_barcodes "
    command=$command"barcode_length=$barcode_length "
    command=$command"molecule_size_mu=$molecule_size_mu "
    command=$command"read_length=$read_length "
    command=$command"sequencing_machine=$sequencing_machine "
    thread_count=1
    depends="NONE"
    slurm "$filename" "$job_name" "$mem" "$tim" "$command" "$thread_count" "$depends"
    simulate_deps=$simulate_deps":"${last_job_id##* }

    for thread_count in 1 2 4 8
    do
        echo "Calib thread count $thread_count"
        job_name="calib.nm_$num_molecules.tc_$thread_count"
        filename="$slurm_path/$job_name.pbs"
        mem="102400"
        tim="05:59:59"
        command="./benchmark calib_log "
        command=$command"log_comment=t_$thread_count "
        command=$command"random_seed=$random_seed "
        command=$command"reference_name=hg38 "
        command=$command"gene_list_name=COSMIC_cancer_genes "
        command=$command"num_molecules=$num_molecules "
        command=$command"num_barcodes=$num_barcodes "
        command=$command"barcode_length=$barcode_length "
        command=$command"molecule_size_mu=$molecule_size_mu "
        command=$command"read_length=$read_length "
        command=$command"sequencing_machine=$sequencing_machine "
        command=$command"thread_count=$thread_count "
        depends="$simulate_deps"
        slurm "$filename" "$job_name" "$mem" "$tim" "$command" "$thread_count" "$depends"
        calib_deps=$calib_deps":"${last_job_id##* }
    done
done

job_name="scalibility_plots"
filename="$slurm_path/$job_name.pbs"
mem="10240"
tim="00:59:59"
command="slurm_scripts/calib_scalability_tests_plotting.py $output_html"
for num_molecules in 100000 1000000 2000000
do
    tsv_path="simulating/datasets"
    tsv_path=$tsv_path"/randS_"$random_seed
    tsv_path=$tsv_path"/barL_"$barcode_length
    tsv_path=$tsv_path".barNum_"$num_barcodes
    tsv_path=$tsv_path"/geneNum_35"
    tsv_path=$tsv_path".refName_hg38"
    tsv_path=$tsv_path".geneList_COSMIC_cancer_genes"
    tsv_path=$tsv_path"/ref_hg38.molMin_"$read_length
    tsv_path=$tsv_path".molMu_"$molecule_size_mu
    tsv_path=$tsv_path".molDev_25"
    tsv_path=$tsv_path".molNum"$num_molecules
    tsv_path=$tsv_path"/pcrC_7.pcrDR_0.6.pcrER_0.00005/seqMach_"$sequencing_machine
    tsv_path=$tsv_path".readL_"$read_length
    tsv_path=$tsv_path"/calib_benchmarks.tsv"
    command="$command $tsv_path"
done
depends="$calib_deps"
slurm "$filename" "$job_name" "$mem" "$tim" "$command" "$thread_count" "$depends"

exit
