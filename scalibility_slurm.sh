#!/bin/bash


echo ./benchmark annotation reference reference_name=hg38 gnu_time
echo make

random_seed=5000
num_barcodes=5000
barcode_length=8
# echo ./benchmark barcodes num_barcodes=$num_barcodes random_seed=$random_seed barcode_length=$barcode_length
# echo ./benchmark panel gene_list_name=COSMIC_cancer_genes reference_name=hg38 random_seed=$random_seed
molecule_size_mu=300
read_length=150
sequencing_machine=HS25
for num_molecules in 100000 1000000 2000000
do
    command="./benchmark simulate make_calib_log_file "
    command=$command"random_seed=$random_seed "
    command=$command"reference_name=hg38 "
    command=$command"gene_list_name=COSMIC_cancer_genes "
    command=$command"num_molecules=$num_molecules "
    command=$command"num_barcodes=$num_barcodes "
    command=$command"barcode_length=$barcode_length "
    command=$command"molecule_size_mu=$molecule_size_mu "
    command=$command"read_length=$read_length "
    command=$command"sequencing_machine=$sequencing_machine "
    # echo -e $command
done
error_tolerance=2
kmer_size=8
minimizers_num=5
minimizers_threshold=2
for num_molecules in 100000 1000000 2000000
do
    for thread_count in 1 4 8
    do
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
        command=$command"barcode_error_tolerance=$error_tolerance "
        command=$command"kmer_size=$kmer_size "
        command=$command"minimizers_num=$minimizers_num "
        command=$command"minimizers_threshold=$minimizers_threshold "
        command=$command"thread_count=$thread_count "
        # echo -e $command
    done
done

random_seed=6000
num_barcodes=5000
barcode_length=8
echo ./benchmark barcodes num_barcodes=$num_barcodes random_seed=$random_seed barcode_length=$barcode_length
echo ./benchmark panel gene_list_name=COSMIC_cancer_genes reference_name=hg38 random_seed=$random_seed
molecule_size_mu=300
read_length=150
sequencing_machine=HS25
for num_molecules in 100000 1000000 2000000
do
    command="./benchmark molecules "
    command=$command"random_seed=$random_seed "
    command=$command"reference_name=hg38 "
    command=$command"gene_list_name=COSMIC_cancer_genes "
    command=$command"num_molecules=$num_molecules "
    command=$command"molecule_size_mu=$molecule_size_mu "
    command=$command"read_length=$read_length "
    echo -e $command
done
for num_molecules in 100000 1000000 2000000
do
    for pcr_cycles in 5 7 9
    do
        command="./benchmark simulate make_calib_log_file "
        command=$command"random_seed=$random_seed "
        command=$command"reference_name=hg38 "
        command=$command"gene_list_name=COSMIC_cancer_genes "
        command=$command"num_molecules=$num_molecules "
        command=$command"num_barcodes=$num_barcodes "
        command=$command"barcode_length=$barcode_length "
        command=$command"molecule_size_mu=$molecule_size_mu "
        command=$command"pcr_cycles=$pcr_cycles "
        command=$command"read_length=$read_length "
        command=$command"sequencing_machine=$sequencing_machine "
        echo -e $command
    done
done

error_tolerance=2
kmer_size=8
minimizers_num=5
minimizers_threshold=2
thread_count=1
for num_molecules in 100000 1000000 2000000
do
    for pcr_cycles in 5 7 9
    do
        command="./benchmark calib_log "
        command=$command"log_comment=c_$pcr_cycles "
        command=$command"random_seed=$random_seed "
        command=$command"reference_name=hg38 "
        command=$command"gene_list_name=COSMIC_cancer_genes "
        command=$command"num_molecules=$num_molecules "
        command=$command"num_barcodes=$num_barcodes "
        command=$command"barcode_length=$barcode_length "
        command=$command"molecule_size_mu=$molecule_size_mu "
        command=$command"pcr_cycles=$pcr_cycles "
        command=$command"read_length=$read_length "
        command=$command"sequencing_machine=$sequencing_machine "
        command=$command"barcode_error_tolerance=$error_tolerance "
        command=$command"kmer_size=$kmer_size "
        command=$command"minimizers_num=$minimizers_num "
        command=$command"minimizers_threshold=$minimizers_threshold "
        command=$command"thread_count=$thread_count "
        # echo -e $command
    done
done

exit
