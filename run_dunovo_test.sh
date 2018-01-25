#!/bin/bash

export PATH=$PATH:~/anaconda3/bin
./benchmark dunovo_conda

for dataset in 4
do
  case $dataset in
    1) echo "Test dataset one"
    num_barcodes=5000
    num_molecules=1000000
    ;;
    2)  echo "Test dataset two"
    num_barcodes=100
    num_molecules=100000
    ;;
    3)  echo "Test dataset three"
    num_barcodes=5000
    num_molecules=500000
    ;;
    4)  echo "Test dataset three"
    num_barcodes=100
    num_molecules=10
    ;;
  esac
  # bargoat_log
  for error_tolerance in 1
  do
    for kmer_size in 8
    do
        ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=$barcode_error_tolerance num_molecules=$num_molecules num_barcodes=$num_barcodes kmer_size=$kmer_size minimizers_num=5 minimizers_threshold=2
        ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=$barcode_error_tolerance num_molecules=$num_molecules num_barcodes=$num_barcodes kmer_size=$kmer_size minimizers_num=4 minimizers_threshold=1
    done
  done
  # du novo
  for dunovo_dist in 1 2 3 4
  do
    ./benchmark dunovo_log dunovo_prefix=~/bin/dunovo reference_name=hg38 bed=RCCPanelV2.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes dunovo_dist=$dunovo_dist
  done
  ./benchmark dunovo_conda_clean
done
exit

# ./benchmark bargoat_log rainbow=~/projects/bargoat/tools/rainbow_2.0.4/rainbow starcode=~/projects/bargoat/tools/starcode/starcode art_illumina=~/bin/art_bin_MountRainier/art_illumina reference_name=hg38 bed=RCCPanelV2.hg38 pcr_error_rate=0.00005 num_molecules=1000000 pcr_duplication_rate=0.6 num_barcodes=5000 barcode_length=8 minimizers_num=5 minimizers_threshold=3 kmer_size=4 rainbow_mismatch?=3
