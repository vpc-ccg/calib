#!/bin/bash

# export PATH=$PATH:~/anaconda3/bin
dunovo_path=~/bin/dunovo
./benchmark dunovo_conda

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
    "tiny")  echo "Test dataset Four"
    num_barcodes=100
    num_molecules=10
    ;;
  esac
  # du novo
  for dunovo_dist in 1 2 3 4
  do
    ./benchmark dunovo_log dunovo_prefix="$dunovo_path" reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes dunovo_dist=$dunovo_dist
  done
  ./benchmark dunovo_conda_clean
done
exit

# ./benchmark calib_log rainbow=~/projects/calib/tools/rainbow_2.0.4/rainbow starcode=~/projects/calib/tools/starcode/starcode art_illumina=~/bin/art_bin_MountRainier/art_illumina reference_name=hg38 bed=RCCPanelV2.hg38 pcr_error_rate=0.00005 num_molecules=1000000 pcr_duplication_rate=0.6 num_barcodes=5000 barcode_length=8 minimizers_num=5 minimizers_threshold=3 kmer_size=4 rainbow_mismatch?=3
