#!/bin/bash

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
  # calib_log
  ./benchmark calib_log reference_name=hg38 bed=Panel.hg38 barcode_error_tolerance=1 kmer_size=8 minimizers_num=5 minimizers_threshold=2 num_molecules=$num_molecules num_barcodes=$num_barcodes
  ./benchmark calib_log reference_name=hg38 bed=Panel.hg38 barcode_error_tolerance=1 kmer_size=4 minimizers_num=5 minimizers_threshold=2 num_molecules=$num_molecules num_barcodes=$num_barcodes
  ./benchmark calib_log reference_name=hg38 bed=Panel.hg38 barcode_error_tolerance=2 kmer_size=8 minimizers_num=5 minimizers_threshold=2 num_molecules=$num_molecules num_barcodes=$num_barcodes
  ./benchmark calib_log reference_name=hg38 bed=Panel.hg38 barcode_error_tolerance=2 kmer_size=4 minimizers_num=5 minimizers_threshold=2 num_molecules=$num_molecules num_barcodes=$num_barcodes

  # starcode_log
  for starcode_dist in 1 2 3 4
  do
    ./benchmark starcode_log reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes starcode_dist=$starcode_dist
  done
  # rainbow_log
  for rainbow_mismatch in 1 2 3 4 5 6 7 8 9
  do
    ./benchmark rainbow_log reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes rainbow_mismatch=$rainbow_mismatch
    ./benchmark rainbow_log reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes rainbow_mismatch=$rainbow_mismatch rainbow_div=true
  done

done
exit

# ./benchmark calib_log rainbow=~/projects/calib/tools/rainbow_2.0.4/rainbow starcode=~/projects/calib/tools/starcode/starcode art_illumina=~/bin/art_bin_MountRainier/art_illumina reference_name=hg38 bed=RCCPanelV2.hg38 pcr_error_rate=0.00005 num_molecules=1000000 pcr_duplication_rate=0.6 num_barcodes=5000 barcode_length=8 minimizers_num=5 minimizers_threshold=3 kmer_size=4 rainbow_mismatch?=3
