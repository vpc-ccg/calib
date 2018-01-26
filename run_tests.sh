#!/bin/bash

for dataset in $@
do
  case $dataset in
    "small") echo "Test dataset one"
    num_barcodes=5000
    num_molecules=1000000
    ;;
    "medium")  echo "Test dataset two"
    num_barcodes=100
    num_molecules=100000
    ;;
    "large")  echo "Test dataset three"
    num_barcodes=5000
    num_molecules=500000
    ;;
    "tiny")  echo "Test dataset three"
    num_barcodes=100
    num_molecules=10
    ;;
  esac
  # bargoat_log
  ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=1 kmer_size=8 minimizers_num=5 minimizers_threshold=2 num_molecules=$num_molecules num_barcodes=$num_barcodes
  ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=1 kmer_size=4 minimizers_num=5 minimizers_threshold=2 num_molecules=$num_molecules num_barcodes=$num_barcodes
  ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=1 kmer_size=8 minimizers_num=4 minimizers_threshold=1 num_molecules=$num_molecules num_barcodes=$num_barcodes
  ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=2 kmer_size=8 minimizers_num=5 minimizers_threshold=2 num_molecules=$num_molecules num_barcodes=$num_barcodes
  ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=2 kmer_size=8 minimizers_num=4 minimizers_threshold=1 num_molecules=$num_molecules num_barcodes=$num_barcodes
  ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=2 kmer_size=8 minimizers_num=3 minimizers_threshold=1 num_molecules=$num_molecules num_barcodes=$num_barcodes

  # for barcode_error_tolerance in 1 2
  # do
  #   for kmer_size in 4 8
  #   do
  #     ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=$barcode_error_tolerance num_molecules=$num_molecules num_barcodes=$num_barcodes kmer_size=$kmer_size minimizers_num=3 minimizers_threshold=1
  #     ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=$barcode_error_tolerance num_molecules=$num_molecules num_barcodes=$num_barcodes kmer_size=$kmer_size minimizers_num=3 minimizers_threshold=2
  #     ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=$barcode_error_tolerance num_molecules=$num_molecules num_barcodes=$num_barcodes kmer_size=$kmer_size minimizers_num=4 minimizers_threshold=1
  #     ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=$barcode_error_tolerance num_molecules=$num_molecules num_barcodes=$num_barcodes kmer_size=$kmer_size minimizers_num=4 minimizers_threshold=2
  #     ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=$barcode_error_tolerance num_molecules=$num_molecules num_barcodes=$num_barcodes kmer_size=$kmer_size minimizers_num=4 minimizers_threshold=3
  #     ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=$barcode_error_tolerance num_molecules=$num_molecules num_barcodes=$num_barcodes kmer_size=$kmer_size minimizers_num=5 minimizers_threshold=2
  #     ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=$barcode_error_tolerance num_molecules=$num_molecules num_barcodes=$num_barcodes kmer_size=$kmer_size minimizers_num=5 minimizers_threshold=3
  #     ./benchmark bargoat_log reference_name=hg38 bed=RCCPanelV2.hg38 barcode_error_tolerance=$barcode_error_tolerance num_molecules=$num_molecules num_barcodes=$num_barcodes kmer_size=$kmer_size minimizers_num=5 minimizers_threshold=4
  #   done
  # done

  # starcode_log
  for starcode_dist in 1 2 3 4
  do
    ./benchmark starcode_log reference_name=hg38 bed=RCCPanelV2.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes starcode_dist=$starcode_dist
  done
  # rainbow_log
  for rainbow_mismatch in 1 2 3 4
  do
    ./benchmark rainbow_log reference_name=hg38 bed=RCCPanelV2.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes rainbow_mismatch=$rainbow_mismatch
    ./benchmark rainbow_log reference_name=hg38 bed=RCCPanelV2.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes rainbow_mismatch=$rainbow_mismatch rainbow_div=true
  done

done
exit

# ./benchmark bargoat_log rainbow=~/projects/bargoat/tools/rainbow_2.0.4/rainbow starcode=~/projects/bargoat/tools/starcode/starcode art_illumina=~/bin/art_bin_MountRainier/art_illumina reference_name=hg38 bed=RCCPanelV2.hg38 pcr_error_rate=0.00005 num_molecules=1000000 pcr_duplication_rate=0.6 num_barcodes=5000 barcode_length=8 minimizers_num=5 minimizers_threshold=3 kmer_size=4 rainbow_mismatch?=3
