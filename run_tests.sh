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
    # Initializing the expirement by creating the log files and dataset
    ./benchmark make_log_files simulate reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes

    # calib_log
    for error_tolerance in 1 2
    do
        for kmer_size in 4 8
        do
            ./benchmark calib_log reference_name=hg38 bed=Panel.hg38 barcode_error_tolerance=1 kmer_size=$kmer_size minimizers_num=3 minimizers_threshold=1 num_molecules=$num_molecules num_barcodes=$num_barcodes
            ./benchmark calib_log reference_name=hg38 bed=Panel.hg38 barcode_error_tolerance=1 kmer_size=$kmer_size minimizers_num=3 minimizers_threshold=2 num_molecules=$num_molecules num_barcodes=$num_barcodes
            ./benchmark calib_log reference_name=hg38 bed=Panel.hg38 barcode_error_tolerance=2 kmer_size=$kmer_size minimizers_num=4 minimizers_threshold=1 num_molecules=$num_molecules num_barcodes=$num_barcodes
            ./benchmark calib_log reference_name=hg38 bed=Panel.hg38 barcode_error_tolerance=2 kmer_size=$kmer_size minimizers_num=4 minimizers_threshold=2 num_molecules=$num_molecules num_barcodes=$num_barcodes
            ./benchmark calib_log reference_name=hg38 bed=Panel.hg38 barcode_error_tolerance=1 kmer_size=$kmer_size minimizers_num=4 minimizers_threshold=3 num_molecules=$num_molecules num_barcodes=$num_barcodes
            ./benchmark calib_log reference_name=hg38 bed=Panel.hg38 barcode_error_tolerance=1 kmer_size=$kmer_size minimizers_num=5 minimizers_threshold=2 num_molecules=$num_molecules num_barcodes=$num_barcodes
            ./benchmark calib_log reference_name=hg38 bed=Panel.hg38 barcode_error_tolerance=2 kmer_size=$kmer_size minimizers_num=5 minimizers_threshold=3 num_molecules=$num_molecules num_barcodes=$num_barcodes
            ./benchmark calib_log reference_name=hg38 bed=Panel.hg38 barcode_error_tolerance=2 kmer_size=$kmer_size minimizers_num=5 minimizers_threshold=4 num_molecules=$num_molecules num_barcodes=$num_barcodes
        done
    done

    # starcode_log
    for starcode_dist in 4 3 2 1
    do
      for starcode_ratio in 1 2 3 4 5
      do
          ./benchmark starcode_log reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes starcode_dist=$starcode_dist starcode_ratio=$starcode_ratio
      done
    done

    # rainbow_log
    for rainbow_mismatch in 1 2 3 4 5 6 7 8 9
    do
      ./benchmark rainbow_log reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes rainbow_mismatch=$rainbow_mismatch
      ./benchmark rainbow_log reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes rainbow_mismatch=$rainbow_mismatch rainbow_div=true
    done

    # cd-hit-est
    for cdhitest_dist in 0.95 0.90 0.85
    do
    ./benchmark cdhitest_log reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes cdhitest_dist=$cdhitest_dist
    done

    # du novo
    for dunovo_dist in 1 2 3 4
    do
    ./benchmark dunovo_log reference_name=hg38 bed=Panel.hg38 num_molecules=$num_molecules num_barcodes=$num_barcodes dunovo_dist=$dunovo_dist
    done

done
exit
