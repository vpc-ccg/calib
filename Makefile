#!/usr/bin/make -rRf
# Usage:
# calib correct reads1=reads1.fq reads2=reads2.fq
# calib simulate
#

# Figuring out where we are
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir :=  $(patsubst %/,%,$(dir $(mkfile_path)))/
clustering_path?=$(current_dir)clustering/
simulating_path?=$(current_dir)simulating/

SHELL=/bin/bash -o pipefail

cc?=g++
python3?= PYTHONHASHSEED=0 python3
cc_files?= $(clustering_path)*.cc $(clustering_path)*.h
cc_flags?=
cc_args?= $(cc_flags) -O3 -std=c++11
art_illumina?=art_illumina


# Common arguments
random_seed?=42

# Simulation common
simulation_datasets_path?=$(simulating_path)datasets/
bed?=
bed_flag=--bed $(references_path)$(bed).bed
bed_prefix=_bed$(bed)
ifeq ($(bed),)
	bed_flag=
	bed_prefix=NA
endif
reference_name?=e_coli
simulation_prefix?=$(simulation_datasets_path)rs_$(random_seed).

references_path?=$(simulating_path)genomes/

## Generating barcodes
num_barcodes?=100
barcode_length?=8
barcodes_params?=bl_$(barcode_length).nb_$(num_barcodes).
simulated_barcodes?=$(simulation_prefix)$(barcodes_params)barcodes.txt

## Generating molecules
molecule_size_mu?=200
molecule_size_dev?=25
num_molecules?=500
read_length?=100
reference?=$(references_path)$(reference_name).fa
molecules_params?=ref_$(reference_name).bed_$(bed_prefix).mu_$(molecule_size_mu).dev_$(molecule_size_dev).nm$(num_molecules).
simulated_molecules?=$(simulation_prefix)$(molecules_params)molecules.fa


## Barcoding molecules
simulated_barcoded_molecules?=$(simulation_prefix)$(barcodes_params)$(molecules_params)barcoded_molecules.fa

## PCR duplication
pcr_cycles?=7
pcr_duplication_rate?=0.6
pcr_error_rate?=0.00005
pcr_params?=pc_$(pcr_cycles).pdr_$(pcr_duplication_rate).per_$(pcr_error_rate).
amplified_barcoded_molecules?=$(simulation_prefix)$(barcodes_params)$(molecules_params)$(pcr_params)amplified_barcoded_molecules.fa

## Generating reads
sequencing_machine?=HS20
simulated_reads?=$(simulation_prefix)$(barcodes_params)$(molecules_params)$(pcr_params)sm_$(sequencing_machine).
simulated_reads_log?=$(simulated_reads)art_illumina.log

# Clustering arguments
input_reads_prefix?=$(simulated_reads)
forward_reads?=$(input_reads_prefix)1.fq
reverse_reads?=$(input_reads_prefix)2.fq
ignored_sequence_prefix_length?=0
minimizers_num?=3
kmer_size?=8
barcode_error_tolerance?=2
minimizers_threshold?=1
silent?=--silent

clustering_params?=bl_$(barcode_length).mn_$(minimizers_num).ks_$(kmer_size).bet_$(barcode_error_tolerance).mt_$(minimizers_threshold).
output_prefix?=$(input_reads_prefix)calib.$(clustering_params)

# Rand Index Accuracy arguments
cluster_file?=$(output_prefix)cluster
input_amplified_molecules?=$(amplified_barcoded_molecules)
output_accuracy_results?=$(output_prefix)accuracy


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: help clean simulate_clean simulate cluster accuracy
calib: $(cc_files)
	$(cc) $(cc_files) $(cc_args) -o $(current_dir)calib

help:
	@echo 'calib: Clustering without alignment using LSH and MinHashing of barcoded reads'
	@echo 'Usage: make [COMMAND]... [PARAMETER=VALUE]...'
	@echo 'Example: make simulate'
	@echo '			make cluster forward_reads=R1.fastq reverse_reads=R2.fastq barcode_length=8'
	@echo 'Checking the possible list in Makefile file. The main ones for cluster command are:'
	@echo '        forward_reads (string)'
	@echo '        reverse_reads (string)'
	@echo '        output_prefix (string)'
	@echo '        barcode_length (int)'
	@echo '        ignored_sequence_prefix_length (int)'
	@echo '        minimizers_num (int)'
	@echo '        kmer_size (int)'
	@echo '        barcode_error_tolerance (int < barcode length)'
	@echo '        minimizers_threshold (int < number of minimizers)'
	@echo '        silent (no value)'
	@echo 'Note that you can run Calib clustering module from the executable in CURRENT_DIR/calib:'
	@echo 'The main ones for simulate command are:'
	@echo '        random_seed (int)'
	@echo '        bed (bed file should be in simulating/genomes/<bed>.bed)'
	@echo '        reference_name (fasta file should be in simulating/genomes/<reference_name>.fa)'
	@echo '        num_barcodes (int)'
	@echo '        barcode_length (int)'
	@echo '        molecule_size_mu (int)'
	@echo '        molecule_size_dev (int)'
	@echo '        num_molecules (int)'
	@echo '        read_length (int)'
	@echo '        pcr_cycles (int)'
	@echo '        pcr_duplication_rate (float (0,1))'
	@echo '        pcr_error_rate (float (0,1))'
	@echo '        sequencing_machine (string from (HS10, HS20, HS25, HSXn, HSXt, MinS, MSv1, MSv3, NS50))'
	@echo 'You can pipeline the whole process in a single run:'
	@echo '        make simulate cluster accuracy'
	@echo 'The results will be stored in simulating/datasets/* with file names for every part of this pipeline'


$(simulated_barcodes):
	$(python3) $(simulating_path)generate_barcodes.py \
		--num-of-barcodes $(num_barcodes) \
		--len-of-one-end-barcode $(barcode_length) \
		--random-seed $(random_seed) \
		--output-barcodes $(simulated_barcodes)

$(simulated_molecules):
	$(python3) $(simulating_path)generate_molecules.py \
		--reference $(reference) \
		--number-of-molecules $(num_molecules) \
		--molecule-size-mean $(molecule_size_mu) \
		--molecule-size-standard-dev $(molecule_size_dev) \
		--min-molecule-size $(read_length) \
		--random-seed $(random_seed) \
		--output-molecules $(simulated_molecules) \
		$(bed_flag)

$(simulated_barcoded_molecules): $(simulated_barcodes) $(simulated_molecules)
	$(python3) $(simulating_path)attach_barcodes_to_molecules.py \
		--input-barcodes $(simulated_barcodes) \
		--input-molecules $(simulated_molecules) \
		--random-seed $(random_seed) \
		--output-barcoded-molecules $(simulated_barcoded_molecules)

$(amplified_barcoded_molecules): $(simulated_barcoded_molecules)
	$(python3) $(simulating_path)pcr_duplication.py \
		--molecules $(simulated_barcoded_molecules) \
		--number-of-cycles $(pcr_cycles) \
		--duplication-rate-per-cycle $(pcr_duplication_rate) \
		--error-rate $(pcr_error_rate) \
		--random-seed $(random_seed) \
		--pcr-product $(amplified_barcoded_molecules)

$(simulated_reads_log): $(amplified_barcoded_molecules)
	$(art_illumina) \
		--seqSys $(sequencing_machine) \
		--amplicon \
		--paired \
		--noALN \
		--in $(amplified_barcoded_molecules) \
		--len $(read_length) \
		--fcov 1 \
		--rndSeed $(random_seed) \
		--out $(simulated_reads) \
		> $(simulated_reads_log);


$(reverse_reads): $(simulated_reads_log)
$(forward_reads): $(reverse_reads)
simulate: $(forward_reads) $(reverse_reads) $(simulated_reads_log)

simulate_clean:
	rm -f $(simulation_datasets_path)*



cluster: calib $(forward_reads) $(reverse_reads)
	$(current_dir)calib \
		--input-forward $(forward_reads) \
		--input-reverse $(reverse_reads) \
		--output-prefix $(output_prefix) \
		--barcode-length $(barcode_length) \
		--ignored-sequence-prefix-length $(ignored_sequence_prefix_length) \
		--minimizer-count $(minimizers_num) \
		--kmer-size $(kmer_size) \
		--error-tolerance $(barcode_error_tolerance) \
		--minimizer-threshold $(minimizers_threshold) \
		$(silent)

accuracy: $(cluster_file)
	$(python3) $(simulating_path)rand_index.py \
		--input-cluster-file $(cluster_file) \
		--input-amplified-molecule $(input_amplified_molecules) \
		--output-accuracy-results $(output_accuracy_results)

benchmark: $(clustering_path)calib.o $(forward_reads) $(reverse_reads)

clean:
	rm -f $(current_dir)calib
