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

CXX?=g++
CXXFLAGS=-std=c++11 -pthread -I clustering/ -lz -Wall -Wextra
SOURCES:=$(wildcard clustering/*.cc)
OBJECTS:=$(SOURCES:.cc=.o)
EXECUTABLE=calib
DBG_OBJECTS:=$(SOURCES:.cc=.dbg)
DEBUGGABLE=$(EXECUTABLE)_dbg

python3?= PYTHONHASHSEED=0 python3
art_illumina?=art_illumina
CONVERT_FASTQ=$(simulating_path)convert_fastq_to_true_cluster.sh

# Simulation common
random_seed?=42
simulation_prefix?=$(simulation_datasets_path)randS_$(random_seed)/
simulation_datasets_path?=$(simulating_path)datasets/
references_path?=$(simulating_path)genomes/
reference_name?=hg38

## Generating panel
annotation?=$(references_path)$(reference_name).gtf
gene_list_name?=
gene_list?=$(references_path)$(reference_name).$(gene_list_name).txt
num_genes?=35
panel_name?=geneNum_$(num_genes).refName_$(reference_name).geneList_$(gene_list_name)
ifeq ($(gene_list_name),)
	panel_name=noPanel
endif
panel_params?=$(panel_name)/
panel_prefix?=$(simulation_prefix)$(panel_params)
panel?=$(panel_prefix)panel.bed


## Generating molecules
bed_flag=--bed $(panel)
ifeq ($(panel_name),noPanel)
	bed_flag=
endif
molecule_size_mu?=200
molecule_size_dev?=25
min_molecule_size?=$(read_length)
num_molecules?=500
read_length?=100
reference?=$(references_path)$(reference_name).fa
molecules_params?=ref_$(reference_name).molMin_$(min_molecule_size).molMu_$(molecule_size_mu).molDev_$(molecule_size_dev).molNum$(num_molecules)/
molecules_prefix?=$(simulation_prefix)$(panel_params)$(molecules_params)
molecules?=$(molecules_prefix)molecules.fa

## Generating barcodes
num_barcodes?=100
barcode_length?=8
barcodes_params?=barL_$(barcode_length).barNum_$(num_barcodes)/
barcodes_prefix?=$(simulation_prefix)$(barcodes_params)
barcodes?=$(barcodes_prefix)barcodes.txt


## Barcoding molecules
barcoded_molecules_prefix?=$(simulation_prefix)$(barcodes_params)$(panel_params)$(molecules_params)
barcoded_molecules?=$(barcoded_molecules_prefix)barcoded_molecules.fa

## PCR duplication
pcr_cycles?=7
pcr_duplication_rate?=0.6
pcr_error_rate?=0.00005
pcr_params?=pcrC_$(pcr_cycles).pcrDR_$(pcr_duplication_rate).pcrER_$(pcr_error_rate)/
amplified_barcoded_molecules_prefix?=$(simulation_prefix)$(barcodes_params)$(panel_params)$(molecules_params)$(pcr_params)
amplified_barcoded_molecules?=$(amplified_barcoded_molecules_prefix)amplified_barcoded_molecules.fa

## Generating reads
sequencing_machine?=HS20
reads_prefix?=$(simulation_prefix)$(barcodes_params)$(panel_params)$(molecules_params)$(pcr_params)seqMach_$(sequencing_machine).readL_$(read_length)/
true_cluster?=$(reads_prefix)true.cluster
reads_log?=$(reads_prefix)art_illumina.log

# Clustering arguments
input_reads_prefix?=$(reads_prefix)
forward_reads?=$(input_reads_prefix)1.fq
reverse_reads?=$(input_reads_prefix)2.fq
ignored_sequence_prefix_length?=0
thread_count?=1
minimizers_num?=-1
kmer_size?=-1
barcode_error_tolerance?=-1
minimizers_threshold?=-1
no_sort?=--no-sort
silent?=--silent

calib_params?=l_$(barcode_length).m_$(minimizers_num).k_$(kmer_size).e_$(barcode_error_tolerance).m_$(minimizers_threshold).t_$(thread_count).
calib_output_prefix?=$(input_reads_prefix)calib.$(calib_params)

# Rand Index Accuracy arguments
cluster_file?=$(calib_output_prefix)cluster
input_amplified_molecules?=$(amplified_barcoded_molecules)
output_accuracy_results?=$(calib_output_prefix)accuracy


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: help clean simulate_clean simulate cluster accuracy
$(EXECUTABLE): $(EXECUTABLE).cc $(OBJECTS) Makefile
	$(CXX) $(OBJECTS) -O2 $< $(CXXFLAGS) -o $@

$(DEBUGGABLE): $(EXECUTABLE).cc $(DBG_OBJECTS) Makefile
	$(CXX) $(CXXFLAGS) $(DBG_OBJECTS) -g -O0 $< -o $@ -DGITINFO="\"$(GITINFO)\""

%.o: %.cc Makefile
	$(CXX) $(CXXFLAGS) -O2 -c $< -o $@

%.dbg: %.cc Makefile
	$(CXX) $(CXXFLAGS) -g -O0 -c $< -o $@

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

$(barcodes):
	@echo "Simulating barcodes"
	mkdir -p $(barcodes_prefix)
	$(python3) $(simulating_path)generate_barcodes.py \
		--num-of-barcodes $(num_barcodes) \
		--len-of-one-end-barcode $(barcode_length) \
		--random-seed $(random_seed) \
		--output-barcodes $(barcodes);
	chmod -w $(barcodes);

$(references_path)hg38.fa:
	@echo 'Downloading hg38 reference genome from UCSC Golden Path'
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
		-O $(references_path)hg38.fa.gz;
	zcat $(references_path)hg38.fa.gz >$(references_path)hg38.fa;
	rm $(references_path)hg38.fa.gz;
	samtools faidx $(references_path)hg38.fa;
	chmod -w $(references_path)hg38.fa;

$(references_path)hg38.gtf:
	@echo 'Downloading hg38 GTF annotation from ENSEMBL database'
	wget ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz \
		-O $(references_path)hg38.gtf.gz;
	zcat $(references_path)hg38.gtf.gz > $(references_path)hg38.gtf;
	rm $(references_path)hg38.gtf.gz;
	chmod -w $(references_path)hg38.gtf;

$(panel): $(annotation)
ifeq ($(gene_list_name),)
	@echo "Making empty Panel"
	mkdir -p $(panel_prefix)
	touch $(panel)
	chmod -w $(panel);
else
	@echo "Making panel"
	mkdir -p $(panel_prefix)
	$(python3) $(simulating_path)generate_panel.py \
		--gene-annotation $(annotation) \
		--number-of-genes $(num_genes) \
		--gene-list $(gene_list) \
		--random-seed $(random_seed) \
		--output-panel $(panel);
	chmod -w $(panel);
endif


$(molecules): $(reference) $(panel)
	@echo $(panel)
	@echo "Simulating molecules"
	mkdir -p $(molecules_prefix)
	$(python3) $(simulating_path)generate_molecules.py \
		--reference $(reference) \
		--number-of-molecules $(num_molecules) \
		--molecule-size-mean $(molecule_size_mu) \
		--molecule-size-standard-dev $(molecule_size_dev) \
		--min-molecule-size $(min_molecule_size) \
		--random-seed $(random_seed) \
		--output-molecules $(molecules) \
		$(bed_flag);
	chmod -w $(molecules);

$(barcoded_molecules): $(barcodes) $(molecules)
	@echo "Simulating barcoded molecules"
	mkdir -p $(barcoded_molecules_prefix)
	$(python3) $(simulating_path)attach_barcodes_to_molecules.py \
		--input-barcodes $(barcodes) \
		--input-molecules $(molecules) \
		--random-seed $(random_seed) \
		--output-barcoded-molecules $(barcoded_molecules);
	chmod -w $(barcoded_molecules);

$(amplified_barcoded_molecules): $(barcoded_molecules)
	@echo "Simulating amplified barcoded molecules"
	mkdir -p $(amplified_barcoded_molecules_prefix)
	$(python3) $(simulating_path)pcr_duplication.py \
		--molecules $(barcoded_molecules) \
		--number-of-cycles $(pcr_cycles) \
		--duplication-rate-per-cycle $(pcr_duplication_rate) \
		--error-rate $(pcr_error_rate) \
		--random-seed $(random_seed) \
		--pcr-product $(amplified_barcoded_molecules);
	chmod -w $(amplified_barcoded_molecules);

$(reads_log): $(amplified_barcoded_molecules)
	@echo "Simulating reads with ART Illumina"
	mkdir -p $(reads_prefix)
	$(art_illumina) \
		--seqSys $(sequencing_machine) \
		--amplicon \
		--paired \
		--noALN \
		--in $(amplified_barcoded_molecules) \
		--len $(read_length) \
		--fcov 1 \
		--rndSeed $(random_seed) \
		--out $(reads_prefix) \
		> $(reads_log);
	cat $(forward_reads) | bash $(CONVERT_FASTQ) > $(true_cluster);
	chmod -w $(true_cluster) $(forward_reads) $(reverse_reads) $(reads_log);


$(reverse_reads): $(reads_log)
$(forward_reads): $(reverse_reads)

simulate: $(forward_reads) $(reverse_reads) $(reads_log)

simulate_clean:
	rm -rf $(simulation_datasets_path)randomSeed_*


panel: $(panel)
annotation: $(annotation)
reference: $(reference)
barcodes: $(barcodes)
molecules: $(molecules)
barcoded_molecules: $(barcoded_molecules)
amplified_barcoded_molecules: $(amplified_barcoded_molecules)


cluster: calib $(forward_reads) $(reverse_reads)
	$(current_dir)calib \
		--input-forward $(forward_reads) \
		--input-reverse $(reverse_reads) \
		--output-prefix $(calib_output_prefix) \
		--barcode-length $(barcode_length) \
		--ignored-sequence-prefix-length $(ignored_sequence_prefix_length) \
		--minimizer-count $(minimizers_num) \
		--kmer-size $(kmer_size) \
		--error-tolerance $(barcode_error_tolerance) \
		--minimizer-threshold $(minimizers_threshold) \
		--threads $(thread_count) \
		$(silent) \
		# $(no_sort)

accuracy: $(cluster_file) $(true_cluster)
	$(python3) $(simulating_path)rand_index.py \
		--input-cluster-file $(cluster_file) \
		--input-amplified-molecule $(input_amplified_molecules) \
		--output-accuracy-results $(output_accuracy_results)

clean:
	rm -f $(current_dir)calib
