# Calib
Calib clusters paired-end reads using their barcodes and sequences. Calib is suitable for amplicon sequencing where a molecule is tagged, then PCR amplified with high depth.

## Reproducing Benchmarks

If you are insterested in reprodcuing our paper's benchmarks, please switch to the paper git branch, and follow instructions at [Paper Branch](https://github.com/vpc-ccg/calib/tree/paper/)


## Prerequisites

Calib read clustering is implemented in C++11 and has been tested on Linux operating system.

The simulatation module of Calib is implemented in Python 3. Calib read simulation has the following prerequisites:

- [pyfaidx](https://pypi.python.org/pypi/pyfaidx)
- [numpy](https://pypi.python.org/pypi/numpy)
- [scipy](https://pypi.python.org/pypi/scipy)
- [scikit-learn](https://pypi.python.org/pypi/scikit-learn)
- [biopython](https://pypi.python.org/pypi/biopython)
- [pandas](https://pypi.python.org/pypi/pandas)
- [ART Illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) (version 2.5.8)

All these prerequisites can easily be satisfied using Anaconda.



## Calib "Hello world..." Run

The simplest way to run Calib is using its Makefile as an executable. You can modify any of the command-line parameters by passing them to the executable. Let us first try simulating an amplicon sequencing run. This can be easily done using the `simulate` command:

```bash
./calib simulate
```

This will generate the following files under `simulating/datasets/` subdirectory:

- `rs_42.bl_8.nb_100.barcodes.txt`: List of 100 barcode tags each of length 8 each in a new line
- `rs_42.ref_e_coli.bed_NA.mu_200.dev_25.nm500.molecules.fa`: FASTA file of 500 molecules with mean length of 200, standard deviation of 25 generated from uniformly random starting positions from E. Coli FASTA file (included with Calib).
- `rs_42.bl_8.nb_100.ref_e_coli.bed_NA.mu_200.dev_25.nm500.barcoded_molecules.fa`: FASTA file of the same molecules prepended and appended by a randomly sampled barcode tag with replacement from `rs_42.bl_8.nb_100.barcodes.txt`.
- `rs_42.bl_8.nb_100.ref_e_coli.bed_NA.mu_200.dev_25.nm500.pc_7.pdr_0.6.per_0.000005.amplified_barcoded_molecules.fa`: FASTA file of the PCR product of the generated molecules with 7 PCR cycles, 0.6 efficiency rate, 0.000005 substitution error rate.
- `rs_42.bl_8.nb_100.ref_e_coli.bed_NA.mu_200.dev_25.nm500.pc_7.pdr_0.6.per_0.000005.sm_HS20.1.fq` and `rs_42.bl_8.nb_100.ref_e_coli.bed_NA.mu_200.dev_25.nm500.pc_7.pdr_0.6.per_0.000005.sm_HS20.2.fq`: Paired-end reads generated from the PCR product by ART Illumina simulating HiSeq 2000.

Note that any of the specified parameters can easily be changed by passing to the Makefile. For example:

```bash
./calib simulate molecule_size_mu=400
```

Will generate molecules with mean size of 400 nucleotides. The main simulating parameters are:

- random_seed (int)'
- bed (bed file should be in simulating/genomes/\<bed\>.bed)'
- reference_name (fasta file should be in simulating/genomes/<reference_name>.fa)'
- num_barcodes (int)'
- barcode_length (int)'
- molecule_size_mu (int)'
- molecule_size_dev (int)'
- num_molecules (int)'
- read_length (int)'
- pcr_cycles (int)'
- pcr_duplication_rate (float (0,1))'
- pcr_error_rate (float (0,1))'
- sequencing_machine (string from (HS10, HS20, HS25, HSXn, HSXt, MinS, MSv1, MSv3, NS50))'



Now let us try clustering the simulated reads:

```bash
./calib cluster
```

By default, `cluster` will run on the simulated reads. This will output the following files:

- `rs_42.bl_8.nb_100.ref_e_coli.bed_NA.mu_200.dev_25.nm500.pc_7.pdr_0.6.per_0.000005.sm_HS20.calib.bl_8.mn_3.ks_8.bet_2.mt_1.cluster`: A tab delimited file of the clusters. The header of each cluster starts with a `#` sign followed without a space with the cluster ID (incremental counter). Each line following that is a read pair belonging to this cluster with the following columns:
-- Node ID: The node in the graph that represented this read
-- Read ID: The order of the read in the input FASTQ files.
-- Read Name 1
-- Read Sequence 1
-- Read Quality 1 (if not kept, then Q1)
-- Read Name 2
-- Read Sequence 2
-- Read Quality 2 (if not kept, then Q2)
- `rs_42.bl_8.nb_100.ref_e_coli.bed_NA.mu_200.dev_25.nm500.pc_7.pdr_0.6.per_0.000005.sm_HS20.calib.bl_8.mn_3.ks_8.bet_2.mt_1.cluster.log`: A log file containing what was printed to std.out in addition to time logging

Calib `cluster` has some parameters. Most important of them are the following:
- forward_reads (string)'
- reverse_reads (string)'
- output_prefix (string)'
- barcode_length (int)'
- ignored_sequence_prefix_length (int)'
- minimizers_num (int)'
- kmer_size (int)'
- barcode_error_tolerance (int < barcode length)'
- minimizers_threshold (int < number of minimizers)'
- silent (no value)'

Finally, if you ran a simulated dataset, you can check Calib `cluster` accuracy using the `accuracy` command:

```bash
./calib accuracy
```

The results of accuracy are stored in `rs_42.bl_8.nb_100.ref_e_coli.bed_NA.mu_200.dev_25.nm500.pc_7.pdr_0.6.per_0.000005.sm_HS20.calib.bl_8.mn_3.ks_8.bet_2.mt_1.accuracy`. Note that accuracy is measured using [justed Rand Index](https://en.wikipedia.org/wiki/Rand_index#Adjusted_Rand_index).
