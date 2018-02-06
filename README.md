# Calib
Calib clusters paired-end reads using their barcodes and sequences. Calib is suitable for amplicon sequencing where a molecule is tagged, then PCR amplified with high depth.

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

This will generate the following files under Data subdirectory:

- `simulated_barcodes.txt`: List of 10000 barcodes each of length 10 each in a new line
- `simulated_molecules.fa`: FASTA file of 1000 molecules with mean length of 250, standard deviation of 25, and minimum length of 150 generated randomly for E. Coli FASTA file (included with Calib).
- `simulated_barcoded_molecules.fa`: FASTA file of the same molecules prepended and appended by a random barcode from `simulated_barcodes.txt`.
- `simulated_reads_1.fq` and `simulated_reads_2.fq`: Paired-end reads generated from `simulated_barcoded_molecules.fa` by ART Illumina

Note that any of the specified parameters can easily be changed by passing to the Makefile. For example:

```bash
./calib simulate molecule_size_mu=400
```

Will generate molecules with mean size of 400 nucleotides.

Now let us try clustering the simulated reads:

```bash
./calib cluster
```

By default, `cluster` will run on the simulated reads. This will output the following files:

- `simulated_reads_l10_m3_k8_e2_t1_q1.0.clusters`: A tab delimited file of the clusters. The header of each cluster starts with a `#` sign followed without a space with the cluster ID (incremental counter). The header then includes the number of nodes, the number of edges, the number of reads, and the graph density (function of number of nodes and edges) of the cluster.
- `simulated_reads_l10_m3_k8_e2_t1_q1.0.log`: A log file containing what was printed to std.out
- `simulated_reads_l10_m3_k8_e2_t1_q1.0.tsv`: A tab delimited file that contains the extracted barcode and minimizers of the reads.

Calib `cluster` has some parameters. Most important of them are the following six:

- `barcode_length`: Barcode length on a single side of the read. Default is 10.

- `minimizers_num`: Number of minimizers to extract per mate. Default is 3.
- `kmer_size`: Size of each minimizer. Default is 3.
- `barcode_error_tolerance`: Max hamming distance between two barcodes (of length `barcode_length` *2) to be counted as similar. Default is 2.
- `minimizers_threshold`: Minimum number of matching minimizers on each mate to be considered similar. Default is 1.

You may note that these parameters are included in the output prefix of all `cluster` output files for convenience.

Finally, if you ran a simulated dataset, you can check Calib `cluster` accuracy using the `accuracy` command:

```bash
./calib accuracy
```

Which will take the default simulated reads clusters produced by `cluster` command. Accuracy is measured using [Rand Index](https://en.wikipedia.org/wiki/Rand_index).

## Reproducing Benchmarks

To reporucede any of the benchmarks in our paper, please switch to the paper git branch, and follow instructions at [Paper Branch](https://github.com/vpc-ccg/calib/tree/paper/)
