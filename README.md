# BarGoat
BarGoat clusters paired-end reads using their barcodes and sequences. BarGoat is suitable for amplicon sequencing where a molecule is tagged, then PCR amplified with high depth.

## Prerequisites

BarGoat is implemented in Python 3 and C++11. Clustering reads requires installing the following dependencies:

- [NumPy](http://www.numpy.org/) 
- [argparse](https://pypi.python.org/pypi/argparse)
- [igraph](http://igraph.org/python/)

Note that igraph installation has some C dependencies. Make sure igraph is working before running BarGoat. We personally recommend using [Anaconda 3](https://docs.anaconda.com/anaconda/install/) to manage all these dependencies.

If you wish to simulate read sets using BarGoat, you will need the following prerequisites:

- [pyfaidx](https://pythonhosted.org/pyfaidx/)
- [argparse](https://pypi.python.org/pypi/argparse)
- [ART Illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) (version 2.5.8)

As of the time of updating this README, Anaconda's ART package is outdated. For that, you either have to compile your own or, even better, use [Linuxbrew](http://linuxbrew.sh/). Once you have Linuxbrew installed, you simply can run:

```bash
brew install homebrew/science/art
```



## BarGoat "Hello world..." Run

The simplest way to run BarGoat is using its Makefile as an executable. You can modify any of the command-line parameters by passing them to the executable. Let us first try simulating an amplicon sequencing run. This can be easily done using the `simulate` command:

```bash
./bargoat simulate
```

This will generate the following files under Data subdirectory:

- `simulated_barcodes.txt`: List of 10000 barcodes each of length 10 each in a new line
- `simulated_molecules.fa`: FASTA file of 1000 molecules with mean length of 250, standard deviation of 25, and minimum length of 150 generated randomly for E. Coli FASTA file (included with BarGoat).
- `simulated_barcoded_molecules.fa`: FASTA file of the same molecules prepended and appended by a random barcode from `simulated_barcodes.txt`.
- `simulated_reads_1.fq` and `simulated_reads_2.fq`: Paired-end reads generated from `simulated_barcoded_molecules.fa` by ART Illumina

Note that any of the specified parameters can easily be changed by passing to the Makefile. For example:

```bash
./bargoat simulate molecule_size_mu=400
```

Will generate molecules with mean size of 400 nucleotides.

Now let us try clustering the simulated reads:

```bash
./bargoat cluster
```

By default, `cluster` will run on the simulated reads. This will output the following files:

- `simulated_reads_l10_m3_k8_e2_t1_q1.0.clusters`: A tab delimited file of the clusters. The header of each cluster starts with a `#` sign followed without a space with the cluster ID (incremental counter). The header then includes the number of nodes, the number of edges, the number of reads, and the graph density (function of number of nodes and edges) of the cluster.
- `simulated_reads_l10_m3_k8_e2_t1_q1.0.log`: A log file containing what was printed to std.out
- `simulated_reads_l10_m3_k8_e2_t1_q1.0.tsv`: A tab delimited file that contains the extracted barcode and minimizers of the reads.

BarGoat `cluster` has some parameters. Most important of them are the following six:

- `barcode_length`: Barcode length on a single side of the read. Default is 10.

- `minimizers_num`: Number of minimizers to extract per mate. Default is 3.
- `kmer_size`: Size of each minimizer. Default is 3.
- `barcode_error_tolerance`: Max hamming distance between two barcodes (of length `barcode_length` *2) to be counted as similar. Default is 2.
- `minimizers_threshold`: Minimum number of matching minimizers on each mate to be considered similar. Default is 1.
- `ratio`: The ratio of the subsample of nodes to consider for each template of barcode similarity checking. Default is 1.0.

You may note that these parameters are included in the output prefix of all `cluster` output files for convenience. 

Finally, if you ran a simulated dataset, you can check BarGoat `cluster` accuracy using the `accuracy` command:

```bash
./bargoat accuracy
```

Which will take the default simulated reads clusters produced by `cluster` command. Accuracy is measured using [Rand Index](https://en.wikipedia.org/wiki/Rand_index).

## Reproducing Benchmarks

To reporucede any of the benchmarks in our report, run a command like this:

```bash
./bargoat simulate cluster accuracy num_molecules=300000 num_barcodes=30000 barcode_length=8 minimizers_num=5 kmer_size=4 barcode_error_tolerance=1 minimizers_threshold=3 
```
Where you may change the value of any parameters to match any of the reported datasets (or a new dataset for that matter).

## Report
Our report is hosted on Overleaf at [https://www.overleaf.com/read/jywsyjsmtrdp](https://www.overleaf.com/read/jywsyjsmtrdp)





