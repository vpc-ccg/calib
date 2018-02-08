# Calib
Calib clusters paired-end reads using their barcodes and sequences. Calib is suitable for amplicon sequencing where a molecule is tagged, then PCR amplified with high depth.

## Prerequisites

### Simulation Prerequisites
The simulatation module of Calib is implemented in Python 3 and requires that the following Python packages to be installed and importable:

- [pyfaidx](https://pypi.python.org/pypi/pyfaidx)
- [numpy](https://pypi.python.org/pypi/numpy)
- [scipy](https://pypi.python.org/pypi/scipy)
- [scikit-learn](https://pypi.python.org/pypi/scikit-learn)
- [biopython](https://pypi.python.org/pypi/biopython)
- [pandas](https://pypi.python.org/pypi/pandas)

In addition, ART Illumina version 2.5.8 need to be in your `PATH`.
- [ART Illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) (version 2.5.8)

Finally, our tests are run on hg38 reference genome. Please download it and have it in `simulaing/genomes` directory. To do this, assuming you are in Calib's directory:

```bash
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz 
zcat hg38.fa.gz > simulating/genomes/hg38.fa
```

Please ensure the correct naming of the reference genome FASTA file is used.

All these prerequisites can easily be satisfied using [Anaconda](https://docs.anaconda.com/anaconda/install/linux).

### Calib, Rainbow, and starcode Prerequisites
We are benchmarking Calib against [Rainbow](https://github.com/ChongLab/rainbow), [starcode](https://github.com/gui11aume/starcode), and [Du Novo](https://github.com/galaxyproject/dunovo). Rainbow and starcode need to be downloaded, compiled, and their executables to be in your `PATH`. Calib was tested using GCC 5.2, but earlier versions supporting C++11 should work.


### Du Novo Prerequistes
Do Novo requires using Python 2.7, and some old version of samtools. We handled this by a script that creates a special Anaconda environment for Du Novo to run in. All what you need is:
- Install [Anaconda](https://docs.anaconda.com/anaconda/install/linux)
- Download [Du Novo](https://github.com/galaxyproject/dunovo) to `~/bin/dunovo` or make sure to modify the line in `run_dunovo_test.sh` to reflect where dunovo is installed. Make sure to not add `/` at the end:
```bash
dunovo_path=~/bin/dunovo
```

## Running Tests

There are 3 different datasets we preconfigured, and one tiny additional dataset. To run any of those just run:
```bash
./run_tests.sh tiny small medium huge
```
Where you can omit any of the dataset names. To run Du Novo's benchmarking, run:
```bash
./run_dunovo_test.sh tiny small medium huge
```

All results will be in `simulating/datasets/*tsv`. Details of clusterings are also present in the same directory, with decriptive filenames.

