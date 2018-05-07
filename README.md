# Calib
Calib clusters paired-end reads using their barcodes and sequences. Calib is suitable for amplicon sequencing where a molecule is tagged, then PCR amplified with high depth.

## Reproducing Benchmarks

If you are insterested in reprodcuing our paper's benchmarks, please switch to the paper git branch, and follow instructions at [Paper Branch](https://github.com/vpc-ccg/calib/tree/paper/)

### Clustering Prerequisites
The only prerequisite for Calib clustering is having GCC with version >5.2. 

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

All these prerequisites can easily be satisfied using [Anaconda](https://docs.anaconda.com/anaconda/install/linux).

## Running Calib clustering
To run Calib, you first need to clone Calib using ``git``:
```bash
git clone https://github.com/vpc-ccg/calib.git
```
Then you need to compile Calib:
```bash
cd calib
make
```
Calib executable should now be inside the cloned directory and will be named ``calib``.
To run calib simply run:
```bash
./calib 
```
And the help will printed for you.

## Running Calib simulating
The simulation module can be run using the Makefile:
```bash
make simulate [parameter=value]
```
The simulation parameters are:
- ``random_seed``: integer
- ``reference_name``: reference genome name (excluding .fa extention) for targeted simulation. The reference genome file must be in simulating/genomes directectory. Note that reference ``hg38`` is automatically downloaded.
- ``bed``: bed name (excluding reference_name.bed extention) for targeted simulation. The bed file must be in `simulating/genomes` directectory. An example bed file is included in the directory: `simulating/genomes/Panel.hg38.bed`
- ``num_barcodes``: integer
- ``barcode_length``: integer for the length of half of the barcode (length of the tag)
- ``molecule_size_mu``: integer for average size of the molecule
- ``molecule_size_dev``: integer for standard deviation of the size of the molecule
- ``num_molecules``: integer for the number of molecules to generate.
- ``pcr_cycles``: integer for number of PCR cycles
- ``pcr_duplication_rate``: float between zero and one, exclusively
- ``pcr_error_rate``: float between zero and one, exclusively
- ``sequencing_machine``: string for ART Illumina sequencing platform

For a complete list of default values for these parameters, check the Makefile.
