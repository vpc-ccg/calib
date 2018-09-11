# Calib Experiments Scripts
We performed four experiments for Calib. The results of the experiments are [here](../experiments/). This directory contains scripts for running those experiments on [Slurm Workload Manager](https://slurm.schedmd.com/). However, the Slurm files generated are still BASH files and can be run using:

```bash
bash <slurm_file_name>
```

Note: **All Slurm scripts must be run from Calib's root directory!**

## Simulated Datasets

To run the simulated dataset tests, you need first to clone the other tools by running:

```bash
cd <CALIB_ROOT_DIRECTORY>
git submodule update --init --recursive
```

Then run the Slurm generating and running script:

```bash
cd <CALIB_ROOT_DIRECTORY>
slurum_scripts/simulated_tests.sh <dataset_name>
```

There are three datasets tests reported in Calib's paper, 

- `small`: Has 100K molecules and 100 barcode tags
- `medium`: Has 1M molecules and 5K barcode tags
- `large` Has 1M molecules and 25K barcode tags

The results will be in TSV files off different tools. For `small`:

> `simulating/datasets/randS_42/barL_8.barNum_100/geneNum_35.refName_hg38.geneList_COSMIC_cancer_genes/ref_hg38.molMin_150.molMu_300.molDev_25.molNum100000/pcrC_7.pcrDR_0.6.pcrER_0.00005/seqMach_HS25.readL_150/*_benchmarks.tsv`

For `medium`:

> `simulating/datasets/randS_42/barL_8.barNum_25000/geneNum_35.refName_hg38.geneList_COSMIC_cancer_genes/ref_hg38.molMin_150.molMu_300.molDev_25.molNum1000000/pcrC_7.pcrDR_0.6.pcrER_0.00005/seqMach_HS25.readL_150/*_benchmarks.tsv`

For `large`:

> `simulating/datasets/randS_42/barL_8.barNum_5000/geneNum_35.refName_hg38.geneList_COSMIC_cancer_genes/ref_hg38.molMin_150.molMu_300.molDev_25.molNum1000000/pcrC_7.pcrDR_0.6.pcrER_0.00005/seqMach_HS25.readL_150/*_benchmarks.tsv`

### Note on running things without Slurm

If you wish to run the Slurm scripts on BASH, you need to run one Slurm file at a time. The Slurm files have dependency on one another, and things can break down if the dependency is not respected. This is true especially in case of simulating the datasets. If two scripts end up trying to generate the same reads at the same time, they will end up inevitably corrupting the simulation outputs.

## Real dataset

To run the simulated dataset tests, you need first to clone the other tools by running:

```bash
cd <CALIB_ROOT_DIRECTORY>
git submodule update --init --recursive
```

Then run the Slurm generating and running script:

```bash
cd <CALIB_ROOT_DIRECTORY>
slurm_scripts/real_tests.sh <R1.fastq> <R2.fastq> <panel.hg19.bed> <output_directory>
```

`<R1.fastq>` and `<R2.fastq>` are the real dataset FASTQ files which can be downloaded from [here (BROKEN LINK)]().  `panel.hg19.bed` is the panel of targeted regions used to pull down material for sequencing which is also included in the download link of the FASTQ files. Finally, output directory is where the results will be put of running the complete pipeline.

### Notes

If you are not going to use Slurm, follow the same note mentioned in the simulated dataset runs above.

The real dataset testing pipeline assumes that `samtools` is in your `$PATH`. If `samtools` is installed somewhere else, please edit the corresponding variable in `slurm_scripts/real_tests.sh`.

## Scalability of multithreading

Simply run:

```bash
cd <CALIB_ROOT_DIRECTORY>
slurm_scripts/calib_scalability_tests.sh
```

Feel free to edit the variables and for loop in the script to try different parameters. The same not applies regarding not running Slurm as above.

## Parameter selection

Simply run:

```bash
cd <CALIB_ROOT_DIRECTORY>
slurm_scripts/calib_parameter_tests.sh
```

Feel free to edit the variables and for loop in the script to try different parameters. The same not applies regarding not running Slurm as above.