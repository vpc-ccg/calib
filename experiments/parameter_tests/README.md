# Calib parameter selection
If no clustering parameters are passed to Calib, Calib will automatically selects its own parameters based on the the barcode length and inferred read length.
The selection of the parameter was done using extensive testing of wide range of possible paramters.
The testing script is available [here](../../slurm_scripts/).

Basically, we tested 9 different dataset types by using:
- Barcode tag lengths of 4, 8, or 12
- Read mate length of 75, 150, or 250
Then, for each dataset type, we generated 10 different datasets, each using a different random seed.
We then ran Calib on each of these sets using:
- Error tolerance of 1 or 2
- k-mer size of 4 or 8
- Minimizer count of:
  - 3 with thresholds of 1 or 2
  - 4 with thresholds of 1, 2, or 3
  - 5 with thresholds of 2, 3, or 4
  - 6 with thresholds of 2, 3, 4, or 5
  - 7 with thresholds of 2, 3, 4, 5, or 6

The complete results of 68 different parameter sets on the 10 different random seeds of each of the dataset types in this directory's TSV files.
The results are plotted in these [plots](https://cdn.rawgit.com/vpc-ccg/calib/dev/experiments/parameter_tests/parameter_sets_tests_plots.html).
In these plots, we filtered out any parameter set that scored less than 0.99/1.00 in any of the 10 random seeds on a given dataset type.

Based on these results, we preconfigured Calib with the following default parameter selection:
```C++
if (barcode_length >= 1 && barcode_length <= 6) {
    if (mean_read_size >= 61 && mean_read_size <= 100) {
        error_tolerance     = 1;
        kmer_size           = 4;
        minimizer_count     = 6;
        minimizer_threshold = 1;
    }
    if (mean_read_size >= 101 && mean_read_size <= 150) {
        error_tolerance     = 1;
        kmer_size           = 8;
        minimizer_count     = 7;
        minimizer_threshold = 2;
    }
    if (mean_read_size >= 151 && mean_read_size <= 250) {
        error_tolerance     = 1;
        kmer_size           = 8;
        minimizer_count     = 7;
        minimizer_threshold = 2;
    }
}
if (barcode_length >= 7 && barcode_length <= 11) {
    if (mean_read_size >= 61 && mean_read_size <= 100) {
        error_tolerance     = 2;
        kmer_size           = 4;
        minimizer_count     = 7;
        minimizer_threshold = 3;
    }
    if (mean_read_size >= 101 && mean_read_size <= 150) {
        error_tolerance     = 2;
        kmer_size           = 8;
        minimizer_count     = 7;
        minimizer_threshold = 2;
    }
    if (mean_read_size >= 151 && mean_read_size <= 250) {
        error_tolerance     = 2;
        kmer_size           = 8;
        minimizer_count     = 7;
        minimizer_threshold = 2;
    }
}
if (barcode_length >= 12) {
    if (mean_read_size >= 61 && mean_read_size <= 100) {
        error_tolerance     = 2;
        kmer_size           = 4;
        minimizer_count     = 7;
        minimizer_threshold = 3;
    }
    if (mean_read_size >= 101 && mean_read_size <= 150) {
        error_tolerance     = 2;
        kmer_size           = 4;
        minimizer_count     = 7;
        minimizer_threshold = 3;
    }
    if (mean_read_size >= 151 && mean_read_size <= 250) {
        error_tolerance     = 2;
        kmer_size           = 8;
        minimizer_count     = 7;
        minimizer_threshold = 2;
    }
}
```

Note that `mean_read_size` if the average read lenght of the first 10,000 reads in the first FASTQ file.
This preconfiguration is what we deemed suitable tradeoff between accuracy, time, and memory use, with more preference on better accuracy than on better time or memory.
Please feel free to select a different set of parameters based on the plot results if that suits your application more.
