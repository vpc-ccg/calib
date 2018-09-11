# Calib parameter selection

Calib's clustering input is paired-end FASTQ files. Calib currently only supports paired-end reads, but support for single-end is planned soon. Calib does not support interleaved FASTQ files (yet). The fixed length barcode tags of each read should be at the beginning of the sequence of either read mates (no adapters should be present in the read sequences).

### Algorithm

Calib clusters reads by defining barcode similarity and sequence similarity between pair of reads. Two reads barcodes are considered similar if they are no more than `e` [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance) from one another. Sequence similarity is defined in terms of the sequence minimizers. Minimizers are

