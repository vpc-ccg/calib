//
// Created by borabi on 19/12/17.
//

#include "extract_barcodes_and_minimizers.h"

int extract_barcodes_and_minimizers() {
    ifstream fastq1;
    ifstream fastq2;
    fastq1.open (input_prefix + "1.fastq");
    fastq1.open (input_prefix + "2.fastq");
    ofstream output(argv[3]);
    unsigned int barcode_length = (unsigned int) atoi(argv[4]);
    unsigned int k_mer_size = (unsigned int) atoi(argv[5]);
    unsigned int minimizers_count = (unsigned int) atoi(argv[6]);

    string r1, s1, q1, r2, s2, q2, trash;
    // Processing FASTQ files one read at a time
    while (getline(fastq1, r1)) {
        getline(fastq1, s1);
        getline(fastq1, trash);
        getline(fastq1, q1);
        getline(fastq2, r2);
        getline(fastq2, s2);
        getline(fastq2, trash);
        getline(fastq2, q2);

        unsigned int s1_length = (unsigned int) s1.size();
        unsigned int s2_length = (unsigned int) s2.size();
        uint64_t min_kmer;

        // Extracting the barcode from the start of both mates
        if (s1_length < barcode_length){
            output << string (barcode_length, 'N');
        } else {
            output << s1.substr(0, barcode_length);
        }
        if (s2_length < barcode_length){
            output << string (barcode_length, 'N') << "\t";
        } else {
            output << s2.substr(0, barcode_length) << "\t";
        }
        s1_length -= barcode_length;
        s2_length -= barcode_length;

        // Splitting the remaining sequence into ~ equally sized segments, and extracting minimizers from each
        int s1_seg_length = s1_length / minimizers_count;
        for (int i = 0; i < minimizers_count; i++){
            min_kmer = minimizer(s1, barcode_length + i*s1_seg_length, s1_seg_length, k_mer_size);
            output << min_kmer << "\t";
        }

        int s2_seg_length = s2_length / minimizers_count;
        for (int i = 0; i < minimizers_count; i++){
            min_kmer = minimizer(s2, barcode_length + i*s2_seg_length, s2_seg_length, k_mer_size);
            output << min_kmer << "\t";
        }
        output << "\n";
    }
}
