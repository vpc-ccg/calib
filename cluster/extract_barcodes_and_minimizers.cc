//
// Created by borabi on 19/12/17.
//

#include "extract_barcodes_and_minimizers.h"
#include "global.h"

using namespace std;

int node_count = 0;
int read_count = 0;
node_map node_to_read;


void extract_barcodes_and_minimizers() {
    ifstream fastq1;
    ifstream fastq2;
    fastq1.open (input_prefix + "1.fq");
    fastq2.open (input_prefix + "2.fq");

    log << "Reading fastq files\n";

    string r1, s1, q1, r2, s2, q2, trash;

    Node current_node;

    cout << "#nodes\t" << node_count << "\n";

    // Processing FASTQ files one read at a time
    while (getline(fastq1, r1)) {
        getline(fastq1, s1);
        getline(fastq1, trash);
        getline(fastq1, q1);
        getline(fastq2, r2);
        getline(fastq2, s2);
        getline(fastq2, trash);
        getline(fastq2, q2);

        int s1_length = s1.size();
        int s2_length = s2.size();

        // Extracting the barcode from the start of both mates
        if (s1_length >= barcode_length && s2_length >= barcode_length){
            current_node.barcode = s1.substr(0, barcode_length) + s2.substr(0, barcode_length);
        } else {
            current_node.barcode = string (barcode_length*2, 'N');
        }

        s1_length -= barcode_length;
        s2_length -= barcode_length;


        // Splitting the remaining sequence into ~ equally sized segments, and extracting minimizers from each
        int s1_seg_length = s1_length / minimizer_count;
        if (s1_seg_length >= kmer_size){
            int start = barcode_length;
            for (int i = 0; i < minimizer_count; i++){
                current_node.minimizers_1[i] = minimizer(s1, start, s1_seg_length);
                start += s1_seg_length;
            }
        } else {
            for (int i = 0; i < minimizer_count; i++){
                current_node.minimizers_1[i] = -1;
            }
        }

        int s2_seg_length = s2_length / minimizer_count;
        if (s2_seg_length >= kmer_size){
            for (int i = 0; i < minimizer_count; i++){
                current_node.minimizers_2[i] = minimizer(s2, barcode_length + i*s2_seg_length, s2_seg_length);
            }
        } else {
            for (int i = 0; i < minimizer_count; i++){
                current_node.minimizers_2[i] = -1;
            }
        }




        if (node_to_read.find(current_node) != node_to_read.end()) {
            node_to_read[current_node].push_back(read_count++);
        } else {
            vector<int> current_vector;
            current_vector.push_back(read_count++);
            node_to_read[current_node] = current_vector;
            current_node = Node();
        }
//
//        cout << current_node.id << "\t" << current_node.barcode << "\t" ;
//        for (int i =0; i < minimizer_count; i++)
//            cout << current_node.minimizers_1[i] << "\t";
//        for (int i =0; i < minimizer_count; i++)
//            cout << current_node.minimizers_2[i] << "\t";
//        cout << "\n";
    }

    delete [] current_node.minimizers_1;
    delete [] current_node.minimizers_2;

}

// Extract the lexicographically minimum k-mer in a given string with start and range
uint64_t minimizer(string& seq, int start, int length){
    // Report -1 if no k-mer fits
    if (length < kmer_size){
        return (uint64_t) - 1;
    }
    uint64_t current_k_mer = 0;
    // Biggest possible k-mer is all 1's
    uint64_t min_k_mer = (uint64_t) - 1;
    // Ignoring leftmost bytes depending on k-mer size
    uint64_t kmer_size_mask = (uint64_t) - 1;
    kmer_size_mask >>= (sizeof(uint64_t)-kmer_size)*8;

    // Building the first k-mer
    for (int i = start; i < start + kmer_size; i++){
        current_k_mer <<= 8;
        current_k_mer |= (uint64_t) seq.at(i);
    }
    min_k_mer = min_k_mer < current_k_mer ? min_k_mer : current_k_mer;

    // Bit shifting to get new k-mers, and masking to keep k-mer size fixed
    for (int i = start + kmer_size; i < start + length; i++){
        current_k_mer <<= 8;
        current_k_mer &= kmer_size_mask;
        current_k_mer |= (uint64_t) seq.at(i);
        min_k_mer = min_k_mer < current_k_mer ? min_k_mer : current_k_mer;
    }

    return min_k_mer;
}
