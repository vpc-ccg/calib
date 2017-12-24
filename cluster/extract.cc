//
// Created by borabi on 19/12/17.
//

#include "extract.h"
#include "global.h"

using namespace std;


#define ASCII_SIZE 128
#define BYTE_SIZE 8
#define ENCODE_SIZE 2

int node_count = 0;
int read_count = 0;
node_map node_to_read;
uint64_t encode [ASCII_SIZE];

kmer_t invalid_kmer = (kmer_t) -1;



void extract_barcodes_and_minimizers() {

    for (int i = 0; i < ASCII_SIZE; i++)
        encode[i] = (uint64_t) -1;
    encode['A'] = 0x0;
    encode['C'] = 0x1;
    encode['G'] = 0x2;
    encode['T'] = 0x3;
    encode['a'] = 0x0;
    encode['c'] = 0x1;
    encode['g'] = 0x2;
    encode['t'] = 0x3;

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
kmer_t minimizer(string& seq, int start, int length){
    int end = start + length;
    // printf("%d\t%d\t%s\n", start, end, seq.substr(start, length).c_str());;
    kmer_t current_k_mer = (kmer_t) -1;
    // Biggest possible k-mer is all 1's
    kmer_t min_k_mer = (kmer_t) -1;
    // printf("MIN\t0x%08X\n", min_k_mer);
    // Ignoring leftmost bits depending on k-mer size
    kmer_t kmer_size_mask = (kmer_t) -1;
    // TODO: Compute this once?
    kmer_size_mask >>= (sizeof(kmer_t)*BYTE_SIZE-kmer_size*ENCODE_SIZE);
    // printf("MASK\t0x%08X\n", kmer_size_mask);

    if (start + kmer_size > end){
        // printf("Oops\n");
        return invalid_kmer--;
    }

    // Building the first k-mer
    // printf("i  SEQ\tCURRENT__\tENCODE___\tMASKED__\n");
    for (int i = start; i < start + kmer_size && i < end; i++){
        // printf("%2d   %c\t", i, seq.at(i));
        // printf("0x%08X\t", current_k_mer);
        // printf("0x%08X\t", encode[(size_t)seq.at(i)]);
        current_k_mer <<= ENCODE_SIZE;
        current_k_mer |= encode[(size_t)seq.at(i)];
        // printf("0x%08X\n", current_k_mer&kmer_size_mask);
        // Hit a non ACTG nucleotide
        if (encode[(size_t)seq.at(i)] == (kmer_t) -1){
            current_k_mer = (kmer_t) -1;
            // printf("HERE\n");
            start = i+1;
            // Hit yet another non nucleotide
            if (start + kmer_size > end){
                // printf("Oops\n");
                return invalid_kmer--;
            }

        }
    }
    current_k_mer &= kmer_size_mask;

    min_k_mer = min_k_mer < current_k_mer ? min_k_mer : current_k_mer;

    // Bit shifting to get new k-mers, and masking to keep k-mer size fixed
    // printf("MIN\t0x%08X\n", min_k_mer);

    // printf("After first minimizer\n");
    // printf("i  SEQ\tCURRENT__\tENCODE___\tMASKED__\tMINIMUM\n");

    for (int i = start + kmer_size; i < end; i++){
        // printf("%2d   %c\t", i, seq.at(i));
        // printf("0x%08X\t", current_k_mer);
        // printf("0x%08X\t", encode[(size_t)seq.at(i)]);

        current_k_mer <<= ENCODE_SIZE;
        current_k_mer &= kmer_size_mask;
        current_k_mer |= encode[(size_t)seq.at(i)];
        // printf("0x%08X\t", current_k_mer&kmer_size_mask);

        // Re-building the first k-mer if we hit a non ACTG nucleotide
        if (encode[(size_t)seq.at(i)] == (kmer_t) -1){
            // printf("HERE\n");
            current_k_mer = (kmer_t) -1;
            i++;
            int j;
            for (j = i; j < i + kmer_size && j < end; j++){
                // printf("j %2d:%c\t", j, seq.at(j));
                // printf("0x%08X\t", current_k_mer);
                // printf("0x%08X\t", encode[(size_t)seq.at(j)]);

                current_k_mer <<= ENCODE_SIZE;
                current_k_mer |= encode[(size_t)seq.at(j)];
                // printf("0x%08X\n", current_k_mer&kmer_size_mask);

                // Hit a non ACTG nucleotide
                if (encode[(size_t)seq.at(j)] == (kmer_t) -1){
                    // printf( "==\tHERE\n");
                    i = j+1;
                }
            }
            current_k_mer &= kmer_size_mask;
            i=j;
        }

        min_k_mer = min_k_mer < current_k_mer ? min_k_mer : current_k_mer;
        // printf("0x%08X\n", min_k_mer);

    }

    return min_k_mer;
}
