//
// Created by borabi on 19/12/17.
//

#include "extract.h"
#include "global.h"
#include <algorithm>

using namespace std;


#define ASCII_SIZE 128
#define BYTE_SIZE 8
#define ENCODE_SIZE 2
#define INVALID_A 0x00
#define INVALID_C 0x15
#define INVALID_G 0x2A
#define INVALID_T 0x3F
#define INVALID_LENGTH 3

node_id_t node_count = 0;
read_id_t read_count = 0;

// node_to_read_id find reads with identical barcodes and minimizers
typedef std::unordered_map<Node, std::vector<read_id_t>, NodeHash, NodeEqual> node_to_read_id_unordered_map;

read_vector reads;
node_vector nodes;
node_id_to_read_id_vector node_to_read_vector;

minimizer_t encode [ASCII_SIZE];
minimizer_t invalid_kmer = (minimizer_t) -1;
vector<bool> invalid_minimizers;

// bool invalid_minimizer(minimizer_t kmer) {
//     return ((((kmer & 0x3F ) == 0x00) || ((kmer & 0x3F ) == 0x15) || ((kmer & 0x3F ) == 0x2A) || ((kmer & 0x3F ) == 0x3F)) ||
//     (((kmer & 0xFC ) == 0x00) || ((kmer & 0xFC ) == 0x54) || ((kmer & 0xFC ) == 0xA8) || ((kmer & 0xFC ) == 0xFC)));
// }

string minimizer_t_to_dna(minimizer_t minimizer, size_t size) {
    string s = "";
    // printf("0x%x\t", minimizer);
    for (int i = 0; i < (int) size; i++) {
        // cout << "\t" << (minimizer & 0x3);

        switch (minimizer & 0x3) {
        case 0x0:
            s = "A" + s;
            break;
        case 0x1:
            s = "C" + s;
            break;
        case 0x2:
            s = "G" + s;
            break;
        case 0x3:
            s = "T" + s;
            break;
        }
        // cout << "-" + s << "\t";
        minimizer = minimizer >> ENCODE_SIZE;
    }
    // cout << "\n";
    return s;
}

void make_invalid_minimizer_vector() {
    minimizer_t max_minimizer = 1;
    for (int i = 0; i < kmer_size*ENCODE_SIZE; i++) {
        max_minimizer *= 2;
    }
    invalid_minimizers = vector<bool>((size_t) max_minimizer, false);
    minimizer_t triplet_mask = (minimizer_t) -1;
    triplet_mask = triplet_mask >> (sizeof(minimizer_t)*BYTE_SIZE - INVALID_LENGTH*ENCODE_SIZE);
    max_minimizer--;
    while (max_minimizer != (minimizer_t)-1) {
        minimizer_t temp_minimizer = max_minimizer;

        for (int i = 0; i < kmer_size - INVALID_LENGTH + 1; i++) {
            // printf("0x%x\t0x%x\t0x%x\n", max_minimizer, temp_minimizer, temp_minimizer & triplet_mask);
            // cout << minimizer_t_to_dna(max_minimizer, kmer_size - i) << "\t" <<  minimizer_t_to_dna(temp_minimizer, kmer_size - i) << "\t" <<  minimizer_t_to_dna(temp_minimizer & triplet_mask, INVALID_LENGTH) << "\n";
            switch (temp_minimizer & triplet_mask) {
            case INVALID_A:
                invalid_minimizers[max_minimizer] = true;
                break;
            case INVALID_C:
                invalid_minimizers[max_minimizer] = true;
                break;
            case INVALID_G:
                invalid_minimizers[max_minimizer] = true;
                break;
            case INVALID_T:
                invalid_minimizers[max_minimizer] = true;
                break;
            }
            temp_minimizer = temp_minimizer >> ENCODE_SIZE;
        }
        max_minimizer--;
    }
}

void extract_barcodes_and_minimizers() {
    node_to_read_id_unordered_map node_to_read_map;

    for (int i = 0; i < ASCII_SIZE; i++) {
        encode[i] = (minimizer_t) -1;
    }
    encode['A'] = 0x0;
    encode['C'] = 0x1;
    encode['G'] = 0x2;
    encode['T'] = 0x3;
    encode['a'] = 0x0;
    encode['c'] = 0x1;
    encode['g'] = 0x2;
    encode['t'] = 0x3;
    if (no_triplets) {
        make_invalid_minimizer_vector();
    } else {
        minimizer_t max_minimizer = 1;
        for (int i = 0; i < kmer_size*ENCODE_SIZE; i++) {
            max_minimizer *= 2;
        }
        invalid_minimizers = vector<bool>((size_t) max_minimizer, false);
    }

    // for (minimizer_t i =0; i < invalid_minimizers.size(); i++) {
    //     printf("%s\t%s\n", minimizer_t_to_dna(i, kmer_size).c_str(), invalid_minimizers[i]? "true" : "false");
    // }
    ifstream fastq1;
    ifstream fastq2;
    fastq1.open (input_1);
    fastq2.open (input_2);

    dog << "Reading fastq files\n";

    string trash;
    Node current_node;
    reads.push_back(Read());

    // Processing FASTQ files one read at a time
    while (getline(fastq1, reads.back().name_1)) {
        getline(fastq1, reads.back().sequence_1);
        getline(fastq1, trash);
        if (keep_qual) {
            getline(fastq1, reads.back().quality_1);
        } else {
            getline(fastq1, trash);
            reads.back().quality_1 = "Q1";
        }
        getline(fastq2, reads.back().name_2);
        getline(fastq2, reads.back().sequence_2);
        getline(fastq2, trash);
        if (keep_qual) {
            getline(fastq2, reads.back().quality_2);
        } else {
            getline(fastq2, trash);
            reads.back().quality_2 = "Q2";
        }

        int s1_length = reads.back().sequence_1.size();
        int s2_length = reads.back().sequence_2.size();

        // Extracting the barcode from the start of both mates
        if (s1_length >= barcode_length && s2_length >= barcode_length) {
            current_node.barcode =
                reads.back().sequence_1.substr(0, barcode_length) +
                reads.back().sequence_2.substr(0, barcode_length);
        } else {
            current_node.barcode = string (barcode_length*2, 'N');
        }

        s1_length -= barcode_length;
        s2_length -= barcode_length;

        s1_length -= ignored_sequence_prefix_length;
        s2_length -= ignored_sequence_prefix_length;

        // Splitting the remaining sequence into ~ equally sized segments, and extracting minimizers from each
        int s1_seg_length = s1_length / minimizer_count;
        if (s1_seg_length >= kmer_size) {
            int start = barcode_length + ignored_sequence_prefix_length;
            for (int i = 0; i < minimizer_count; i++) {
                current_node.minimizers_1[i] = minimizer(reads.back().sequence_1, start, s1_seg_length);
                start += s1_seg_length;
                if (debug) {
                    if (current_node.minimizers_1[i] != invalid_kmer) {
                        reads.back().sequence_1 += "-" + minimizer_t_to_dna(current_node.minimizers_1[i], kmer_size);
                    } else {
                        reads.back().sequence_1 += "-" + current_node.minimizers_1[i];
                    }
                }
            }
        } else {
            for (int i = 0; i < minimizer_count; i++) {
                current_node.minimizers_1[i] = -1;
            }
        }

        int s2_seg_length = s2_length / minimizer_count;
        if (s2_seg_length >= kmer_size) {
            int start = barcode_length + ignored_sequence_prefix_length;
            for (int i = 0; i < minimizer_count; i++) {
                current_node.minimizers_2[i] = minimizer(reads.back().sequence_2, start, s2_seg_length);
                start += s2_seg_length;
                if (debug) {
                    if (current_node.minimizers_2[i] != invalid_kmer) {
                        reads.back().sequence_2 += "-" + minimizer_t_to_dna(current_node.minimizers_2[i], kmer_size);
                    } else {
                        reads.back().sequence_2 += "-" + current_node.minimizers_2[i];
                    }
                }

            }
        } else {
            for (int i = 0; i < minimizer_count; i++) {
                current_node.minimizers_2[i] = -1;
            }
        }


        if (node_to_read_map.find(current_node) != node_to_read_map.end()) {
            node_to_read_map[current_node].push_back(reads.size()-1);
        } else {
            node_to_read_map.emplace(current_node, vector<node_id_t>{reads.size()-1});
            // vector<node_id_t> current_vector;
            // current_vector.push_back(read_count++);
            // node_to_read[current_node] = current_vector;
            current_node = Node();
        }
        reads.push_back(Read());
//        cout << current_node.id << "\t" << current_node.barcode << "\t" ;
//        for (int i =0; i < minimizer_count; i++)
//            cout << current_node.minimizers_1[i] << "\t";
//        for (int i =0; i < minimizer_count; i++)
//            cout << current_node.minimizers_2[i] << "\t";
//        cout << "\n";
    }

    reads.pop_back();

    read_count = reads.size();
    node_count = node_to_read_map.size();

    if (!silent) {
        cout << "Read count: " << read_count << "\n";
    }
    dog << "Read count: " << read_count << "\n";
    if (!silent) {
        cout << "Node count: " << node_count << "\n";
    }
    dog << "Node count: " << node_count << "\n";

    nodes.reserve(node_count);
    node_to_read_vector.reserve(node_count);
    for (auto kv : node_to_read_map) {
        nodes.push_back(move(kv.first));
        node_to_read_vector.push_back(move(kv.second));
    }

    node_to_read_map.clear();
}

// Extract the lexicographically minimum k-mer in a given string with start and range
minimizer_t minimizer(string& seq, int start, int length){
    int end = start + length;
    // printf("%d\t%d\t%s\n", start, end, seq.substr(start, length).c_str());;
    minimizer_t current_k_mer = (minimizer_t) -1;
    // Biggest possible k-mer is all 1's
    minimizer_t min_k_mer = (minimizer_t) -1;
    // printf("MIN\t0x%08X\n", min_k_mer);
    // Ignoring leftmost bits depending on k-mer size
    minimizer_t kmer_size_mask = (minimizer_t) -1;
    // TODO: Compute this once?
    kmer_size_mask >>= (sizeof(minimizer_t)*BYTE_SIZE-kmer_size*ENCODE_SIZE);
    // printf("MASK\t0x%08X\n", kmer_size_mask);

    if (start + kmer_size > end) {
        // printf("Oops\n");
        return --invalid_kmer;
    }

    // Building the first k-mer
    // printf("i  SEQ\tCURRENT__\tENCODE___\tMASKED__\n");
    for (int i = start; i < start + kmer_size && i < end; i++) {
        // printf("%2d   %c\t", i, seq.at(i));
        // printf("0x%08X\t", current_k_mer);
        // printf("0x%08X\t", encode[(size_t)seq.at(i)]);
        current_k_mer <<= ENCODE_SIZE;
        current_k_mer |= encode[(size_t)seq.at(i)];
        // printf("0x%08X\n", current_k_mer&kmer_size_mask);
        // Hit a non ACTG nucleotide
        if (encode[(size_t)seq.at(i)] == (minimizer_t) -1) {
            current_k_mer = (minimizer_t) -1;
            // printf("HERE\n");
            start = i+1;
            // Hit yet another non nucleotide
            if (start + kmer_size > end) {
                // printf("Oops\n");
                return --invalid_kmer;
            }

        }
    }
    current_k_mer &= kmer_size_mask;

    min_k_mer = ((min_k_mer < current_k_mer) || invalid_minimizers[current_k_mer]) ? min_k_mer : current_k_mer;
    // min_k_mer = min_k_mer < current_k_mer ? min_k_mer : current_k_mer;

    // Bit shifting to get new k-mers, and masking to keep k-mer size fixed
    // printf("MIN\t0x%08X\n", min_k_mer);

    // printf("After first minimizer\n");
    // printf("i  SEQ\tCURRENT__\tENCODE___\tMASKED__\tMINIMUM\n");

    for (int i = start + kmer_size; i < end; i++) {
        // printf("%2d   %c\t", i, seq.at(i));
        // printf("0x%08X\t", current_k_mer);
        // printf("0x%08X\t", encode[(size_t)seq.at(i)]);

        current_k_mer <<= ENCODE_SIZE;
        current_k_mer &= kmer_size_mask;
        current_k_mer |= encode[(size_t)seq.at(i)];
        // printf("0x%08X\t", current_k_mer&kmer_size_mask);

        // Re-building the first k-mer if we hit a non ACTG nucleotide
        if (encode[(size_t)seq.at(i)] == (minimizer_t) -1) {
            // printf("HERE\n");
            current_k_mer = (minimizer_t) -1;
            i++;
            int j;
            for (j = i; j < i + kmer_size && j < end; j++) {
                // printf("j %2d:%c\t", j, seq.at(j));
                // printf("0x%08X\t", current_k_mer);
                // printf("0x%08X\t", encode[(size_t)seq.at(j)]);

                current_k_mer <<= ENCODE_SIZE;
                current_k_mer |= encode[(size_t)seq.at(j)];
                // printf("0x%08X\n", current_k_mer&kmer_size_mask);

                // Hit a non ACTG nucleotide
                if (encode[(size_t)seq.at(j)] == (minimizer_t) -1) {
                    // printf( "==\tHERE\n");
                    i = j+1;
                }
            }
            current_k_mer &= kmer_size_mask;
            i=j;
        }

        min_k_mer = ((min_k_mer < current_k_mer) || invalid_minimizers[current_k_mer]) ? min_k_mer : current_k_mer;
        // min_k_mer = min_k_mer < current_k_mer ? min_k_mer : current_k_mer;
        // printf("0x%08X\n", min_k_mer);

    }

    return min_k_mer != (minimizer_t) -1 ? min_k_mer : --invalid_kmer;
    // return min_k_mer;
}
