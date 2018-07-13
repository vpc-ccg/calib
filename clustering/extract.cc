//
// Created by borabi on 19/12/17.
//

#include "extract.h"
#include "global.h"
#include <algorithm>

// Debug includes
#include <sstream>
#include <iomanip>

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
barcode_id_t barcode_count = 0;

// node_to_read_id finds reads with identical barcodes and minimizers
typedef unordered_map<Node, vector<read_id_t>, NodeHash, NodeEqual> node_to_read_id_unordered_map;
// barcode_to_node_id finds nodes with identical barcodes
typedef unordered_map<barcode_t, vector<node_id_t> > barcode_to_node_id_unordered_map;

read_id_to_node_id_vector read_to_node_vector;
node_id_to_minimizers_vector node_to_minimizers;
barcode_vector barcodes;
barcode_id_to_node_ids_vector barcode_to_nodes_vector;

minimizer_t encode [ASCII_SIZE];
minimizer_t invalid_kmer = (minimizer_t) -1;

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

    ifstream fastq1;
    ifstream fastq2;
    fastq1.open (input_1);
    fastq2.open (input_2);

    dog << "Reading fastq files\n";

    string name_1, sequence_1, name_2, sequence_2, trash;
    Node current_node;
    read_count = 0;
    // Processing FASTQ files one read at a time
    while (getline(fastq1, name_1)) {
        getline(fastq1, sequence_1);
        getline(fastq1, trash);
        getline(fastq1, trash);
        getline(fastq2, name_2);
        getline(fastq2, sequence_2);
        getline(fastq2, trash);
        getline(fastq2, trash);

        int s1_length = sequence_1.size();
        int s2_length = sequence_2.size();

        // Extracting the barcode from the start of both mates
        if (s1_length >= barcode_length && s2_length >= barcode_length) {
            current_node.barcode =
                sequence_1.substr(0, barcode_length) +
                sequence_2.substr(0, barcode_length);
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
                current_node.minimizers[i] = minimizer(sequence_1, start, s1_seg_length);
                start += s1_seg_length;
            }
        } else {
            for (int i = 0; i < minimizer_count; i++) {
                current_node.minimizers[i] = --invalid_kmer;
            }
        }

        int s2_seg_length = s2_length / minimizer_count;
        if (s2_seg_length >= kmer_size) {
            int start = barcode_length + ignored_sequence_prefix_length;
            for (int i = minimizer_count; i < minimizer_count*2; i++) {
                current_node.minimizers[i] = minimizer(sequence_2, start, s2_seg_length);
                start += s2_seg_length;
            }
        } else {
            for (int i = minimizer_count; i < minimizer_count*2; i++) {
                current_node.minimizers[i] = --invalid_kmer;
            }
        }

        if (node_to_read_map.find(current_node) != node_to_read_map.end()) {
            node_to_read_map[current_node].push_back(read_count);
        } else {
            node_to_read_map.emplace(current_node, vector<read_id_t>{read_count});
            current_node = Node();
        }
        read_count++;
    }
    read_to_node_vector.reserve(read_count);
    node_to_minimizers.reserve(node_to_read_map.size());

    node_count = 0;
    barcode_to_node_id_unordered_map barcode_to_node_map;
    for (auto kv : node_to_read_map) {
        node_to_minimizers.push_back(move(kv.first.minimizers));
        barcode_to_node_map[kv.first.barcode].push_back(node_count);
        for (read_id_t rid : kv.second) {
            read_to_node_vector[rid] = node_count;
        }
        node_count++;
    }
    node_to_read_map.clear();

    barcode_count = barcode_to_node_map.size();
    barcodes.reserve(barcode_count);
    barcode_to_nodes_vector.reserve(barcode_count);

    for (auto kv : barcode_to_node_map) {
        barcodes.push_back(move(kv.first));
        sort(kv.second.begin(), kv.second.end());
        barcode_to_nodes_vector.push_back(move(kv.second));
    }
    barcode_to_node_map.clear();

    if (!silent) {
        cout << "Read count: " << read_count << "\n";
        cout << "Node count: " << node_count << "\n";
        cout << "Barcode count: " << barcode_count << "\n";
    }
    dog << "Read count: " << read_count << "\n";
    dog << "Node count: " << node_count << "\n";
    dog << "Barcode count: " << barcode_count << "\n";

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

    min_k_mer = (min_k_mer < current_k_mer)  ? min_k_mer : current_k_mer;

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

        min_k_mer = (min_k_mer < current_k_mer) ? min_k_mer : current_k_mer;
        // min_k_mer = min_k_mer < current_k_mer ? min_k_mer : current_k_mer;
        // printf("0x%08X\n", min_k_mer);

    }

    return min_k_mer != (minimizer_t) -1 ? min_k_mer : --invalid_kmer;
    // return min_k_mer;
}
