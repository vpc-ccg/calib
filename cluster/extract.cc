//
// Created by borabi on 19/12/17.
//

#include "extract.h"
#include "global.h"

using namespace std;


#define ASCII_SIZE 128
#define BYTE_SIZE 8
#define ENCODE_SIZE 2

node_id_t node_count = 0;
read_id_t read_count = 0;

// node_to_read_id find reads with identical barcodes and minimizers
typedef std::unordered_map<Node, std::vector<read_id_t>, NodeHash, NodeEqual> node_to_read_id_unordered_map;

read_vector reads;
node_vector nodes;
node_id_to_read_id_vector node_to_read_vector;

minimizer_t encode [ASCII_SIZE];
minimizer_t invalid_kmer = (minimizer_t) -1;


void extract_barcodes_and_minimizers() {
    node_to_read_id_unordered_map node_to_read_map;

    for (int i = 0; i < ASCII_SIZE; i++)
        encode[i] = (minimizer_t) -1;
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

    dog << "Reading fastq files\n";

    string trash;
    Node current_node;
    reads.push_back(Read());

    cout << "#nodes\t" << node_count << "\n";

    // Processing FASTQ files one read at a time
    while (getline(fastq1, reads.back().name_1)) {
        getline(fastq1, reads.back().sequence_1);
        getline(fastq1, trash);
        getline(fastq1, reads.back().quality_1);
        getline(fastq2, reads.back().name_2);
        getline(fastq2, reads.back().sequence_2);
        getline(fastq2, trash);
        getline(fastq2, reads.back().quality_2);

        int s1_length = reads.back().sequence_1.size();
        int s2_length = reads.back().sequence_2.size();

        // Extracting the barcode from the start of both mates
        if (s1_length >= barcode_length && s2_length >= barcode_length){
            current_node.barcode =
                reads.back().sequence_1.substr(0, barcode_length) +
                reads.back().sequence_2.substr(0, barcode_length);
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
                current_node.minimizers_1[i] = minimizer(reads.back().sequence_1, start, s1_seg_length);
                start += s1_seg_length;
            }
        } else {
            for (int i = 0; i < minimizer_count; i++){
                current_node.minimizers_1[i] = -1;
            }
        }

        int s2_seg_length = s2_length / minimizer_count;
        if (s2_seg_length >= kmer_size){
            int start = barcode_length;
            for (int i = 0; i < minimizer_count; i++){
                current_node.minimizers_2[i] = minimizer(reads.back().sequence_2, start, s2_seg_length);
                start += s2_seg_length;
            }
        } else {
            for (int i = 0; i < minimizer_count; i++){
                current_node.minimizers_2[i] = -1;
            }
        }


        if (node_to_read_map.find(current_node) != node_to_read_map.end()) {
            node_to_read_map[current_node].push_back(reads.size());
        } else {
            node_to_read_map.emplace(current_node, vector<node_id_t>{reads.size()});
            // vector<node_id_t> current_vector;
            // current_vector.push_back(read_count++);
            // node_to_read[current_node] = current_vector;
            current_node = Node();
        }
        reads.push_back(Read());

//
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
    cout << node_to_read_map.size() << "\t" << "\t" << node_count << "\n";

    nodes.reserve(node_count);
    node_to_read_vector.reserve(node_count);
    for (auto kv : node_to_read_map){
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
        if (encode[(size_t)seq.at(i)] == (minimizer_t) -1){
            current_k_mer = (minimizer_t) -1;
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
        if (encode[(size_t)seq.at(i)] == (minimizer_t) -1){
            // printf("HERE\n");
            current_k_mer = (minimizer_t) -1;
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
                if (encode[(size_t)seq.at(j)] == (minimizer_t) -1){
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
