//
// Created by borabi on 20/12/17.
//

#include "cluster.h"

#include <stack>
#include <algorithm>
#include <time.h>
// Debug includes
#include <sstream>
#include <iomanip>

using namespace std;

char masked_barcode_buffer[150];
#define ASCII_SIZE 256
bool valid_base [ASCII_SIZE];



void cluster(){
    node_id_to_node_id_vector_of_vectors adjacency_lists(node_count);
    time_t start;

    if (!silent) {
        cout << "Adding edges due to barcode barcode similarity\n";
    }
    start = time(NULL);
    barcode_similarity(adjacency_lists);
    if (!silent) {
        cout << "Edges due to barcodes took: " << difftime(time(NULL), start) << "\n";
    }
    dog << "Edges due to barcodes took: " << difftime(time(NULL), start) << "\n";

    if (!silent) {
        cout << "Removing edges of unmatched minimizers\n";
    }
    start = time(NULL);
    remove_edges_of_unmatched_minimizers(adjacency_lists);
    dog << "Removing edges due to minimizers: " << difftime(time(NULL), start)<< "\n";
    if (!silent) {
        cout << "Removing edges due to minimizers: " << difftime(time(NULL), start)<< "\n";
    }

    if (!silent) {
        cout << "Extracting and outputting clusters\n";
    }
    start = time(NULL);
    extract_clusters(adjacency_lists);
    if (!silent) {
        cout << "Extracting clusters and outputting took: " << difftime(time(NULL), start) << "\n";
    }
    dog << "Extracting clusters and outputting took: " << difftime(time(NULL), start) << "\n";

}

void barcode_similarity(node_id_to_node_id_vector_of_vectors &adjacency_lists){
    for (int i = 0; i < ASCII_SIZE; i++) {
        valid_base[i] = false;
    }
    valid_base['A'] = true;
    valid_base['C'] = true;
    valid_base['G'] = true;
    valid_base['T'] = true;
    valid_base['a'] = true;
    valid_base['c'] = true;
    valid_base['g'] = true;
    valid_base['t'] = true;

    vector<bool> mask(barcode_length*2, false);
    std::fill(mask.begin() + error_tolerance, mask.end(), true);
    masked_barcode_buffer[barcode_length*2-error_tolerance] = '\0';

    string template_barcode;
    for (char c = 'A'; c < 'A' + barcode_length*2; c++) {
        template_barcode += c;
    }
    time_t start;
    time_t build_time = 0, process_time = 0;
    // int template_id
    do {
        start = time(NULL);
        masked_barcode_to_node_id_unordered_map lsh;
        if (!silent) {
            string current_mask_bin;
            for (bool p: mask) {
                current_mask_bin += p ? "1" : "0";
            }
            cout << current_mask_bin << "\n";
        }
        dog << mask_barcode(string(template_barcode), mask) << "\n";
        string masked_barcode;
        for (node_id_t i = 0; i < node_count; i++) {
            masked_barcode = mask_barcode(nodes[i].barcode, mask);
            if (masked_barcode != "0"){
                lsh[masked_barcode].push_back(i);
            }
        }
        build_time += difftime(time(NULL), start);
        if (!silent) {
            cout << "Building LSH took: " << difftime(time(NULL), start) << "\n";
        }
        dog << "Building LSH took: " << difftime(time(NULL), start) << "\n";
        start = time(NULL);
        for (auto bucket : lsh) {
            sort(bucket.second.begin(), bucket.second.end());
            for (node_id_t node : bucket.second) {
                vector<node_id_t> result;
                set_union(adjacency_lists[node].begin(), adjacency_lists[node].end(),
                          bucket.second.begin(), bucket.second.end(),
                          back_inserter(result)
                          );
                adjacency_lists[node] = move(result);
            }
        }
        process_time += difftime(time(NULL), start);
        if (!silent) {
            cout << "Processing LSH took: " << difftime(time(NULL), start) << "\n";
        }
        dog << "Processing LSH took: " << difftime(time(NULL), start) << "\n";
    } while (std::next_permutation(mask.begin(), mask.end()));
    if (!silent) {
        cout << "Building all LSH took: " << build_time << "\n";
    }
    dog << "Building all LSH took: " << build_time << "\n";
    if (!silent) {
        cout << "Processing all LSH took: " << process_time << "\n";
    }
    dog << "Processing all LSH took: " << process_time << "\n";
}


string mask_barcode(const string& barcode, const vector<bool>& mask){
    int pos = 0;
    for (int i = 0; i < barcode_length*2; i++) {
        if (mask[i]) {
            if (valid_base[(uint8_t) barcode.at(i)] == false) {
                return "0";
            }
            masked_barcode_buffer[pos] = barcode.at(i);
            pos++;
        }
    }
    return string(masked_barcode_buffer);
}

void remove_edges_of_unmatched_minimizers(node_id_to_node_id_vector_of_vectors &adjacency_lists){
    for (node_id_t node = 0; node < node_count; node++) {
        vector<node_id_t> good_neighbors;
        for (node_id_t neighbor : adjacency_lists[node]) {
            if (node != neighbor && !unmatched_minimimizers(node, neighbor)) {
                good_neighbors.push_back(neighbor);
            }
        }
        adjacency_lists[node] = move(good_neighbors);
    }
}


bool unmatched_minimimizers(node_id_t node_id, node_id_t neighbor_id){
    int matched_minimimizers_1 = 0;
    int matched_minimimizers_2 = 0;
    for (int i =0; i < minimizer_count; i++) {
        matched_minimimizers_1 += nodes[node_id].minimizers_1[i] == nodes[neighbor_id].minimizers_1[i];
        matched_minimimizers_2 += nodes[node_id].minimizers_2[i] == nodes[neighbor_id].minimizers_2[i];
    }
    if (debug) {
        int hamming_distance = 0;
        for (int i = 0; i < barcode_length*2; i++) {
            hamming_distance += nodes[node_id].barcode.at(i) != nodes[neighbor_id].barcode.at(i);
        }
        dog << "M\t";
        dog << !(matched_minimimizers_1 >= minimizer_threshold && matched_minimimizers_2 >= minimizer_threshold) << "\t";
        dog << matched_minimimizers_1 << "\t" << matched_minimimizers_2 << "\t" << hamming_distance <<"\n";

        dog << "M1\t";
        for (int i =0; i < minimizer_count; i++) {
            dog << nodes[node_id].minimizers_1[i] << "\t";
        }
        for (int i =0; i < minimizer_count; i++) {
            dog << nodes[node_id].minimizers_2[i] << "\t";
        }
        dog << nodes[node_id].barcode << "\t";
        dog << reads[node_to_read_vector[node_id].front()].sequence_1 << "\t" << reads[node_to_read_vector[node_id].front()].sequence_2 <<"\n";

        dog << "M2\t";
        for (int i =0; i < minimizer_count; i++) {
            dog << nodes[neighbor_id].minimizers_1[i] << "\t";
        }
        for (int i =0; i < minimizer_count; i++) {
            dog << nodes[neighbor_id].minimizers_2[i] << "\t";
        }
        dog << nodes[neighbor_id].barcode << "\t";
        dog << reads[node_to_read_vector[neighbor_id].front()].sequence_1 << "\t" << reads[node_to_read_vector[neighbor_id].front()].sequence_2 <<"\n";

    }
    return !(matched_minimimizers_1 >= minimizer_threshold && matched_minimimizers_2 >= minimizer_threshold);
}


// void process_lsh(masked_barcode_to_node_id_unordered_map &lsh,
//                  node_id_to_node_id_vector_of_vectors adjacency_lists,
//                  size_t reminder,
//                  size_t divisor){
//
//
// }


void extract_clusters(node_id_to_node_id_vector_of_vectors &adjacency_lists){
    vector<bool> pushed(node_count, false);
    stack<node_id_t> opened;
    size_t cluster_count = 0;
    ofstream clusters;
    if (bc_format) {
        clusters = ofstream(output_prefix + "bc");
    } else {
        clusters = ofstream(output_prefix + "cluster");
    }

    ofstream cluster_debug;
    ofstream cluster_R1;
    ofstream cluster_R2;
    if (debug) {
        cluster_debug = ofstream(output_prefix + "cluster.debug");
        cluster_R1 = ofstream(output_prefix + "cluster.R1.fastq");
        cluster_R2 = ofstream(output_prefix + "cluster.R2.fastq");
    }
    for (node_id_t node = 0; node < node_count; node++) {
        if (!pushed[node]) {
            if (bc_format) {

            } else {
                clusters << "# "<< cluster_count <<"\n";
            }
            if (debug) {
                cluster_debug << "#\t" << cluster_count << "\n";
            }
            opened.push(node);
            pushed[node] = true;
            while(!opened.empty()) {
                for (node_id_t neighbor: adjacency_lists[opened.top()]) {
                    if (!pushed[neighbor]) {
                        opened.push(neighbor);
                        pushed[neighbor] = true;
                    }
                }
                for (read_id_t read : node_to_read_vector[opened.top()]) {
                    if (bc_format) {
                        clusters << cluster_count << "\t" << reads[read].sequence_1 << "\t" << reads[read].sequence_2 << "\n";
                    } else {
                        clusters << opened.top() << "\t" << read << "\t";
                        clusters << reads[read].name_1 << "\t" << reads[read].sequence_1 << "\t" << reads[read].quality_1 <<
                        "\t";
                        clusters << reads[read].name_2 << "\t" << reads[read].sequence_2 << "\t" << reads[read].quality_2 <<
                        "\n";
                    }
                }
                if (debug) {
                    node_id_t current_node = opened.top();
                    double left_matching_minimizers  = 0.0;
                    double right_matching_minimizers = 0.0;
                    double hamming_distance = 0.0;

                    stringstream stream_neighbors;

                    for (node_id_t neighbor: adjacency_lists[current_node]) {
                        double current_left_matching_minimizers = 0.0;
                        double current_right_matching_minimizers = 0.0;
                        double current_hamming_distance = 0.0;

                        stream_neighbors << neighbor << "_";

                        for (int i = 0; i < minimizer_count; i++) {
                            current_left_matching_minimizers  += nodes[current_node].minimizers_1[i] == nodes[neighbor].minimizers_1[i];
                            current_right_matching_minimizers += nodes[current_node].minimizers_2[i] == nodes[neighbor].minimizers_2[i];
                        }
                        stream_neighbors << current_left_matching_minimizers << "_" << current_right_matching_minimizers << "_";

                        for (int i = 0; i < barcode_length*2; i++) {
                            current_hamming_distance += nodes[current_node].barcode.at(i) != nodes[neighbor].barcode.at(i);
                        }
                        stream_neighbors << current_hamming_distance << ", ";

                        hamming_distance += current_hamming_distance;
                        left_matching_minimizers  += current_left_matching_minimizers;
                        right_matching_minimizers += current_right_matching_minimizers;

                    }

                    double average_left_connectivity;
                    double average_right_connectivity;
                    double average_hamming;

                    if (adjacency_lists[current_node].size() == 0) {
                        average_left_connectivity  = 0.0;
                        average_right_connectivity = 0.0;
                        average_hamming = 0.0;
                    } else {
                        average_left_connectivity  =  left_matching_minimizers/adjacency_lists[current_node].size();
                        average_right_connectivity = right_matching_minimizers/adjacency_lists[current_node].size();
                        average_hamming = hamming_distance/adjacency_lists[current_node].size();
                    }
                    stringstream stream;
                    stream.precision(4);
                    stream << fixed;

                    stream << cluster_count << "\t";
                    stream << std::setfill (' ') << std::setw (10) << current_node << "\t";
                    stream << std::setfill (' ') << std::setw (3) << adjacency_lists[current_node].size() << "\t";
                    stream << std::setfill (' ') << std::setw (2) << node_to_read_vector[current_node].size() << "\t";

                    stream << average_left_connectivity  << "\t";
                    stream << average_right_connectivity << "\t";
                    stream << average_hamming << "\t";
                    stream << stream_neighbors.str();

                    cluster_debug << stream.str() << "\n";

                    for (read_id_t read : node_to_read_vector[current_node]) {
                        cluster_R1 << "@" << cluster_count << "_" << current_node << "_" << read << "/1\n";
                        cluster_R1 << reads[read].sequence_1 << "\n";
                        cluster_R1 << "+\n";
                        cluster_R1 << string(reads[read].sequence_1.size(), 'K') << "\n";

                        cluster_R2 << "@" << cluster_count << "_" << current_node << "_" << read << "/2\n";
                        cluster_R2 << reads[read].sequence_2 << "\n";
                        cluster_R2 << "+\n";
                        cluster_R2 << string(reads[read].sequence_2.size(), 'K') << "\n";
                    }

                }
                opened.pop();
            }
            cluster_count++;
        }
    }
    clusters.close();
}

void print_node(node_id_t node_id){
    dog <<  node_id << "\t" << nodes[node_id].barcode << "\t";
    for (int i =0; i < minimizer_count; i++)
        dog << nodes[node_id].minimizers_1[i] << "\t";
    for (int i =0; i < minimizer_count; i++)
        dog << nodes[node_id].minimizers_2[i] << "\t";
    dog << "\n";

    for (read_id_t read: node_to_read_vector[node_id])
        dog << "\t" << read << "\t" << reads[read].sequence_1 << "\t" << reads[read].sequence_2 << "\n";
}
