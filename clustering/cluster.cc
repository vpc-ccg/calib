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
        masked_barcode_to_barcode_id_unordered_map lsh;
        if (!silent) {
            string current_mask_bin;
            for (bool p: mask) {
                current_mask_bin += p ? "1" : "0";
            }
            cout << current_mask_bin << "\n";
        }
        dog << mask_barcode(string(template_barcode), mask) << "\n";
        string masked_barcode;
        for (barcode_id_t i = 0; i < barcode_count; i++) {
            masked_barcode = mask_barcode(barcodes[i], mask);
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
            for (barcode_id_t bid: bucket.second){
                for (node_id_t node : barcode_to_nodes_vector[bid]) {
                    for (barcode_id_t bid_o: bucket.second){
                        if (bid == bid_o){
                            continue;
                        }
                        vector<node_id_t> result;
                        set_union(adjacency_lists[node].begin(), adjacency_lists[node].end(),
                                  barcode_to_nodes_vector[bid_o].begin(), barcode_to_nodes_vector[bid_o].end(),
                                  back_inserter(result)
                                  );
                        adjacency_lists[node] = move(result);
                    }
                }
            }
        }
        process_time += difftime(time(NULL), start);
        if (!silent) {
            cout << "Processing LSH took: " << difftime(time(NULL), start) << "\n";
        }
        dog << "Processing LSH took: " << difftime(time(NULL), start) << "\n";
    } while (std::next_permutation(mask.begin(), mask.end()));

    for (barcode_id_t i = 0; i < barcode_count; i++) {
        for (node_id_t node : barcode_to_nodes_vector[i]) {
            vector<node_id_t> result;
            set_union(adjacency_lists[node].begin(), adjacency_lists[node].end(),
                      barcode_to_nodes_vector[i].begin(), barcode_to_nodes_vector[i].end(),
                      back_inserter(result)
                      );
            adjacency_lists[node] = move(result);
        }
    }
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
        matched_minimimizers_1 += node_to_minimizers[node_id][i] == node_to_minimizers[neighbor_id][i];
        matched_minimimizers_2 += node_to_minimizers[node_id][i+minimizer_count] == node_to_minimizers[neighbor_id][i+minimizer_count];
    }
    return !(matched_minimimizers_1 >= minimizer_threshold && matched_minimimizers_2 >= minimizer_threshold);
}


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

    ofstream cluster_R1;
    ofstream cluster_R2;
    for (node_id_t node = 0; node < node_count; node++) {
        if (!pushed[node]) {
            if (bc_format) {

            } else {
                clusters << "# "<< cluster_count <<"\n";
            }
            opened.push(node);
            pushed[node] = true;
            while(!opened.empty()) {
                node_id_t current_node = opened.top();
                opened.pop();
                for (node_id_t neighbor: adjacency_lists[current_node]) {
                    if (!pushed[neighbor]) {
                        opened.push(neighbor);
                        pushed[neighbor] = true;
                    }
                }
                for (read_id_t read : node_to_reads_vector[current_node]) {
                    if (bc_format) {
                        clusters << cluster_count << "\t" << reads[read].sequence_1 << "\t" << reads[read].sequence_2 << "\n";
                    } else {
                        clusters << current_node << "\t" << read << "\t";
                        clusters << reads[read].name_1 << "\t" << reads[read].sequence_1 << "\t" << reads[read].quality_1 <<
                        "\t";
                        clusters << reads[read].name_2 << "\t" << reads[read].sequence_2 << "\t" << reads[read].quality_2 <<
                        "\n";
                    }
                }
            }
            cluster_count++;
        }
    }
    clusters.close();
}
