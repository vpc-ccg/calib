//
// Created by borabi on 20/12/17.
//

#include "cluster.h"

#include <stack>
#include <algorithm>
#include <time.h>

using namespace std;

char masked_barcode_buffer[150];


void cluster(){
    node_id_to_node_id_vector_of_vectors adjacency_lists(node_count);
    time_t start;

    cout << "Adding edges due to barcode barcode similarity\n";
    start = time(NULL);
    barcode_similarity(adjacency_lists);
    cout << "Edges due to barcodes took: " << difftime(time(NULL), start) << "\n";
    dog << "Edges due to barcodes took: " << difftime(time(NULL), start) << "\n";

    cout << "Removing edges of unmatched minimizers\n";
    start = time(NULL);
    remove_edges_of_unmatched_minimizers(adjacency_lists);
    cout << "Removing edges due to minimizers: " << difftime(time(NULL), start)<< "\n";
    dog << "Removing edges due to minimizers: " << difftime(time(NULL), start)<< "\n";

    cout << "Extracting and outputting clusters\n";
    start = time(NULL);
    extract_clusters(adjacency_lists);
    cout << "Extracting clusters and outputting took: " << difftime(time(NULL), start) << "\n";
    dog << "Extracting clusters and outputting took: " << difftime(time(NULL), start) << "\n";

}

void barcode_similarity(node_id_to_node_id_vector_of_vectors &adjacency_lists){
    vector<bool> mask(barcode_length*2, false);
    std::fill(mask.begin() + error_tolerance, mask.end(), true);
    masked_barcode_buffer[barcode_length*2-error_tolerance] = '\0';

    string template_barcode;
    for (char c = 'A'; c < 'A' + barcode_length*2; c++) {
        template_barcode += c;
    }
    time_t start;
    time_t build_time = 0, process_time = 0;
    do {
        start = time(NULL);
        masked_barcode_to_node_id_unordered_map lsh;
        cout << mask_barcode(string(template_barcode), mask) << "\n";
        dog << mask_barcode(string(template_barcode), mask) << "\n";
        for (node_id_t i = 0; i < node_count; i++) {
            lsh[mask_barcode(nodes[i].barcode, mask)].push_back(i);
        }
        build_time += difftime(time(NULL), start);
        cout << "Building LSH took: " << difftime(time(NULL), start) << "\n";
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
        cout << "Processing LSH took: " << difftime(time(NULL), start) << "\n";
        dog << "Processing LSH took: " << difftime(time(NULL), start) << "\n";
    } while (std::next_permutation(mask.begin(), mask.end()));
    cout << "Building all LSH took: " << build_time << "\n";
    dog << "Building all LSH took: " << build_time << "\n";
    cout << "Processing all LSH took: " << process_time << "\n";
    dog << "Processing all LSH took: " << process_time << "\n";
}


string mask_barcode(const string& barcode, const vector<bool>& mask){
    int pos = 0;
    for (int i = 0; i < barcode_length*2; i++) {
        if (mask[i]) {
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
            if (!unmatched_minimimizers(node, neighbor)) {
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
    return !(matched_minimimizers_1 >= minimizer_threshold && matched_minimimizers_2 >= minimizer_threshold);
}


void process_lsh(masked_barcode_to_node_id_unordered_map &lsh,
                 node_id_to_node_id_vector_of_vectors adjacency_lists,
                 size_t reminder,
                 size_t divisor){


}

void extract_clusters(node_id_to_node_id_vector_of_vectors &adjacency_lists){
    vector<bool> pushed(node_count, false);
    stack<node_id_t> opened;
    size_t cluster_count = 0;
    ofstream clusters;
    clusters = ofstream(output_prefix + ".cluster");

    for (node_id_t node = 0; node < node_count; node++) {
        if (!pushed[node]) {
            // cout << "cluster # " << cluster_count << ":\t";
            clusters << "# "<< cluster_count <<"\n";
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
                    clusters << opened.top() << "\t" << read << "\t";
                    clusters << reads[read].name_1 << "\t" << reads[read].sequence_1 << "\t" << reads[read].quality_1 <<
                        "\t";
                    clusters << reads[read].name_2 << "\t" << reads[read].sequence_2 << "\t" << reads[read].quality_2 <<
                        "\n";
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
