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
    cout << "# of nodes: "<< node_count << "\n";
    cout << "# of nodes: "<< node_to_read_vector.size() << "\n";
    // for (node_id_t i = 0; i < node_count; i++){
    //     print_node(i);
    // }

    // adjacency_sets per nodes
    node_id_to_node_id_vector_of_vectors adjacency_lists(node_count);

    time_t start = time(NULL);
    barcode_similarity(adjacency_lists);
    cout << "Edges due to barcodes took: " << -difftime(start, time(NULL)) << "\n";

    dog << "Removing edges of unmatched minimizers\n";
    cout << "Removing edges of unmatched minimizers\n";

    start = time(NULL);
    remove_edges_of_unmatched_minimizers(adjacency_lists);
    cout << "Removing edges due to minimizers: " << -difftime(start, time(NULL)) << "\n";

    dog << "Extracting clusters\n";
    cout << "Extracting clusters\n";

    start = time(NULL);
    extract_clusters(adjacency_lists);
    cout << "Extracting clusters + output took: " << -difftime(start, time(NULL)) << "\n";

    // cout <<  adjacency_sets.size() << "\n";
    // for (size_t i = 0; i < node_count; i++){
    //     dog << i << "\t" << adjacency_sets[i].size();
    //     cout << i << "\t" << adjacency_sets[i].size();
    //     for (auto neighbor : adjacency_sets[i]){
    //         dog << neighbor << "\t";
    //         cout << neighbor << "\t";
    //     }
    //     dog << "\n";
    //     cout << "\n";
    // }
}

void extract_clusters(node_id_to_node_id_vector_of_vectors &adjacency_lists){
    vector<bool> pushed(node_count, false);
    stack<node_id_t> opened;
    size_t cluster_count = 0;
    for (node_id_t node = 0; node < node_count; node++) {
        if (!pushed[node]) {
            // cout << "cluster # " << cluster_count << ":\t";
            dog << "# "<< cluster_count <<"\n";
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
                    dog << opened.top() << "\t" << read << "\t";
                    dog << reads[read].name_1 << "\t" << reads[read].sequence_1 << "\t" << reads[read].quality_1 <<
                        "\t";
                    dog << reads[read].name_2 << "\t" << reads[read].sequence_2 << "\t" << reads[read].quality_2 <<
                        "\n";
                }
                opened.pop();
            }
            cluster_count++;
        }
    }
}

bool unmatched_minimimizers(node_id_t node_id, node_id_t neighbor_id){
    int matched_minimimizers_1 = 0;
    int matched_minimimizers_2 = 0;
    // std::cout << node_id << " ---- " << neighbor_id << '\n';
    for (int i =0; i < minimizer_count; i++) {
        // std::cout << nodes[node_id].minimizers_1[i] << " ---- " << nodes[neighbor_id].minimizers_1[i] << '\n';
        // std::cout << nodes[node_id].minimizers_2[i] << " ---- " << nodes[neighbor_id].minimizers_2[i] << '\n';
        matched_minimimizers_1 += nodes[node_id].minimizers_1[i] == nodes[neighbor_id].minimizers_1[i];
        matched_minimimizers_2 += nodes[node_id].minimizers_2[i] == nodes[neighbor_id].minimizers_2[i];
    }
    // cout << unmatched_minimimizers_1 + unmatched_minimimizers_2 << "\n";
    return !(matched_minimimizers_1 >= minimizer_threshold && matched_minimimizers_2 >= minimizer_threshold);
}

void remove_edges_of_unmatched_minimizers(node_id_to_node_id_vector_of_vectors &adjacency_lists){
    // int removed_count = 0;
    for (node_id_t node = 0; node < node_count; node++) {
        // std::cout << "NODE: " << node << "\twith neighborhood of\t" << adjacency_sets[node].size()<< '\n';
        vector<node_id_t> good_neighbors;
        for (node_id_t neighbor : adjacency_lists[node]) {
            // std::cout << "==>" << *neighbor << '\n';
            if (!unmatched_minimimizers(node, neighbor)) {
                good_neighbors.push_back(neighbor);
            }
        }
        adjacency_lists[node] = move(good_neighbors);
    }
    // cout << "Removed " << removed_count << " edges\n";
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
        for (node_id_t i = 0; i < node_count; i++) {
            lsh[mask_barcode(nodes[i].barcode, mask)].push_back(i);
        }
        build_time += difftime(time(NULL), start);
        cout << "Building LSH took: " << difftime(time(NULL), start) << "\n";

        //int buckets_processed = 0;
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
            //buckets_processed++;
        }
//        lsh.clear();
        process_time += difftime(time(NULL), start);
        cout << "Processing LSH took: " << difftime(time(NULL), start) << "\n";
    } while (std::next_permutation(mask.begin(), mask.end()));
    cout << "Building all LSH took: " << build_time << "\n";
    cout << "Processing all LSH took: " << process_time << "\n";
    // for (auto neighborhood : adjacency_lists){
    //   adjacency_sets.push_back(  unordered_set<node_id_t>(make_move_iterator(neighborhood.begin()), make_move_iterator(neighborhood.end()) )    );
    // }

}

void process_lsh(masked_barcode_to_node_id_unordered_map &lsh,
                 node_id_to_node_id_vector_of_vectors adjacency_lists,
                 size_t reminder,
                 size_t divisor){


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
