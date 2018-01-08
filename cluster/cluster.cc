//
// Created by borabi on 20/12/17.
//

#include "cluster.h"
#include <igraph.h>

using namespace std;

char masked_barcode_buffer[150];

void cluster(){
    cout << "# of nodes: "<< node_count << "\n";
    cout << "# of nodes: "<< node_to_read_vector.size() << "\n";
    // for (node_id_t i = 0; i < node_count; i++){
    //     print_node(i);
    // }

    // adjacency_sets per nodes
    node_id_to_node_id_vector adjacency_sets(node_count);

    barcode_similarity(adjacency_sets);

    dog << "Removing edged of unmatched minimizers\n";
    cout << "Removing edged of unmatched minimizers\n";

    remove_edges_of_unmatched_minimizers(adjacency_sets);

    cout <<  adjacency_sets.size() << "\n";
    for (size_t i = 0; i < node_count; i++){
        dog << i << "\t" << adjacency_sets[i].size();
        cout << i << "\t" << adjacency_sets[i].size();
        for (auto neighbor : adjacency_sets[i]){
            dog << neighbor << "\t";
            cout << neighbor << "\t";
        }
        dog << "\n";
        cout << "\n";
    }
    std::cout << "HERE WE ARE!" << '\n';
}

bool unmatched_minimimizers(node_id_t node_id, node_id_t neighbor_id){
    int unmatched_minimimizers_1 = 0;
    int unmatched_minimimizers_2 = 0;


    std::cout << node_id << " ---- " << neighbor_id << '\n';
    for (int i =0; i < minimizer_count; i++){
        std::cout << nodes[node_id].minimizers_1[i] << " ---- " << nodes[neighbor_id].minimizers_1[i] << '\n';
        std::cout << nodes[node_id].minimizers_2[i] << " ---- " << nodes[neighbor_id].minimizers_2[i] << '\n';
        unmatched_minimimizers_1 += nodes[node_id].minimizers_1[i] != nodes[neighbor_id].minimizers_1[i];
        unmatched_minimimizers_2 += nodes[node_id].minimizers_2[i] != nodes[neighbor_id].minimizers_2[i];
    }
    cout << unmatched_minimimizers_1 + unmatched_minimimizers_2 << "\n";
    return unmatched_minimimizers_1 > minimizer_threshold && unmatched_minimimizers_2 > minimizer_threshold;
}



void remove_edges_of_unmatched_minimizers(node_id_to_node_id_vector &adjacency_sets){
    int removed_count = 0;
    for (node_id_t node = 0; node < adjacency_sets.size(); node++){
        adjacency_sets[node].erase(node);
        std::cout << "NODE: " << node << "\twith neighborhood of\t" << adjacency_sets[node].size()<< '\n';

        for (auto neighbor = adjacency_sets[node].begin(); neighbor != adjacency_sets[node].end(); ){
            std::cout << "==>" << *neighbor << '\n';
            if (unmatched_minimimizers(node, *neighbor)){
                removed_count++;
                adjacency_sets[*neighbor].erase(node);
                neighbor = adjacency_sets[node].erase(neighbor);
            } else {
                neighbor++;
            }
        }
    }
    cout << "Removed " << removed_count << " edges\n";
}

void barcode_similarity(node_id_to_node_id_vector &adjacency_sets){
    masked_barcode_to_node_id_unordered_map lsh;

    vector<bool> mask(barcode_length);
    std::fill(mask.begin() + error_tolerance, mask.end(), true);
    masked_barcode_buffer[barcode_length-error_tolerance] = '\0';

    do {
        for (char c = 'A'; c < 'A' + barcode_length-error_tolerance + 1; c++){
            masked_barcode_buffer[c-'A'] = c;
        }
        cout << mask_barcode(string(masked_barcode_buffer), mask) << "\n";
        for (node_id_t i = 0; i < node_count; i++){
            lsh[mask_barcode(nodes[i].barcode, mask)].push_back(i);
        }
    } while (std::next_permutation(mask.begin(), mask.end()));

    for (auto bucket : lsh){
        // cout << bucket.first << "\twith size\t" << bucket.second.size() << "\n";

        for (node_id_t node : bucket.second){

            // cout << "==> " << node << "\n";

            adjacency_sets[node].insert(bucket.second.begin(), bucket.second.end());


            // for (auto n : adjacency_sets[node]){
                // cout << "====> " << n << "\n";
            // }
        }
    }
    // cout << "adjacency_sets.size() : "  << adjacency_sets.size() << "\n";

}


string mask_barcode(const string& barcode, const vector<bool>& mask){
    int pos = 0;
    for (int i = 0; i < barcode_length; i++){
        if (mask[i]){
            masked_barcode_buffer[pos] = barcode.at(i);
            pos++;
        }
    }
    return string(masked_barcode_buffer);
}


void print_node(node_id_t node_id){
    dog <<  node_id << "\t" << nodes[node_id].barcode << "\t" ;
    for (int i =0; i < minimizer_count; i++)
        dog << nodes[node_id].minimizers_1[i] << "\t";
    for (int i =0; i < minimizer_count; i++)
        dog << nodes[node_id].minimizers_2[i] << "\t";
    dog << "\n";

    for (read_id_t read: node_to_read_vector[node_id])
        dog << "\t" << read << "\t" << reads[read].sequence_1 << "\t" << reads[read].sequence_2 << "\n";

}
