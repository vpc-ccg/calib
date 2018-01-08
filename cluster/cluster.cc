//
// Created by borabi on 20/12/17.
//

#include "cluster.h"

#include <stack>
#include <algorithm>

using namespace std;

char masked_barcode_buffer[150];

void cluster(){
    cout << "# of nodes: "<< node_count << "\n";
    cout << "# of nodes: "<< node_to_read_vector.size() << "\n";
    // for (node_id_t i = 0; i < node_count; i++){
    //     print_node(i);
    // }

    // adjacency_sets per nodes
    node_id_to_node_id_vector_of_sets adjacency_sets;

    barcode_similarity(adjacency_sets);

    dog << "Removing edges of unmatched minimizers\n";
    cout << "Removing edges of unmatched minimizers\n";

    remove_edges_of_unmatched_minimizers(adjacency_sets);

    dog << "Extracting clusters\n";
    cout << "Extracting clusters\n";

    extract_clusters(adjacency_sets);
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

void extract_clusters(node_id_to_node_id_vector_of_sets &adjacency_sets){
    vector<bool> pushed(node_count, false);
    stack<node_id_t> opened;
    size_t cluster_count = 0;
    for (node_id_t node = 0; node < node_count; node++){
        if (!pushed[node]){
            // cout << "cluster # " << cluster_count << ":\t";
            dog << "# "<< cluster_count <<"\n";
            opened.push(node);
            pushed[node] = true;
            while(!opened.empty()){
                for (node_id_t neighbor: adjacency_sets[opened.top()]){
                    if (!pushed[neighbor]){
                        opened.push(neighbor);
                        pushed[neighbor] = true;
                    }
                }
                for (read_id_t read : node_to_read_vector[opened.top()]){
                    dog << opened.top() << "\t" << read << "\t";
                    dog << reads[read].name_1 << "\t" << reads[read].sequence_1 << "\t" << reads[read].quality_1 << "\t";
                    dog << reads[read].name_2 << "\t" << reads[read].sequence_2 << "\t" << reads[read].quality_2 << "\n";
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
    for (int i =0; i < minimizer_count; i++){
        // std::cout << nodes[node_id].minimizers_1[i] << " ---- " << nodes[neighbor_id].minimizers_1[i] << '\n';
        // std::cout << nodes[node_id].minimizers_2[i] << " ---- " << nodes[neighbor_id].minimizers_2[i] << '\n';
        matched_minimimizers_1 += nodes[node_id].minimizers_1[i] == nodes[neighbor_id].minimizers_1[i];
        matched_minimimizers_2 += nodes[node_id].minimizers_2[i] == nodes[neighbor_id].minimizers_2[i];
    }
    // cout << unmatched_minimimizers_1 + unmatched_minimimizers_2 << "\n";
    return !(matched_minimimizers_1 >= minimizer_threshold && matched_minimimizers_2 >= minimizer_threshold);
}



void remove_edges_of_unmatched_minimizers(node_id_to_node_id_vector_of_sets &adjacency_sets){
    int removed_count = 0;
    for (node_id_t node = 0; node < adjacency_sets.size(); node++){
        adjacency_sets[node].erase(node);
        // std::cout << "NODE: " << node << "\twith neighborhood of\t" << adjacency_sets[node].size()<< '\n';

        for (auto neighbor = adjacency_sets[node].begin(); neighbor != adjacency_sets[node].end(); ){
            // std::cout << "==>" << *neighbor << '\n';
            if (unmatched_minimimizers(node, *neighbor)){
                removed_count++;
                adjacency_sets[*neighbor].erase(node);
                neighbor = adjacency_sets[node].erase(neighbor);
            } else {
                neighbor++;
            }
        }
    }
    // cout << "Removed " << removed_count << " edges\n";
}

void barcode_similarity(node_id_to_node_id_vector_of_sets &adjacency_sets){
    node_id_to_node_id_vector_of_vectors adjacency_lists(node_count);

    vector<bool> mask(barcode_length, false);
    std::fill(mask.begin() + error_tolerance, mask.end(), true);
    masked_barcode_buffer[barcode_length-error_tolerance] = '\0';

    string template_barcode;
    for (char c = 'A'; c < 'A' + barcode_length; c++){
        template_barcode += c;
    }

    do {
        cout << mask_barcode(string(template_barcode), mask) << "\n";
        masked_barcode_to_node_id_unordered_map lsh;
        for (node_id_t i = 0; i < node_count; i++){
            lsh[mask_barcode(nodes[i].barcode, mask)].push_back(i);
        }
        int buckets_processed = 0;
        for (auto bucket : lsh){
          sort(bucket.second.begin(), bucket.second.end());
          for (node_id_t node : bucket.second){
            vector<node_id_t> result;
            set_union(adjacency_lists[node].begin(), adjacency_lists[node].end(),
                        bucket.second.begin(), bucket.second.end(),
                        back_inserter(result)
                      );
            adjacency_lists[node] = move(result);

          }
          buckets_processed++;
        }
    } while (std::next_permutation(mask.begin(), mask.end()));

    for (auto neighborhood : adjacency_lists){
      adjacency_sets.push_back(  unordered_set<node_id_t>(make_move_iterator(neighborhood.begin()), make_move_iterator(neighborhood.end()) )    );
    }

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
