//
// Created by borabi on 20/12/17.
//

#include "cluster.h"
#include <igraph.h>

using namespace std;


void cluster(){
    cout << "# of nodes: "<< node_count << "\n";
    cout << "# of nodes: "<< node_to_read_vector.size() << "\n";
    for (node_id_t i = 0; i < node_count; i++){
        print_node(i);
    }
    cout << "HERE!\n";

    // adjacency_sets per nodes
    node_id_to_node_id_vector adjacency_sets(node_count);
    cout << "HERE!\n";

    for (int i = 0; i < 3; i++){
        dog << "LSH of id " << i << "\n";
        cout << "Mask ID " << i << "\n";

        barcode_similarity(i, adjacency_sets);
    }

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
        std::cout << "NODE: " << node << "\twith neighborhood of\t" << adjacency_sets[node].size()<< '\n';
        for (node_id_t neighbor : adjacency_sets[node]){
            std::cout << "==>" << neighbor << '\n';
            if (unmatched_minimimizers(node, neighbor)){
                removed_count++;
                adjacency_sets[node].erase(neighbor);
                adjacency_sets[neighbor].erase(node);
            }
        }
    }
    cout << "Removed " << removed_count << " edges\n";
}

void barcode_similarity(int mask_id, node_id_to_node_id_vector &adjacency_sets){
    cout << "Processing mask ID " << mask_id << " for " << adjacency_sets.size() << " sets\n";
    masked_barcode_to_node_id_unordered_map lsh;
    for (node_id_t i = 0; i < node_count; i++){
        lsh[mask(nodes[i].barcode, mask_id)].push_back(i);
    }

    for (auto bucket : lsh){
        cout << bucket.first << "\twith size\t" << bucket.second.size() << "\n";

        for (node_id_t node : bucket.second){

            cout << "==> " << node << "\n";

            adjacency_sets[node].insert(bucket.second.begin(), bucket.second.end());


            for (auto n : adjacency_sets[node]){
                cout << "====> " << n << "\n";
            }
        }
    }
    cout << "adjacency_sets.size() : "  << adjacency_sets.size() << "\n";

}

string mask(string barcode, int mask_id){
    return barcode;
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
