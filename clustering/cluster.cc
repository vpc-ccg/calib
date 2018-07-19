//
// Created by borabi on 20/12/17.
//

#include "cluster.h"

// #include <pthread.h>
#include <thread>
#include <mutex>
#include <stack>
#include <algorithm>
#include <time.h>
// Debug includes
#include <sstream>
#include <iomanip>

using namespace std;


// extern variables declarations
cluster_id_t cluster_count = 0;
node_id_to_cluster_id_vector node_to_cluster_vector;

// locally global variables
#define ASCII_SIZE 256
bool valid_base [ASCII_SIZE];
node_id_to_node_id_vector_of_vectors* graph_ptr;
mutex graph_lock;
vector<vector<bool> > all_masks;

void cluster(){
    time_t start;
    graph_ptr = new node_id_to_node_id_vector_of_vectors();
    if (!silent) {
        cout << "Adding edges due to barcode barcode similarity\n";
    }
    start = time(NULL);
    barcode_similarity();
    if (!silent) {
        cout << "Adding edges due to barcodes similarity took: " << difftime(time(NULL), start) << "\n";
    }
    if (!silent){
        cout << "Memory after adding edges:\n\t" << get_memory_use() << "MB\n";
    }

    if (!silent) {
        cout << "Extracting clusters\n";
    }
    start = time(NULL);
    extract_clusters();
    if (!silent) {
        cout << "Extracting clusters took: " << difftime(time(NULL), start) << "\n";
    }
    if (!silent){
        cout << "Memory extracting clusters:\n\t" << get_memory_use() << "MB\n";
    }
    delete graph_ptr;
    if (!silent){
        cout << "Memory after releasing graph:\n\t" << get_memory_use() << "MB\n";
    }

    if (!silent) {
        cout << "Outputting clusters\n";
    }
    start = time(NULL);
    output_clusters();
    if (!silent) {
        cout << "Outputting clusters took: " << difftime(time(NULL), start) << "\n";
    }

}

// struct thread_data {
//    vector<vector<bool> > all_masks;
//    size_t mask_remainder;
// };
//
// void *lsh_mask_pthread(void* args) {
//     struct thread_data *my_thread_data;
//     my_thread_data = (struct thread_data *) args;
//     lsh_mask(my_thread_data->all_masks, my_thread_data->mask_remainder);
// }

void process_lsh(masked_barcode_to_barcode_id_unordered_map* lsh_ptr,
                    node_id_to_node_id_vector_of_vectors* local_graph_ptr) {
    for (auto kv = (*lsh_ptr).begin(); kv!= (*lsh_ptr).end(); kv++) {
        for (barcode_id_t bid: kv->second){
            for (node_id_t node : (*barcode_to_nodes_vector_ptr)[bid]) {
                for (barcode_id_t bid_o: kv->second){
                    if (bid == bid_o){
                        continue;
                    }
                    vector<node_id_t> good_neighbors = get_good_neighbors(node, (*barcode_to_nodes_vector_ptr)[bid_o]);
                    vector<node_id_t> result;
                    set_union((*local_graph_ptr)[node].begin(), (*local_graph_ptr)[node].end(),
                        good_neighbors.begin(), good_neighbors.end(),
                        back_inserter(result)
                    );
                    (*local_graph_ptr)[node] = move(result);
                }
            }
        }
    }
}

void process_identical_barcode_nodes(uint8_t barcode_id_remainder) {
    if (!silent) {
        stringstream stream;
        stream << "Adding edges between nodes of identical barcodes with thread " << (int)barcode_id_remainder << "\n";
        cout << stream.str();
    }
    for (barcode_id_t i = barcode_id_remainder; i < barcode_count; i+=thread_count) {
        for (node_id_t node : (*barcode_to_nodes_vector_ptr)[i]) {
            vector<node_id_t> good_neighbors = get_good_neighbors(node, (*barcode_to_nodes_vector_ptr)[i]);
            vector<node_id_t> result;
            set_union((*graph_ptr)[node].begin(), (*graph_ptr)[node].end(),
                      good_neighbors.begin(), good_neighbors.end(),
                      back_inserter(result)
                      );
            (*graph_ptr)[node] = move(result);
        }
    }
}

void lsh_mask(size_t mask_remainder) {
    time_t start;
    time_t build_time = 0, process_time = 0;
    node_id_to_node_id_vector_of_vectors local_graph(node_count);
    char masked_barcode_buffer[150];
    masked_barcode_buffer[barcode_length*2-error_tolerance] = '\0';
    stringstream stream;
    for (int i = 0; i < all_masks.size(); i++) {
        if (i % thread_count != mask_remainder) {
            continue;
        }
        start = time(NULL);
        vector<bool> mask = all_masks[i];
        masked_barcode_to_barcode_id_unordered_map lsh;
        if (!silent) {
            string current_mask_bin;
            for (bool p: mask) {
                current_mask_bin += p ? "1" : "0";
            }
            stream = stringstream();
            stream << current_mask_bin << " is assigned to thread "<< mask_remainder << "\n";
            cout << stream.str();
        }
        string masked_barcode;
        for (barcode_id_t i = 0; i < barcode_count; i++) {
            masked_barcode = mask_barcode((*barcodes_ptr)[i], mask, masked_barcode_buffer);
            if (masked_barcode != "0"){
                lsh[masked_barcode].push_back(i);
            }
        }
        build_time += difftime(time(NULL), start);
        if (!silent) {
            stream = stringstream();
            stream << "Thread "<< mask_remainder <<" built LSH in: " << difftime(time(NULL), start) << "\n";
            cout << stream.str();
        }
        start = time(NULL);
        process_lsh(&lsh, &local_graph);
        process_time += difftime(time(NULL), start);
        if (!silent) {
            stream = stringstream();
            stream << "Thread "<< mask_remainder <<" processed LSH in: " << difftime(time(NULL), start) << "\n";
            cout << stream.str();
        }
    }

    if (!silent) {
        stream = stringstream();
        stream << "On thread " << mask_remainder << " building all LSH took: " << build_time << "\n";
        stream << "On thread " << mask_remainder << " processing all LSH took: " << process_time << "\n";
        cout << stream.str();
    }
    start = time(NULL);
    if (!silent) {
        stream = stringstream();
        stream << "On thread " << mask_remainder << " merging local graph with global graph\n";
        cout << stream.str();
    }
    merge_graphs(&local_graph);
    node_id_to_node_id_vector_of_vectors().swap(local_graph);
    if (!silent) {
        stream = stringstream();
        stream << "On thread " << mask_remainder << " merging took " << difftime(time(NULL), start) <<"\n";
        cout << stream.str();
    }

}

void merge_graphs(node_id_to_node_id_vector_of_vectors* local_graph_ptr) {
    graph_lock.lock();
    if ((*graph_ptr).size() == 0) {
        (*graph_ptr) = move(*local_graph_ptr);
        graph_lock.unlock();
        return;
    }
    for (node_id_t node = 0; node < node_count; node++) {
        vector<node_id_t> result;
        set_union((*local_graph_ptr)[node].begin(), (*local_graph_ptr)[node].end(),
                    (*graph_ptr)[node].begin(), (*graph_ptr)[node].end(),
                    back_inserter(result)
        );
        (*graph_ptr)[node] = move(result);
    }
    graph_lock.unlock();
}

void barcode_similarity(){
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

    size_t mask_count = 1;
    vector<bool> mask(barcode_length*2, false);
    std::fill(mask.begin() + error_tolerance, mask.end(), true);
    for (int i = barcode_length*2; i > barcode_length*2 - error_tolerance; i--) {
        mask_count *= i;
        mask_count /= barcode_length*2 - i + 1;
    }
    all_masks.reserve(mask_count);
    do{
        all_masks.push_back(mask);

    } while (std::next_permutation(mask.begin(), mask.end()));
    if (!silent) {
        cout << "Number of masks is " << all_masks.size() << "\n";
    }

    // pthread_t* thread_array;
    // thread_data* thread_data_array;
    // pthread_attr_t attr;
    // void *status;
    // pthread_attr_init(&attr);
    // pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    thread* thread_array;
    time_t start = time(NULL);
    if (thread_count > 1) {
        // thread_array = new pthread_t[thread_count];
        // thread_data_array = new thread_data[thread_count];
        thread_array = new thread[thread_count];
        for (size_t t_id = 0; t_id < thread_count; t_id++) {
            // thread_data_array[t_id].all_masks = all_masks;
            // thread_data_array[t_id].mask_remainder = t_id;
            // pthread_create(&thread_array[t_id], &attr, lsh_mask_pthread, (void *)&thread_data_array[t_id]);
            thread_array[t_id] = thread(lsh_mask, t_id);
            stringstream stream;
            stream << "Created thread " << t_id << "\n";
            cout << stream.str();
        }
        for (size_t t_id = 0; t_id < thread_count; t_id++) {
            // pthread_join(thread_array[t_id], &status);
            thread_array[t_id].join();
            stringstream stream;
            stream << "Joined thread " << t_id << "\n";
            cout << stream.str();
        }
    } else {
        lsh_mask(0);
    }
    cout << "Building the graph on " << thread_count << " thread(s) took " <<  difftime(time(NULL), start) << "\n";
    // barcodes are no longer needed
    delete barcodes_ptr;
    if (thread_count > 1) {
        for (size_t t_id = 0; t_id < thread_count; t_id++) {
            thread_array[t_id] = thread(process_identical_barcode_nodes, t_id);
        }
        for (size_t t_id = 0; t_id < thread_count; t_id++) {
            thread_array[t_id].join();
            stringstream stream;
            stream << "Joined thread " << t_id << "\n";
            cout << stream.str();
        }
        delete [] thread_array;
    } else {
        process_identical_barcode_nodes(0);
    }
    // barcode id to node id's is no longer needed
    delete barcode_to_nodes_vector_ptr;
    delete node_to_minimizers_ptr;
}

string mask_barcode(const string& barcode, const vector<bool>& mask, char* masked_barcode_buffer){
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

vector<node_id_t> get_good_neighbors(node_id_t node, const vector<node_id_t> &neighbors){
    vector<node_id_t> good_neighbors;
    for (node_id_t neighbor : neighbors) {
        if (node != neighbor && !unmatched_minimimizers(node, neighbor)) {
            good_neighbors.push_back(neighbor);
        }
    }
    return good_neighbors;
}

bool unmatched_minimimizers(node_id_t node_id, node_id_t neighbor_id){
    int matched_minimimizers_1 = 0;
    int matched_minimimizers_2 = 0;
    for (int i =0; i < minimizer_count; i++) {
        matched_minimimizers_1 += (*node_to_minimizers_ptr)[node_id][i] == (*node_to_minimizers_ptr)[neighbor_id][i];
        matched_minimimizers_2 += (*node_to_minimizers_ptr)[node_id][i+minimizer_count] == (*node_to_minimizers_ptr)[neighbor_id][i+minimizer_count];
    }
    return !(matched_minimimizers_1 >= minimizer_threshold && matched_minimimizers_2 >= minimizer_threshold);
}

void extract_clusters(){
    vector<bool> pushed(node_count, false);
    stack<node_id_t> opened;
    node_to_cluster_vector.reserve(node_count);
    cluster_count = 0;

    for (node_id_t node = 0; node < node_count; node++) {
        if (!pushed[node]) {
            opened.push(node);
            pushed[node] = true;
            while(!opened.empty()) {
                node_id_t current_node = opened.top();
                opened.pop();
                node_to_cluster_vector[current_node] = cluster_count;
                for (node_id_t neighbor: (*graph_ptr)[current_node]) {
                    if (!pushed[neighbor]) {
                        opened.push(neighbor);
                        pushed[neighbor] = true;
                        node_to_cluster_vector[neighbor] = cluster_count;
                    }
                }
            }
            cluster_count++;
        }
    }
}

void output_clusters(){
    read_id_t current_read = 0;
    ofstream clusters;
    ifstream fastq1;
    ifstream fastq2;
    fastq1.open (input_1);
    fastq2.open (input_2);
    string name_1, quality_1, sequence_1, name_2, quality_2, sequence_2, trash;

    clusters = ofstream(output_prefix + "cluster");
    while (getline(fastq1, name_1)) {
        getline(fastq1, sequence_1);
        getline(fastq1, trash);
        getline(fastq1, trash);
        getline(fastq2, name_2);
        getline(fastq2, sequence_2);
        getline(fastq2, trash);
        getline(fastq2, trash);

        node_id_t current_read_node = read_to_node_vector[current_read];
        clusters << node_to_cluster_vector[current_read_node] << "\t" << current_read_node << "\t" << current_read << "\t";
        clusters << name_1 << "\t" << sequence_1 << "\t" << quality_1 << "\t";
        clusters << name_2 << "\t" << sequence_2 << "\t" << quality_2 << "\n";
        current_read++;
    }
}
