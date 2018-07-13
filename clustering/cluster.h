//
// Created by borabi on 20/12/17.
//
#include "global.h"
#include <unordered_map>
#include <vector>
#include <unordered_set>

#ifndef CLUSTER_H
#define CLUSTER_H

typedef uint32_t read_id_t;
typedef uint32_t barcode_id_t;

typedef std::string barcode_t;
typedef uint32_t minimizer_t;

extern node_id_t node_count;
extern read_id_t read_count;
extern barcode_id_t barcode_count;

struct Read {
    std::string name_1;
    std::string sequence_1;
    std::string quality_1;

    std::string name_2;
    std::string sequence_2;
    std::string quality_2;
};

struct Node {
    std::string barcode;
    std::vector<minimizer_t> minimizers;     //(minimizer_count);
    Node(){
        barcode = "";
        minimizers = std::vector<minimizer_t> (minimizer_count*2);
    }
};

struct NodeHash {
    size_t operator()(const Node& node) const {
        size_t result = std::hash<std::string>()(node.barcode);
        for (int i = 0; i < minimizer_count*2; i++) {
            result ^= std::hash<int>()(node.minimizers[i]) << i;
        }
        return result;
    }
};

struct NodeEqual {
    bool operator()(const Node& lhs, const Node& rhs) const {
        bool result = lhs.barcode == rhs.barcode;
        for (int i = 0; i < minimizer_count*2; i++) {
            result = result && (lhs.minimizers[i] == rhs.minimizers[i]);
        }
        return result;
    }
};

// reads will store the actual paired-end fastq files contents
typedef std::vector<Read> read_vector;
extern read_vector reads;
// to find read id's of a node ID, use this
typedef std::vector<std::vector<read_id_t> > node_id_to_read_ids_vector;
extern node_id_to_read_ids_vector node_to_reads_vector;
// to find nodes of a unique barcode
typedef std::vector<barcode_t> barcode_vector;
extern barcode_vector barcodes;
typedef std::vector<std::vector<node_id_t> > barcode_id_to_node_ids_vector;
extern barcode_id_to_node_ids_vector barcode_to_nodes_vector;
// to find minimizers of a node ID, use this
typedef std::vector<std::vector<minimizer_t> > node_id_to_minimizers_vector;
extern node_id_to_minimizers_vector node_to_minimizers;
// Use this for the adjacency lists representation of the graph
typedef std::vector<std::vector<node_id_t> > node_id_to_node_id_vector_of_vectors;
// Use this as an LSH dictionary to find similar barcodes
typedef std::unordered_map<barcode_t, std::vector<barcode_id_t> > masked_barcode_to_barcode_id_unordered_map;

void cluster();
void barcode_similarity(node_id_to_node_id_vector_of_vectors &adjacency_lists);
std::string mask_barcode(const std::string& barcode, const std::vector<bool>& mask);
void remove_edges_of_unmatched_minimizers(node_id_to_node_id_vector_of_vectors &adjacency_lists);
bool unmatched_minimimizers(node_id_t node_id, node_id_t neighbor_id);
void extract_clusters(node_id_to_node_id_vector_of_vectors &adjacency_lists);
void print_node(node_id_t node_id);
void process_lsh(masked_barcode_to_barcode_id_unordered_map &lsh,
                 node_id_to_node_id_vector_of_vectors &adjacency_lists,
                 size_t reminder);

#endif //CLUSTER_H
