//
// Created by borabi on 20/12/17.
//
#include "global.h"
#include <unordered_map>
#include <vector>
// #include <bitset>
#include <unordered_set>
#include <algorithm>


#ifndef CLUSTER_H
#define CLUSTER_H


typedef size_t node_id_t;
typedef size_t read_id_t;
typedef uint64_t minimizer_t;

extern size_t node_count;
extern size_t read_count;

struct Read{
    std::string name_1;
    std::string sequence_1;
    std::string quality_1;

    std::string name_2;
    std::string sequence_2;
    std::string quality_2;
};

struct Node{
    std::string barcode;
    std::vector<minimizer_t> minimizers_1; //(minimizer_count, (minimizer_t) -1);
    std::vector<minimizer_t> minimizers_2; //(minimizer_count);
    Node(){
        barcode = "";
        minimizers_1 = std::vector<minimizer_t> (minimizer_count);
        minimizers_2 = std::vector<minimizer_t> (minimizer_count);
    }
};

struct NodeHash {
    size_t operator()(const Node& node) const{
        size_t result = std::hash<std::string>()(node.barcode);
        for (int i = 0; i < minimizer_count; i++){
            result ^= std::hash<int>()(node.minimizers_1[i]) << i;
        }
        for (int i = 0; i < minimizer_count; i++){
            result ^= std::hash<int>()(node.minimizers_2[i]) << i;
        }
        return result;
    }
};

struct NodeEqual {
    bool operator()(const Node& lhs, const Node& rhs) const{
        bool result = lhs.barcode == rhs.barcode;
        for (int i = 0; i < minimizer_count; i++){
            result = result && (lhs.minimizers_1[i] == rhs.minimizers_1[i]);
        }
        for (int i = 0; i < minimizer_count; i++){
            result = result && (lhs.minimizers_2[i] == rhs.minimizers_2[i]);
        }
        return result;
    }
};

// reads will store the actual paired-end fastq files contents
typedef std::vector<Read> read_vector;
// node_vector and node_id_to_read_id will replace node_to_read_id once reading fastq files is done
typedef std::vector<Node> node_vector;
typedef std::vector<std::vector<read_id_t>> node_id_to_read_id_vector;

// node_id_to_node_id maps a node to its neighbors
typedef std::vector<std::unordered_set<node_id_t>> node_id_to_node_id_vector;
// masked_barcode_to_node_id is an LSH dictionary
typedef std::unordered_map<std::string, std::vector<node_id_t>> masked_barcode_to_node_id_unordered_map;

extern read_vector reads;
extern node_vector nodes;
extern node_id_to_read_id_vector node_to_read_vector;

void cluster();
void remove_edges_of_unmatched_minimizers(node_id_to_node_id_vector &adjacency_sets);
void barcode_similarity(node_id_to_node_id_vector &adjacency_sets);
std::string mask_barcode(const std::string& barcode, const std::vector<bool>& mask);
void print_node(node_id_t node_id);

#endif //CLUSTER_H
