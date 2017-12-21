//
// Created by borabi on 20/12/17.
//
#include "global.h"
#include <unordered_map>
#include <vector>



#ifndef CLUSTER_H
#define CLUSTER_H

extern int node_count;
extern int read_count;

struct Node{
    Node() {
        barcode = "";
        minimizers_1 = new uint64_t [minimizer_count];
        minimizers_2 = new uint64_t [minimizer_count];
        id = node_count++;
    }
    std::string barcode;
    uint64_t *minimizers_1;
    uint64_t *minimizers_2;
    int id;
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

typedef std::unordered_map<Node, std::vector<int>, NodeHash, NodeEqual> node_map;
extern node_map node_to_read;


void cluster();
void print_node(Node node, std::vector<int> reads);


#endif //CLUSTER_H
