//
// Created by borabi on 19/12/17.
//

#include <iostream>
#include <fstream>
#include <stdio.h>


#ifndef BARGOAT_GLOBAL_H
#define BARGOAT_GLOBAL_H




extern std::string input_prefix;
extern std::string output_prefix;
extern int barcode_length;
extern int minimizer_count;
extern int kmer_size;
extern int error_tolerance;
extern int minimizer_threshold;
extern std::ofstream log;




struct Node{
    std::string barcode;
    uint64_t *minimizers_1;
    uint64_t *minimizers_2;
    int id;
};

struct NodeHash {
    size_t operator()(const Node *node) const{
        size_t result = std::hash<std::string>()(node->barcode);
        for (int i = 0; i < minimizer_count; i++){
            result ^= std::hash<int>()(node->minimizers_1[i]) << i;
        }
        for (int i = 0; i < minimizer_count; i++){
            result ^= std::hash<int>()(node->minimizers_2[i]) << i;
        }
        return result;
    }
};

struct NodeEqual {
    bool operator()(const Node *lhs, const Node *rhs) const{
        bool result = lhs->barcode == rhs->barcode;
        for (int i = 0; i < minimizer_count; i++){
            result = result && (lhs->minimizers_1[i] == rhs->minimizers_1[i]);
        }
        for (int i = 0; i < minimizer_count; i++){
            result = result && (lhs->minimizers_2[i] == rhs->minimizers_2[i]);
        }
        return result;
    }
};


#endif //BARGOAT_GLOBAL_H

