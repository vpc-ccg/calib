//
// Created by borabi on 20/12/17.
//
#include "global.h"
#include <unordered_map>
#include <vector>


extern int node_count;
extern int read_count;
extern std::unordered_map<Node* , std::vector<int>, NodeHash, NodeEqual> node_to_read;

void cluster();

#ifndef BARGOAT_CLUSTER_H
#define BARGOAT_CLUSTER_H

#endif //BARGOAT_CLUSTER_H
