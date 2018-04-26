//
// Created by borabi on 19/12/17.
//

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>

#ifndef CALIB_GLOBAL_H
#define CALIB_GLOBAL_H

typedef size_t node_id_t;


extern std::string input_1;
extern std::string input_2;
extern std::string output_prefix;
extern int barcode_length;
extern int ignored_sequence_prefix_length;
extern int minimizer_count;
extern int kmer_size;
extern int error_tolerance;
extern int minimizer_threshold;
extern int thread_count;
extern bool silent;
extern bool keep_qual;
extern bool bc_format;
extern bool debug;
extern bool no_triplets;
extern std::ofstream dog;
extern std::ofstream node_dog;
extern std::vector<node_id_t> debug_nodes;



#endif //CALIB_GLOBAL_H
