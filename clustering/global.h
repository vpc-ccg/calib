//
// Created by borabi on 19/12/17.
//

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <string.h>


#ifndef CALIB_GLOBAL_H
#define CALIB_GLOBAL_H

typedef uint32_t node_id_t;


extern std::string input_1;
extern std::string input_2;
extern std::string output_prefix;
extern int barcode_length_1;
extern int barcode_length_2;
extern int ignored_sequence_prefix_length;
extern int minimizer_count;
extern int kmer_size;
extern int error_tolerance;
extern int minimizer_threshold;
extern int thread_count;
extern bool silent;
extern bool no_sort;

int get_memory_use();

#endif //CALIB_GLOBAL_H
