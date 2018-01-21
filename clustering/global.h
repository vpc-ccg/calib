//
// Created by borabi on 19/12/17.
//

#include <iostream>
#include <fstream>
#include <stdio.h>

#ifndef BARGOAT_GLOBAL_H
#define BARGOAT_GLOBAL_H


extern std::string input_1;
extern std::string input_2;
extern std::string output_prefix;
extern int barcode_length;
extern int minimizer_count;
extern int kmer_size;
extern int error_tolerance;
extern int minimizer_threshold;
extern int thread_count;
extern bool silent;
extern bool keep_qual;
extern std::ofstream dog;




#endif //BARGOAT_GLOBAL_H
