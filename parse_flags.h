//
// Created by borabi on 19/12/17.
//
#include <string>
#include <stdlib.h>
#include <iostream>
using namespace std;

#ifndef BARGOAT_PARSE_FLAGS_H
#define BARGOAT_PARSE_FLAGS_H

extern string input_prefix;
extern string output_prefix;
extern int barcode_length;
extern int minimizer_count;
extern int error_tolerance;
extern int minimizer_threshold;


int parse_flags(int argc, char *argv[]);



#endif //BARGOAT_PARSE_FLAGS_H
