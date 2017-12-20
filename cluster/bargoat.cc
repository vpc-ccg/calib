//
// Created by borabi on 19/12/17.
//

#include "string"
#include <iostream>

#include "global.h"
#include "parse_flags.h"
#include "extract_barcodes_and_minimizers.h"
#include "cluster.h"

using namespace std;

// Parameter definitions
string input_prefix;
string output_prefix;
int barcode_length;
int minimizer_count;
int error_tolerance;
int minimizer_threshold;
int kmer_size;
ofstream log;


int main(int argc, char *argv[]){
    parse_flags(argc, argv);
    log = ofstream(output_prefix + ".log");
    log << "Parameters:\n";
    log << "\tinput_prefix:\t" << input_prefix << "\n";
    log << "\toutput_prefix:\t" << output_prefix << "\n";
    log << "\tbarcode_length:\t" << barcode_length << "\n";
    log << "\tminimizer_count:\t" << minimizer_count << "\n";
    log << "\tkmer_size:\t" << kmer_size << "\n";
    log << "\terror_tolerance:\t" << error_tolerance << "\n";
    log << "\tminimizer_threshold:\t" << minimizer_threshold << "\n";

    extract_barcodes_and_minimizers();

    cluster();

}