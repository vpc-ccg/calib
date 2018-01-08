//
// Created by borabi on 19/12/17.
//

#include "parse_flags.h"

using namespace std;

// Parameter definitions
string input_prefix = "";
string output_prefix = "";
int barcode_length = -1;
int minimizer_count = -1;
int error_tolerance = -1;
int minimizer_threshold = -1;
int kmer_size = -1;

void parse_flags(int argc, char *argv[]){

    for (int i = 0; i < argc; i++){
        string current_param(argv[i]);
        if (current_param == "-i" || current_param == "--input-prefix") {
            input_prefix = string(argv[i+1]);
        }

        if (current_param == "-o" || current_param == "--output-prefix") {
            output_prefix = string(argv[i+1]);
        }

        if (current_param == "-l" || current_param == "--barcode-length") {
            barcode_length = atoi(argv[i+1]);
        }

        if (current_param == "-m" || current_param == "--minimizer-count") {
            minimizer_count = atoi(argv[i+1]);
        }

        if (current_param == "-k" || current_param == "--kmer-size") {
            kmer_size = atoi(argv[i+1]);
        }

        if (current_param == "-e" || current_param == "--error-tolerance") {
            error_tolerance = atoi(argv[i+1]);
        }

        if (current_param == "-t" || current_param == "--minimizer-threshold") {
            minimizer_threshold = atoi(argv[i+1]);
        }
    }

    if (barcode_length < 0 || minimizer_count < 0 || error_tolerance < 0 || minimizer_threshold < 0 || kmer_size < 0){
        cout << "Missing parameters!\n";
        exit(-1);
    }
    if (input_prefix == "" || output_prefix == ""){
        cout << "Missing parameters!\n";
        exit(-1);
    }
}
