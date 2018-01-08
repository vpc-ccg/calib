//
// Created by borabi on 19/12/17.
//

#include "string"
#include <iostream>

#include "global.h"
#include "parse_flags.h"
#include "extract.h"
#include "cluster.h"

using namespace std;

ofstream dog;

int main(int argc, char *argv[]){
    parse_flags(argc, argv);
    dog = ofstream(output_prefix + ".log");
    dog << "Parameters:\n";
    dog << "\tinput_prefix:\t" << input_prefix << "\n";
    dog << "\toutput_prefix:\t" << output_prefix << "\n";
    dog << "\tbarcode_length:\t" << barcode_length << "\n";
    dog << "\tminimizer_count:\t" << minimizer_count << "\n";
    dog << "\tkmer_size:\t" << kmer_size << "\n";
    dog << "\terror_tolerance:\t" << error_tolerance << "\n";
    dog << "\tminimizer_threshold:\t" << minimizer_threshold << "\n";

    extract_barcodes_and_minimizers();
    cout << "Done extracting\n";

    cluster();

}
