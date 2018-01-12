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
    dog = ofstream(output_prefix + ".clusters");
    dog << "Parameters:\n";
    dog << "\tinput_1:\t" << input_1 << "\n";
    dog << "\tinput_2:\t" << input_2 << "\n";
    dog << "\toutput_prefix:\t" << output_prefix << "\n";
    dog << "\tbarcode_length:\t" << barcode_length << "\n";
    dog << "\tminimizer_count:\t" << minimizer_count << "\n";
    dog << "\tkmer_size:\t" << kmer_size << "\n";
    dog << "\terror_tolerance:\t" << error_tolerance << "\n";
    dog << "\tminimizer_threshold:\t" << minimizer_threshold << "\n";
    dog << "\tthreads:\t" << thread_count << "\n";

    extract_barcodes_and_minimizers();
    cout << "Done extracting\n";

    cluster();

}
