//
// Created by borabi on 19/12/17.
//

#include "string"
#include <iostream>

#include "global.h"
#include "commandline.h"
#include "extract.h"
#include "cluster.h"

using namespace std;

int main(int argc, char *argv[]){
    parse_flags(argc, argv);

    if (!silent) {
        cout << "Extracting minimizers and barcodes...\n";
    }
    extract_barcodes_and_minimizers();
    if (!silent && print_mem){
        cout << "Memory after exiting extract_barcodes_and_minimizers():\n\t" << get_memory_use() << "MB\n";
    }

    if (!silent) {
        cout << "Clustering...\n";
    }
    cluster();
    if (!silent) {
        cout << "All done! Have a good day!\n";
    }
}
