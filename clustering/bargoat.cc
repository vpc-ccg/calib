//
// Created by borabi on 19/12/17.
//

#include "string"
#include <iostream>
#include <sstream>

#include "global.h"
#include "commandline.h"
#include "extract.h"
#include "cluster.h"

using namespace std;

ofstream dog;

int main(int argc, char *argv[]){
    parse_flags(argc, argv);

    stringstream ss;
    if (silent) {
        cout.rdbuf(ss.rdbuf());
    }


    dog = ofstream(output_prefix + ".cluster.log");
    print_flags(dog);
    // print_flags(&cout);

    cout << "Extracting minimizers and barcodes...\n";
    extract_barcodes_and_minimizers();

    cout << "Clustering...\n";
    cluster();

    cout << "All done! Have good day!\n";
}
