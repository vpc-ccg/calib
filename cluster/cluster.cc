//
// Created by borabi on 20/12/17.
//

#include "cluster.h"

using namespace std;


void cluster(){
    cout << minimizer_count << "\n";
//    log << "node->id" << "\t" << "node->barcode" << "\t" << string("mini_1\t", minimizer_count) << string("mini_2\t", minimizer_count) ;
//    log << "\n";
    cout << node_count << "\t" << read_count << "\n";
    for (auto it : node_to_read){
        Node* node = it.first;

        log << node->id << "\t" << node->barcode << "\t" ;
        for (int i =0; i < minimizer_count; i++)
            log << node->minimizers_1[i] << "\t";
        for (int i =0; i < minimizer_count; i++)
            log << node->minimizers_2[i] << "\t";
        log << "\n";

        for (auto read : it.second)
            log << read << "\n";
    }

}
