//
// Created by borabi on 20/12/17.
//

#include "cluster.h"

using namespace std;


void cluster(){
    for (auto it : node_to_read){
        print_node(it.first, it.second);
    }

    

}


void print_node(Node node, vector<int> reads){
    log << node.id << "\t" << node.barcode << "\t" ;
    for (int i =0; i < minimizer_count; i++)
        log << node.minimizers_1[i] << "\t";
    for (int i =0; i < minimizer_count; i++)
        log << node.minimizers_2[i] << "\t";
    for (auto read: reads)
        log << read << "\t";
    log << "\n";

}
