//
// Created by borabi on 19/12/17.
//

#include "global.h"
#include "cluster.h"

#ifndef BARGOAT_EXTRACT_H
#define BARGOAT_EXTRACT_H


minimizer_t minimizer(std::string& seq, int start, int length);
void extract_barcodes_and_minimizers();

#endif //BARGOAT_EXTRACT_H
