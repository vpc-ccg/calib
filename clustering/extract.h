//
// Created by borabi on 19/12/17.
//

#include "global.h"
#include "cluster.h"

#ifndef CALIB_EXTRACT_H
#define CALIB_EXTRACT_H

std::string minimizer_t_to_dna(minimizer_t minimizer, size_t size);
void make_invalid_minimizer_vector();
minimizer_t minimizer(std::string& seq, int start, int length);
void extract_barcodes_and_minimizers();

#endif //CALIB_EXTRACT_H
