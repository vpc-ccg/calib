//
// Created by borabi on 19/12/17.
//

#include <string>
#include <stdlib.h>
#include <iostream>
#include "global.h"



#ifndef CALIB_COMMANDLINE_H
#define CALIB_COMMANDLINE_H

void parse_flags(int argc, char *argv[]);
void print_flags(std::ofstream &out);
void print_help();

#endif //CALIB_COMMANDLINE_H
