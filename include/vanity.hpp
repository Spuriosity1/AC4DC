#ifndef VANITY_CXX_H
#define VANITY_CXX_H

#include <fstream>
#include <iostream>

void print_file(const char* fname){
    std::ifstream ifs(fname, std::ifstream::in);

    char c = ifs.get();
    while (ifs.good()) {
        std::cout << c;
        c = ifs.get();
    }
    ifs.close();
}

#endif