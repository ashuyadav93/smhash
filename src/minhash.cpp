//
// Created by Saraj Munjal on 10/25/17.
//
#include <iostream>
#include <string>
#include <vector>
#include "murmurhash3/murmurhash3.h"
#include <cstring>
#include <limits>
#include "CountEstimator.cpp"

void log(std::string s);


int main(int argc, char **argv) {
    CountEstimator *ce1 = new CountEstimator(5, 3, true, 0, false);
    ce1->add_sequence("TCAGTTATATTAGCAT");
    ce1->print_sketch();

    CountEstimator *ce2 = new CountEstimator(5, 3, true, 0, false);
    ce2->add_sequence("TCAGTTATATTAGCTA");
    ce2->print_sketch();

    std::cout << "Jaccard distance: " << ce1->calc_jaccard_distance(ce2);
    return 0;
}


void log(std::string s) {
    //TODO replace this with a logging lib later, or write custom outfile interface
    std::cout << s;
}