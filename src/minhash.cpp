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
#include <functional>

void log(std::string s);


int main(int argc, char **argv) {
    std::hash <std::string> hash_fn;
    CountEstimator *ce1 = new CountEstimator(5, 3, true, 0, false, hash_fn);
    ce1->add_sequence("ATTATTATTATT");
    ce1->print_sketch();

//    std::vector<unsigned long> v;
//    v.push_back(4954146047334);
//    v.push_back(18446744073709551615);
//    v.push_back(18446744073709551615);
//    v.push_back(18446744073709551615);
//    v.push_back(18446744073709551615);
////    int idx = ce1->bin_search(v, 0, 4, 3);
////    v.insert(v.begin() + idx, 3);
//    std::cout << "Bin search res: " << ce1->bin_search(v, 0, 4, 4898558728737) << "\n";

    CountEstimator *ce2 = new CountEstimator(5, 3, true, 0, false, hash_fn);
    ce2->add_sequence("TTATTATTATTA");
    ce2->print_sketch();

    std::cout << "Jaccard distance: " << ce1->calc_jaccard_distance(ce2);
    return 0;
}


void log(std::string s) {
    //TODO replace this with a logging lib later, or write custom outfile interface
    std::cout << s;
}