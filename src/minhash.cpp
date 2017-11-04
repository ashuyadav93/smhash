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

std::vector <std::string> get_all_kmers(std::string input, unsigned int k_size);

int main(int argc, char **argv) {
    std::vector <std::string> kmers = get_all_kmers("hellosaraj", 3);
    for (std::vector<std::string>::iterator it = kmers.begin(); it != kmers.end(); ++it) {
        log(*it + "\n");
    }
    char *inp_buf = "hello";
    void *out_buf = malloc(strlen(inp_buf) * sizeof(char));
    MurmurHash3_x64_128(inp_buf, strlen(inp_buf), 123213, out_buf);
    long hash_val = (long) out_buf;
    std::cout << "Hash is " << hash_val;
    CountEstimator *ce = new CountEstimator(10, 3, true, 0, false);
    ce->add("llo");
    ce->print_sketch();
    return 0;
}


std::vector <std::string> get_all_kmers(std::string input, unsigned int k_size) {
    std::vector <std::string> kmers;
    if (input.length() == 0) {
        log("Can't find kmers for empty input");
        return kmers;
    }
    int len = input.length();
    int num_kmers = len - k_size + 1;
    for (int i = 0; i < num_kmers; i++) {
        kmers.push_back(input.substr(i, k_size));
    }
    return kmers;
}

void log(std::string s) {
    //TODO replace this with a logging lib later, or write custom outfile interface
    std::cout << s;
}