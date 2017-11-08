//
// Created by Vaibhav Rustagi on 11/04/17.
//
#include <iostream>
#include <string>
#include <vector>
#include "../../murmurhash3/murmurhash3.h"
#include <cstring>
#include <limits>
#include "../../common/CountEstimator.cpp"
#include "bloom_filter.hpp"
#include <cmath>
#include <set>
#include <ctime>
using namespace std;

void log(std::string s);

std::vector <std::string> get_all_kmers(std::string input, unsigned int k_size);

int main(int argc, char **argv) {
    std::string alpha = "ACTG";

    std::string seq1 = "";
    std::string seq2 = "";
    std::string common_string = "";
    
    srand(time(0));
    for(unsigned int i = 0; i < 1000; ++i)
    {
    seq1 += alpha[rand() % 3];
    }

    for(unsigned int i = 0; i < 15; ++i)
    {
    seq2 += alpha[rand() % 3];
    }

    for(unsigned int i = 0; i < 15; ++i)
    {
    common_string += alpha[rand() % 3];
    }

    // common_string += 
    // common_string = ''.join(np.random.choice(['A', 'C', 'T', 'G'], i_size))
    // seq1 = ''.join(np.random.choice(['A', 'C', 'T', 'G'], 1000)) + common_string
    // # Make seq2 a smaller sequence than seq1
    // seq2 = ''.join(np.random.choice(['A', 'C', 'T', 'G'], 15)) + common_string

    std::hash <std::string> hash_fn;
    CountEstimator *ce1 = new CountEstimator(50, 5, true, 0, false, hash_fn);
    ce1->add_sequence(seq1);
    ce1->print_sketch();
   
    CountEstimator *ce2 = new CountEstimator(50, 5, true, 0, false, hash_fn);
    ce2->add_sequence(seq2);
    ce2->print_sketch();

    bloom_parameters parameters;

    //How many elements roughly do we expect to insert?
    parameters.projected_element_count = 1000;

    //Maximum tolerable false positive probability? (0,1)
    parameters.false_positive_probability = 0.001; // 1 in 10000

    //Simple randomizer (optional)
    // parameters.random_seed = 0xA5A5A5A5;

    if (!parameters)
    {
        std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
        return 1;
    }

    parameters.compute_optimal_parameters();

    //Instantiate Bloom Filter
    bloom_filter filter(parameters);
    std::vector<std::string>str_lista = get_all_kmers(seq1, 5);
    std::vector<std::string>str_listb = ce2->get_kmers();
    vector<string> kmer_lista = get_all_kmers(seq1, 5);
    vector<string> kmer_listb = get_all_kmers(seq2, 5);
    set<string>kmer_a (kmer_lista.begin(), kmer_lista.end());
    set<string>kmer_b (kmer_listb.begin(), kmer_listb.end());
    vector<string>intersection,union_set;

    set_intersection(kmer_a.begin(),kmer_a.end(),kmer_b.begin(),kmer_b.end(), std::back_inserter(intersection));
    set_union(kmer_a.begin(),kmer_a.end(),kmer_b.begin(),kmer_b.end(), std::back_inserter(union_set));
    double true_jaccard = intersection.size()/double(union_set.size());

    int intersect_est = 0; // intersection estimate
    int kmer_len = 0;
   
    // Insert into Bloom Filter
    
    for (unsigned int i = 0; i < str_lista.size(); ++i)
    {   
        if(!filter.contains(str_lista[i])) {
            kmer_len += 1;
            filter.insert(str_lista[i]);
        }
    }

    for (unsigned int i = 0; i < str_listb.size(); ++i)
    {
        std::cout << str_listb[i] << " ";
        if (filter.contains(str_listb[i]))
        {
            std::cout << "BF contains: " << str_listb[i] << std::endl;
            intersect_est+=1;
        }
    }

    int h = ce1->get_hashsize();
    // cout<< intersect_est << " " << h;
    intersect_est -= round(parameters.false_positive_probability*h);  // adjust for the false positive rate
    cout << intersect_est << " ";
    double containment_est = intersect_est / double(h);
    double jaccard_est = str_listb.size() * containment_est / (kmer_len + str_listb.size() - str_listb.size() * containment_est);  // estimate of the containment index
    std::cout<< containment_est << std::endl;
    std::cout<< jaccard_est << std::endl;
    std::cout<< ce1->calc_jaccard_distance(ce2) << std::endl;
    std::cout<< true_jaccard << endl;

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
