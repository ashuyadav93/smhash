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

int i_size(double j, int k, int n) {
    return int(n*j/double(1-j) + k);
}

void log(std::string s);

std::vector <std::string> get_all_kmers(std::string input, unsigned int k_size);

int main(int argc, char **argv) {
    std::string alpha = "ACTG";

    std::vector< int >ir;
    std::vector<double> true_jaccard,estimate_jaccards,containment_jaccards;
    std::vector<int> sequence_len1,sequence_len2;
    int index = 0;

    for(double i=0; i<=1; i += 0.005) {
        ir.push_back(i_size(i, 11, 10000));
    }

    for(int k=0; k<ir.size();k++) {
        std::string seq1 = "";
        std::string seq2 = "";
        std::string common_string = "";
        //cout<<ir[k]<<endl;
        srand(time(0));
        for(unsigned int i = 0; i < 10000; ++i)
        {
        seq1 += alpha[rand() % 3];
        }

        for(unsigned int i = 0; i < 15; ++i)
        {
        seq2 += alpha[rand() % 3];
        }

        for(unsigned int i = 0; i < ir[k]; i++)
        {
        common_string += alpha[rand() % 3];
        }

        seq1 += common_string;
        seq2 += common_string;

        std::hash <std::string> hash_fn;
        CountEstimator *ce1 = new CountEstimator(100, 11, true, 0, false, hash_fn);
        ce1->add_sequence(seq1);
        // ce1->print_sketch();
       
        CountEstimator *ce2 = new CountEstimator(100, 11, true, 0, false, hash_fn);
        ce2->add_sequence(seq2);
        // ce2->print_sketch();
        
        bloom_parameters parameters;

        //How many elements roughly do we expect to insert?
        parameters.projected_element_count = 10000 + ir[k];

        //Maximum tolerable false positive probability? (0,1)
        parameters.false_positive_probability = 0.001; // 1 in 1000

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
        std::vector<std::string>str_lista = get_all_kmers(seq1, 11);
        std::vector<std::string>str_listb = ce2->get_kmers();
        vector<string> kmer_lista = get_all_kmers(seq1, 11);
        vector<string> kmer_listb = get_all_kmers(seq2, 11);
        set<string>kmer_a (kmer_lista.begin(), kmer_lista.end());
        set<string>kmer_b (kmer_listb.begin(), kmer_listb.end());
        vector<string>intersection,union_set;

        set_intersection(kmer_a.begin(),kmer_a.end(),kmer_b.begin(),kmer_b.end(), std::back_inserter(intersection));
        set_union(kmer_a.begin(),kmer_a.end(),kmer_b.begin(),kmer_b.end(), std::back_inserter(union_set));

        true_jaccard.push_back(intersection.size()/double(union_set.size()));
        int intersect_est = 0; // intersection estimate
        int kmer_len = 0;
        
        // Insert into Bloom Filter
        set<string>::iterator iter = kmer_a.begin();
        for (iter = kmer_a.begin(); iter != kmer_a.end(); iter++)
        {   
            if(!filter.contains(*iter)) {
                kmer_len += 1;
                filter.insert(*iter);
            }
        }

        for (unsigned int i = 0; i < str_listb.size(); ++i)
        {
            if (filter.contains(str_listb[i]))
            {
                intersect_est+=1;
            }
        }
        int h = ce1->get_hashsize();
        intersect_est -= round(parameters.false_positive_probability*h);  // adjust for the false positive rate
        
        double containment_est = intersect_est / double(h);
        double jaccard_est = (kmer_b.size() * containment_est) / (kmer_len + kmer_b.size() - (kmer_b.size() * containment_est));  // estimate of the containment index
        
        estimate_jaccards.push_back(ce1->calc_jaccard_distance(ce2));
        containment_jaccards.push_back(jaccard_est);
        sequence_len1.push_back(seq1.size());
        sequence_len2.push_back(seq2.size());

        index+=1;

    }

        for(int i=0;i<index;i++) {
            cout<<sequence_len1[i] << "   ";
            cout<<sequence_len2[i] << "   ";
            cout<<true_jaccard[i] << "   ";
            cout<<estimate_jaccards[i] << "   ";
            cout<<containment_jaccards[i] << "   ";
            cout << "\n";
        }
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
