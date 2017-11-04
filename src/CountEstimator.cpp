//
// Created by Saraj Munjal on 11/3/17.
//
#include <utility>
#include <algorithm>
#include <cmath>
#include "murmurhash3/murmurhash3.h"

class CountEstimator {
private:
    const long max_prime = 9999999999971L;
    const uint32_t SEED = 123213;
    int ksize;
    int n;
    //FILE* input_file_name=None,
    bool save_kmers = false;
    std::vector<long> *predef_hash_list = NULL;
    bool rev_comp = false;
    std::vector<long> sketch_counts;
    std::vector<long> sketch_hashes;
    std::vector <std::string> sketch_kmers;

    long hash(std::string kmer) {
        const char *inp_buf = kmer.c_str();
        void *out_buf = malloc(strlen(inp_buf) * sizeof(char));
        MurmurHash3_x64_128(inp_buf, strlen(inp_buf), SEED, out_buf);
        return (long) out_buf;
    }

    bool is_prime(long num) {
        //Check if a number is prime.
        if (num < 2) {
            return false;
        }
        if (num == 2) {
            return true;
        }
        if (num % 2 == 0) {
            return false;
        }
        long sqrt = (long) std::floor(std::sqrt(num));
        for (int i = 3; i < sqrt + 1; i += 2) {
            if (num % i == 0) {
                return false;
            }
        }
        return true;
    }


    long get_prime_lt_x(long target) {
        /* Backward-find a prime smaller than (or equal to) target.
        Step backwards until a prime number (other than 2) has been
        found.

                Arguments: target -- the number to step backwards from */

        if (target <= 1) {
            return 1;
        }
        long i = target;
        if (i % 2) {
            i -= 1;
        }
        while (i > 0) {
            if (is_prime(i)) {
                return i;
            }
            i -= 2;
        }
        if (i <= 0) {
            return -1;
        }
    }


    std::string reverse_complement(std::string str) {
        std::string rev = str;
        std::reverse(rev.begin(), rev.end());
        int len = rev.length();
        for (int i = 0; i < len; i++) {
            switch (rev[i]) {
                case 'G':
                    rev[i] = 'C';
                    break;

                case 'g':
                    rev[i] = 'c';
                    break;

                case 'C':
                    rev[i] = 'G';
                    break;

                case 'c':
                    rev[i] = 'g';
                    break;

                case 'T':
                    rev[i] = 'A';
                    break;

                case 't':
                    rev[i] = 'a';
                    break;

                case 'A':
                    rev[i] = 'T';
                    break;

                case 'a':
                    rev[i] = 't';
                    break;
            }
        }
        return rev;
    }

    void resize_vects_if_reqd() {
        if (sketch_hashes.size() > n) {
            sketch_hashes.resize(n);
            sketch_counts.resize(n);
            if (save_kmers) {
                sketch_kmers.resize(n);
            }
        }
    }


    void insert_vect_single_long(std::vector<long> *v, int pos, long val) {
        std::vector<long>::iterator it;
        it = v->begin() + pos;
        v->insert(it, val);
    }

    void insert_vect_single_str(std::vector <std::string> *v, int pos, std::string val) {
        std::vector<std::string>::iterator it;
        it = v->begin() + pos;
        v->insert(it, val);
    }

    void print_vec(std::vector<long> *v) {
        int len = v->size();
        for (int i = 0; i < len; i++) {
            std::cout << (*v)[i] << " ";
        }
        std::cout << "\n";
    }

public:
    CountEstimator(int _n, int _ksize, bool _save_kmers, std::vector<long> *predef_hash_list, bool rev_comp) {
        this->n = _n;
        this->predef_hash_list = predef_hash_list;
        this->save_kmers = _save_kmers;
        this->ksize = _ksize;
        this->rev_comp = rev_comp;
        if (save_kmers) {
            this->sketch_kmers.resize(n);
        }
        sketch_counts.resize(n);
        sketch_hashes.resize(n);
        for (int i = 0; i < n; i++) {
            sketch_counts[i] = 0;
            sketch_hashes[i] = std::numeric_limits<long>::max();
        }
    }

    int get_ksize() {
        return ksize;
    }

    bool is_save_kmers() {
        return save_kmers;
    }

    std::vector<long> *get_predef_hash_list() {
        return predef_hash_list;
    }

    bool is_rev_comp() {
        return rev_comp;
    }

    int bin_search(std::vector<long> v, int left, int right, int x) {
        if (x < v[left]) {
            return left - 1;
        }
        if (x > v[right]) {
            return right + 1;
        }
        if (left > right) {
            return 0;
        }
        int mid = (left + right + 1) / 2;
        if (right == mid) {
            // consecutive number situation, or equal
            if (x <= v[left]) {
                return left;
            }
            return right;
        }
        if (v[mid] < x) {
            return bin_search(v, mid, right, x);
        } else {
            return bin_search(v, left, mid, x);
        }
    }

    void add(std::string kmer) {
        if (kmer.length() == 0) {
            return;
        }
        long h;
        if (rev_comp) {
            long h1 = hash(kmer);
            long h2 = hash(reverse_complement(kmer));
            h = std::min(h1, h2);
            if (h == h2) {
                kmer = reverse_complement(kmer);
            }
        } else {
            h = hash(kmer);
        }
        std::cout << "Hash is : " << h << " \n";
        if (predef_hash_list != 0) {
            if (std::find(predef_hash_list->begin(), predef_hash_list->end(), h) == predef_hash_list->end()) {
                /* predef hash list does not contain x */
                std::cout << "Returning as predef Hash list doesn't contain value.: " << h << " \n";
                return;
            }
        }
//        long p = get_prime_lt_x(max_prime);
        long p = max_prime;
        std::cout << "Prime is: " << p << " \n";
        h = h % p;
        int idx = bin_search(sketch_hashes, 0, n - 1, h);
        std::cout << "Idx is: " << idx << " \n";
        if (idx > n - 1) {
            // hash rank is greater than n
            return;
        }
        if (idx < 0) {
            // make this the first hash, shift everything else right by 1
            insert_vect_single_long(&sketch_hashes, 0, h);
            insert_vect_single_long(&sketch_counts, 0, 1);
            if (save_kmers) {
                insert_vect_single_str(&sketch_kmers, 0, kmer);
            }
            // check if size is greater than n. if yes, resize
            print_vec(&sketch_hashes);
            print_vec(&sketch_counts);
            resize_vects_if_reqd();
            print_sketch();
            return;
        }
        // now is the middle case
        // check if idx value actually matches h
        if (sketch_hashes[idx] == h) {
            sketch_counts[idx]++;
            return;
        } else {
            // insert here
            insert_vect_single_long(&sketch_hashes, idx, h);
            insert_vect_single_long(&sketch_counts, idx, 1);
            if (save_kmers) {
                insert_vect_single_str(&sketch_kmers, idx, kmer);
            }
            // check if size is greater than n. if yes, resize
            resize_vects_if_reqd();
        }
        print_sketch();
    }

    void print_sketch() {
        int len = sketch_counts.size();
        std::cout << "Sketch: \n";
        for (int i = 0; i < len; i++) {
            std::cout << "Kmer # " << i << " ; ";
            if (save_kmers) {
                std::cout << "Value : " << sketch_kmers[i] << " ;  ";
            }
            std::cout << "Hash : " << sketch_hashes[i] << " ;  ";
            std::cout << "Count : " << sketch_counts[i] << " ;  ";
            std::cout << "\n";
        }
    }


};