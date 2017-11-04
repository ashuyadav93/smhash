//
// Created by Saraj Munjal on 11/3/17.
//
#include <utility>

class CountEstimator {
private:
    const long max_prime = 9999999999971;
    int ksize;
    int n;
    //FILE* input_file_name=None,
    bool save_kmers = false;
    std::vector<long> *predef_hash_list = NULL;
    bool rev_comp = false;
    std::vector<long> sketch_counts;
    std::vector<long> sketch_hashes;
    std::vector <std::pair<int, int> > sketch_kmers;
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
            sketch_hashes[i] = std::numeric_limits<long>::min();
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
            return left -1;
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
};