//
// Created by karl on 27.05.19.
//

#ifndef LINKERNINGHANALGORITHM_SIGNEDPERMUTATION_H
#define LINKERNINGHANALGORITHM_SIGNEDPERMUTATION_H

#include <vector>

using number_t = long; // The data type used for numbers in signed permutations

// ============================================ SignedPermutation class ================================================

// This class represents a signed permutation and provides algorithms to sort it by reversals.
// The algorithms are described here: https://www.sciencedirect.com/science/article/pii/S0166218X04003440

class SignedPermutation {
private:
    std::vector<number_t> permutation;

public:
    // Initialize the signed permutation with permutation
    // Expects a vector of length n where for every 0 <= i < n either i or -i appears exactly once and where
    // the first element is 0 and the last is n-1
    explicit SignedPermutation(std::vector<number_t> permutation);

    // TODO: Implement Bergeron's algorithms and test them
};


#endif //LINKERNINGHANALGORITHM_SIGNEDPERMUTATION_H
