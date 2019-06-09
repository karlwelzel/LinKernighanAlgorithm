//
// Created by Karl Welzel on 05/27/19.
//

#ifndef LINKERNINGHANALGORITHM_SIGNEDPERMUTATION_H
#define LINKERNINGHANALGORITHM_SIGNEDPERMUTATION_H

#include <vector>

using number_t = size_t; // The data type used for numbers in permutations

// ============================================ SignedPermutation class ================================================

// This class represents a signed permutation and provides algorithms to sort it by reversals.
// The algorithm implemented here is an naive greedy approach:
// The permutation is sorted in a way such that after the i-th step the i-th first elements are 1 to i in the correct
// order. To do this the reversal (permutation[i], i) is applied (disregarding the sign) and then the reversal (i, i)
// is used to make the sign of i positive, if necessary. Obviously at worst 2n reversals are used and the algorithm
// number of reversal steps is not minimal.

// A better algorithm is described here:
// https://www.sciencedirect.com/science/article/pii/S0166218X04003440

class SignedPermutation {
private:
    std::vector<number_t> permutation;
    std::vector<bool> positiveSign;

public:
    // Initialize the signed permutation with permutation
    // Expects a vector of length n where for every 0 <= i < n either i or -i appears exactly once
    explicit SignedPermutation(std::vector<number_t> permutation, std::vector<bool> positiveSign);

    // Computes reversal steps so that when these signed reversals are applied in this order the permutation will
    // become the identity permutation
    std::vector<std::pair<number_t, number_t>> computeReversalSteps();
};


#endif //LINKERNINGHANALGORITHM_SIGNEDPERMUTATION_H
