//
// Created by Karl Welzel on 27.05.19.
//

#ifndef LINKERNINGHANALGORITHM_SIGNEDPERMUTATION_H
#define LINKERNINGHANALGORITHM_SIGNEDPERMUTATION_H

#include <vector>

using number_t = size_t; // The data type used for numbers in permutations

// ============================================ SignedPermutation class ================================================

// This class represents a signed permutation and provides algorithms to sort it by reversals.

// A signed permutation of length n is a permutation of the numbers 0 to n-1 where each number has a sign (+ or -). It
// is represented by an std::pair<number_t, bool>, where the the first element is the number 0 to n-1 and the second
// elements represents the sign (true = +, false = -)
// A signed reversal reverses the order of consecutive elements in the permutation and changes their signs. The goal is
// to find a number of signed reversals that transform a given signed permutation into the signed identity permutation
// +0, +1, +2, ..., +n. A signed reversal is represented by a pair of indices i <= j that indicate that the part of
// permutation between i and j should be reversed.

// The algorithm implemented here is an naive greedy approach:
// The permutation is sorted in a way such that after the i-th step the i first elements are +1 to +i in the correct
// order. To do this the reversal (permutation[i], i) is applied (disregarding the sign) and then the reversal (i, i)
// is used to make the sign of i positive, if necessary. Obviously at worst 2n reversals are used and the number of
// reversal steps is not minimal.

class SignedPermutation {
private:
    std::vector<std::pair<number_t, bool>> permutation;
    std::vector<size_t> indices;

public:
    // Initialize the signed permutation with permutation
    // Expects a vector of length n where for every 0 <= i < n either (i, true) or (i, false) appears exactly once
    explicit SignedPermutation(std::vector<std::pair<number_t, bool>> permutation);

    // Returns the i-th element of the signed permutation
    std::pair<number_t, bool> getElementAt(size_t i) const;

    // Checks whether this signed permutation is the identity permutation, i.e. +0, +1, ..., +n
    bool isIdentityPermutation() const;

    // Computes the next reversal to transform this signed permutation into the identity permutation
    // After a finite number of performReversal(nextReversal()) the permutation will become the identity permutation
    std::pair<size_t, size_t> nextReversal() const;

    // Performs a signed reversal on the signed permutation
    // Expects step.first <= step.second
    void performReversal(std::pair<size_t, size_t> step);
};


#endif //LINKERNINGHANALGORITHM_SIGNEDPERMUTATION_H
