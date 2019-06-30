//
// Created by Karl Welzel on 27.05.19.
//

#include <algorithm>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>
#include "SignedPermutation.h"

SignedPermutation::SignedPermutation(std::vector<std::pair<number_t, bool>> permutation) :
        permutation(std::move(permutation)) {
    // Calculate the index of each element
    indices.resize(this->permutation.size());
    for (std::size_t i = 0; i < this->permutation.size(); ++i) {
        indices[this->permutation[i].first] = i;
    }
}

std::pair<number_t, bool> SignedPermutation::getElementAt(std::size_t i) const {
    return permutation[i];
}

bool SignedPermutation::isIdentityPermutation() const {
    for (std::size_t i = 0; i < permutation.size(); ++i) {
        // Check that the i-th element is +i
        if (permutation[i] != std::make_pair(i, true)) {
            return false;
        }
    }
    return true;
}

std::pair<std::size_t, std::size_t> SignedPermutation::nextReversal() const {
    for (std::size_t i = 0; i < permutation.size(); ++i) {
        if (permutation[i].first != i) {
            // Return reversal that swaps permutation[i] and i
            // This corresponds to the reversal (indices[permutation[i]], indices[i]) = (i, indices[i])
            return std::make_pair(i, indices[i]);
        }
        if (!permutation[i].second) {
            // Return reversal that reverses only i to change the sign
            // This corresponds to the reversal (indices[i], indices[i]) = (i, i)
            return std::make_pair(i, i);
        }
    }

    throw std::runtime_error("The current signed permutation already is the identity permutation.");
}


void SignedPermutation::performReversal(std::pair<std::size_t, std::size_t> step) {
    number_t i = step.first;
    number_t j = step.second;
    for (std::size_t k = 0; k < (j - i + 1) / 2; ++k) {
        // Swap elements and update indices
        std::swap(permutation[i + k], permutation[j - k]);
        std::swap(indices[permutation[i + k].first], indices[permutation[j - k].first]);

        // Change the signs
        permutation[i + k].second = !permutation[i + k].second;
        permutation[j - k].second = !permutation[j - k].second;
    }
    if ((j - i) % 2 == 0) {
        // Change sign of the element in the middle that does not change its position
        permutation[(i + j) / 2].second = !permutation[(i + j) / 2].second;
    }

}

