//
// Created by Karl Welzel on 05/27/19.
//

#include <utility>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "SignedPermutation.h"

SignedPermutation::SignedPermutation(std::vector<std::pair<number_t, bool>> permutation) : permutation(
        std::move(permutation)) {
    // DEBUG: Check whether the permutation is valid
    std::vector<bool> found(this->permutation.size(), false);
    for (std::pair<number_t, bool> element : this->permutation) {
        found[element.first] = true;
    }
    if (!std::all_of(found.begin(), found.end(), [](bool v) { return v; })) {
        throw std::runtime_error("permutation is not a permutation");
    }

    // Calculate the index of each element
    indices.assign(this->permutation.size(), 0);
    for (size_t i = 0; i < this->permutation.size(); ++i) {
        indices[this->permutation[i].first] = i;
    }
}

std::pair<number_t, bool> SignedPermutation::getElementAt(size_t i) const {
    return permutation[i];
}

bool SignedPermutation::isIdentityPermutation() const {
    for (size_t i = 0; i < permutation.size(); ++i) {
        // Check that the i-th element is +i
        if (permutation[i].first != i or !permutation[i].second) {
            return false;
        }
    }
    return true;
}

std::pair<size_t, size_t> SignedPermutation::nextReversal() const {
    for (size_t i = 0; i < permutation.size(); ++i) {
        if (permutation[i].first != i) {
            // Return reversal that reverses permutation[i] and i
            size_t j = indices[i];
            return std::make_pair(i, j);
        }
        if (!permutation[i].second) {
            // Return reversal that reverses i to change the sign (permutation[i] = i)
            return std::make_pair(i, i);
        }
    }

    throw std::runtime_error("The current signed permutation already is the identity permutation.");
}


void SignedPermutation::performReversal(std::pair<size_t, size_t> step) {
    number_t i = step.first;
    number_t j = step.second;
    for (size_t k = 0; k < (j - i + 1) / 2; ++k) {
        std::swap(permutation[i + k], permutation[j - k]);
        std::swap(indices[permutation[i + k].first], indices[permutation[j - k].first]);
        permutation[i + k].second = !permutation[i + k].second;
        permutation[j - k].second = !permutation[j - k].second;
    }
    if ((j - i) % 2 == 0) {
        // Change sign for the element in the middle, that does not change its position
        permutation[(i + j) / 2].second = !permutation[(i + j) / 2].second;
    }

}

