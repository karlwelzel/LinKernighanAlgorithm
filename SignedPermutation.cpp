#include <utility>
#include <cmath>
#include <stdexcept>
#include <algorithm>

//
// Created by karl on 27.05.19.
//

#include "SignedPermutation.h"

SignedPermutation::SignedPermutation(std::vector<number_t> permutation) : permutation(std::move(permutation)) {
    // Check whether permutation is valid
    if (permutation.at(0) != 0 or permutation.at(permutation.size() - 1) != (number_t) permutation.size() - 1) {
        throw std::runtime_error("The given permutation does not have the correct format.");
    }

    std::vector<bool> found(permutation.size(), false);
    for (number_t element : permutation) {
        found.at(std::abs(element)) = true;
    }

    if (!std::all_of(found.begin(), found.end(), [](bool v) { return v; })) {
        throw std::runtime_error("The given permutation does not have the correct format.");
    }
}
