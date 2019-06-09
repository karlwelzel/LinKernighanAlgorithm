//
// Created by Karl Welzel on 05/27/19.
//

#include <utility>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "SignedPermutation.h"

SignedPermutation::SignedPermutation(std::vector<number_t> permutation, std::vector<bool> positiveSign)
        : permutation(std::move(permutation)), positiveSign(std::move(positiveSign)) {
    // DEBUG: Check whether permutation is valid
    if (this->permutation.size() != this->positiveSign.size()) {
        throw std::runtime_error("permutation and positiveSign do not have the same size");
    }

    std::vector<bool> found(this->permutation.size(), false);
    for (number_t element : this->permutation) {
        found.at(element) = true;
    }

    if (!std::all_of(found.begin(), found.end(), [](bool v) { return v; })) {
        throw std::runtime_error("permutation is not a permutation");
    }

}

std::vector<std::pair<number_t, number_t>> SignedPermutation::computeReversalSteps() {
    std::vector<number_t> currentPermutation = permutation;
    std::vector<bool> currentPositiveSign = positiveSign;
    std::vector<std::pair<number_t, number_t>> reversalSteps;

    // Calculate the index of each element
    std::vector<number_t> indices(currentPermutation.size());
    for (size_t i = 0; i < currentPermutation.size(); ++i) {
        indices[currentPermutation[i]] = i;
    }

    for (size_t i = 0; i < currentPermutation.size(); ++i) {
        // Perform reversal (permutation[i], i) if necessary
        if (currentPermutation[i] != i) {
            reversalSteps.emplace_back(currentPermutation[i], i);
            size_t j = indices[i];
            for (size_t k = 0; k < (j - i + 1) / 2; ++k) {
                std::swap(currentPermutation[i + k], currentPermutation[j - k]);
                std::swap(indices[currentPermutation[i + k]], indices[currentPermutation[j - k]]);
                std::swap(currentPositiveSign[i + k], currentPositiveSign[j - k]);
                currentPositiveSign[i + k] = !currentPositiveSign[i + k];
                currentPositiveSign[j - k] = !currentPositiveSign[j - k];
            }
            if ((j - i) % 2 == 0) {
                // Change sign for the element in the middle, that does not change its position
                currentPositiveSign[(i + j) / 2] = !currentPositiveSign[(i + j) / 2];
            }
        }
        // Perform reversal (i, i) if necessary
        if (!currentPositiveSign[i]) {
            reversalSteps.emplace_back(i, i);
            currentPositiveSign[i] = !currentPositiveSign[i];
        }
    }
    return reversalSteps;
}
