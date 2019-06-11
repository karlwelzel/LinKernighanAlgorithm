//
// Created by Karl Welzel on 05/27/19.
//

#include <utility>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "SignedPermutation.h"

SignedPermutation::SignedPermutation(std::vector<std::pair<number_t, bool>> permutation)
        : permutation(std::move(permutation)) {
    // DEBUG: Check whether permutation is valid
    std::vector<bool> found(this->permutation.size(), false);
    for (std::pair<number_t, bool> element : this->permutation) {
        found[element.first] = true;
    }

    if (!std::all_of(found.begin(), found.end(), [](bool v) { return v; })) {
        throw std::runtime_error("permutation is not a permutation");
    }

}

std::vector<std::pair<std::pair<number_t, bool>, std::pair<number_t, bool>>> SignedPermutation::reversalSteps() {
    std::vector<std::pair<number_t, bool>> currentPermutation = permutation;
    std::vector<std::pair<std::pair<number_t, bool>, std::pair<number_t, bool>>> reversals;

    // Calculate the index of each element
    std::vector<number_t> indices(currentPermutation.size());
    for (size_t i = 0; i < currentPermutation.size(); ++i) {
        indices[currentPermutation[i].first] = i;
    }

    for (size_t i = 0; i < currentPermutation.size(); ++i) {
        // Perform reversal (permutation[i], i) if necessary
        if (currentPermutation[i].first != i) {
            size_t j = indices[i];
            reversals.emplace_back(currentPermutation[i], currentPermutation[j]);
            for (size_t k = 0; k < (j - i + 1) / 2; ++k) {
                std::swap(currentPermutation[i + k], currentPermutation[j - k]);
                std::swap(indices[currentPermutation[i + k].first], indices[currentPermutation[j - k].first]);
                currentPermutation[i + k].second = !currentPermutation[i + k].second;
                currentPermutation[j - k].second = !currentPermutation[j - k].second;
            }
            if ((j - i) % 2 == 0) {
                // Change sign for the element in the middle, that does not change its position
                currentPermutation[(i + j) / 2].second = !currentPermutation[(i + j) / 2].second;
            }
        }
        // Perform reversal (i, i) if necessary
        if (!currentPermutation[i].second) {
            reversals.emplace_back(currentPermutation[i], currentPermutation[i]);
            currentPermutation[i].second = !currentPermutation[i].second;
        }
    }
    return reversals;
}
