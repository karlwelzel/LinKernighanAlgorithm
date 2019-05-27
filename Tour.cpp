//
// Created by Karl Welzel on 19/04/2019.
//

#include <utility>
#include <limits>
#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
#include "Tour.h"


// =================================================== Tour class ======================================================

std::vector<dimension_t> BaseTour::inversePermutation(const std::vector<dimension_t> &permutation) {
    std::vector<dimension_t> result(permutation.size());
    for (dimension_t i = 0; i < permutation.size(); ++i) {
        result[permutation[i]] = i;
    }
    return result;
}

std::vector<vertex_t> BaseTour::getNeighbors(vertex_t vertex) {
    return std::vector<vertex_t>({predecessor(vertex), successor(vertex)});
}

bool BaseTour::containsEdge(vertex_t vertex1, vertex_t vertex2) const {
    return predecessor(vertex1) == vertex2 or successor(vertex1) == vertex2;
}

std::vector<dimension_t> BaseTour::cyclicPermutation(const std::vector<vertex_t> &alternatingWalk) const {
    // For every out-edge in the alternating walk choose the vertex that is encountered first on the tour
    // This cuts the number of vertices that need to be sorted in half
    std::vector<dimension_t> outEdges; // contains all indices of the chosen vertices
    vertex_t start = 0; // the vertex from which the order is determined
    for (dimension_t i = 0; i < alternatingWalk.size(); i += 2) {
        // start must be different from alternatingWalk[i] and alternatingWalk[i+1]
        while (start == alternatingWalk[i] or start == alternatingWalk[i + 1]) {
            start = (start + 1) % getDimension();
        }
        if (isBetween(start, alternatingWalk[i], alternatingWalk[i + 1])) {
            outEdges.push_back(i);
        } else {
            outEdges.push_back(i + 1);
        }
    }

    // Sort these vertices by the order the appear on the tour
    std::sort(outEdges.begin(), outEdges.end(),
              [this, alternatingWalk](dimension_t i, dimension_t j) {
                  return isBetween(0, alternatingWalk[i], alternatingWalk[j]);
              });

    // Add all other vertices to get a complete permutation
    std::vector<vertex_t> permutation;
    for (dimension_t i : outEdges) {
        permutation.push_back(i);
        permutation.push_back((i % 2 == 0) ? i + 1 : i - 1);
    }

    return permutation;
}

bool BaseTour::isTourAfterExchange(const std::vector<vertex_t> &alternatingWalk) const {
    std::vector<dimension_t> permutation = cyclicPermutation(alternatingWalk);
    std::vector<dimension_t> indices = inversePermutation(permutation);
    // {alternatingWalk[permutation[2i]], alternatingWalk[permutation[2i+1]]} are the out-edges in alternatingWalk

    dimension_t i = alternatingWalk.size() - 1;
    dimension_t j;
    dimension_t count = 0;
    do {
        j = ((permutation[i] % 2 == 0) ? permutation[i] - 1 : permutation[i] + 1) % alternatingWalk.size();
        // Walk along the in-edge {alternatingWalk[permutation[i]], alternatingWalk[j]}
        i = ((indices[j] % 2 == 0) ? indices[j] - 1 : indices[j] + 1) % alternatingWalk.size();
        // Walk along the segment alternatingWalk[j] -> alternatingWalk[permutation[i]]
        // because with k = (indices[j] % 2 == 0) ? indices[j] + 1 : indices[j] - 1
        // {alternatingWalk[j], alternatingWalk[permutation[k]]} would be an out-edge
        count++;
    } while (i != alternatingWalk.size() - 1);

    return 2 * count == alternatingWalk.size();
}

void BaseTour::exchange(const std::vector<vertex_t> &alternatingWalk) {
    // TODO: Implement an exchange step using Bergeron's algorithm and flip
}


// ================================================ ArrayTour class ====================================================

dimension_t ArrayTour::getDimension() const {
    return sequence.size();
}

void ArrayTour::setVertices(const std::vector<vertex_t> &vertexList) {
    sequence = vertexList;
    indices = inversePermutation(sequence);
}

ArrayTour::ArrayTour(const std::vector<vertex_t> &vertexList) {
    setVertices(vertexList);
}

vertex_t ArrayTour::predecessor(vertex_t vertex) const {
    return sequence[(indices[vertex] - 1) % getDimension()];
}

vertex_t ArrayTour::successor(vertex_t vertex) const {
    return sequence[(indices[vertex] + 1) % getDimension()];
}

bool ArrayTour::isBetween(vertex_t before, vertex_t vertex, vertex_t after) const {
    dimension_t distanceToVertex = (indices[vertex] - indices[before]) % getDimension();
    dimension_t distanceToAfter = (indices[after] - indices[before]) % getDimension();
    return distanceToVertex < distanceToAfter;
}

void ArrayTour::flip(vertex_t a, vertex_t b, vertex_t c, vertex_t d) {
    // Calculate the length of the two segments to decide which of them will be reversed
    dimension_t distanceSegment1 = (indices[c] - indices[a]) % getDimension();
    dimension_t distanceSegment2 = (indices[b] - indices[d]) % getDimension();
    vertex_t segmentStart, segmentEnd;
    dimension_t segmentLength;
    if (distanceSegment1 <= distanceSegment2) {
        // The segment a-c will be reversed
        segmentStart = a;
        segmentEnd = c;
        segmentLength = distanceSegment1;
    } else {
        // The segment d-b will be reversed
        segmentStart = d;
        segmentEnd = b;
        segmentLength = distanceSegment2;
    }
    for (dimension_t i = 0; i < segmentLength / 2; ++i) {
        // TODO: Check if this does only swap the contents of vertex1 and vertex2
        vertex_t &vertex1 = sequence[(segmentStart + i) % getDimension()];
        vertex_t &vertex2 = sequence[(segmentEnd - i) % getDimension()];
        std::swap(vertex1, vertex2);
        std::swap(indices[vertex1], indices[vertex2]);
    }

}


// =============================================== TourWalker class ====================================================


TourWalker::TourWalker(const BaseTour &tour, vertex_t first) : tour(tour), current(first) {}

vertex_t TourWalker::getNextVertex() {
    const vertex_t previous = current;
    current = tour.successor(current);
    return previous;
}

std::ostream &operator<<(std::ostream &out, const BaseTour &tour) {
    std::string output;
    TourWalker tourWalker(tour, 0);
    vertex_t currentVertex = tourWalker.getNextVertex();
    do {
        output += std::to_string(currentVertex + 1) + ", ";
        currentVertex = tourWalker.getNextVertex();
    } while (currentVertex != 0);
    out << output.substr(0, output.length() - 2);
    return out;
}
