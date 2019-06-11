//
// Created by Karl Welzel on 19/04/2019.
//

#include <utility>
#include <limits>
#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include "Tour.h"
#include "SignedPermutation.h"


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
    vertex_t start = 0; // the vertex from which the order of each out-edge is determined
    for (dimension_t i = 0; i < alternatingWalk.size() - 1; i += 2) {
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

    dimension_t size = permutation.size();
    dimension_t i = permutation[1]; // index for alternatingWalk
    dimension_t j;                  // index for permutation
    dimension_t count = 0;

    // Traverse the in-edges on the tour and skip the segments in between
    // Count is a counter for the number of in-edges that are visited. The exchange produces a tour if and only if
    // count is half as big as the number of edges in alternatingWalk.
    do {
        j = ((indices[i] % 2 == 0) ? indices[i] + size - 1 : indices[i] + 1) % size;
        // alternatingWalk[i] -> alternatingWalk[permutation[j]] is a segment on the tour
        // because with k = (indices[i] % 2 == 0) ? indices[i] + 1 : indices[i] - 1
        // {alternatingWalk[i], alternatingWalk[permutation[k]]} would be an out-edge

        i = ((permutation[j] % 2 == 0) ? permutation[j] + size - 1 : permutation[j] + 1) % size;
        // {alternatingWalk[permutation[j]], alternatingWalk[i]} is an in-edge on the tour

        count++;
    } while (i != permutation[1]);

    return 2 * count == size;
}

void BaseTour::exchange(const std::vector<vertex_t> &alternatingWalk) {
    std::vector<dimension_t> permutation = cyclicPermutation(alternatingWalk);
    std::vector<dimension_t> indices = inversePermutation(permutation);
    // {alternatingWalk[permutation[2i]], alternatingWalk[permutation[2i+1]]} are the out-edges in alternatingWalk

    std::vector<std::pair<vertex_t, vertex_t>> segments;
    std::unordered_map<vertex_t, number_t> segmentNumber;
    std::vector<std::pair<number_t, bool>> segmentPermutation;

    dimension_t size = permutation.size();
    dimension_t i = permutation[1]; // index for alternatingWalk
    dimension_t j;                  // index for permutation

    // Traverse the segments on the new tour and number them
    do {
        j = ((indices[i] % 2 == 0) ? indices[i] + size - 1 : indices[i] + 1) % size;
        // alternatingWalk[i] -> alternatingWalk[permutation[j]] is a segment on the tour
        // because with k = (indices[i] % 2 == 0) ? indices[i] + 1 : indices[i] - 1
        // {alternatingWalk[i], alternatingWalk[permutation[k]]} would be an out-edge

        segments.emplace_back(alternatingWalk[i], alternatingWalk[permutation[j]]);
        segmentNumber[alternatingWalk[i]] = segments.size() - 1;
        // segmentNumber[segments[a].first] == a

        i = ((permutation[j] % 2 == 0) ? permutation[j] + size - 1 : permutation[j] + 1) % size;
        // {alternatingWalk[permutation[j]], alternatingWalk[i]} is an in-edge on the tour
    } while (i != permutation[1]);

    // Traverse the segments on the current tour and create a signed permutation
    for (j = 1; j < size; j += 2) {
        // alternatingWalk[permutation[j]] -> alternatingWalk[permutation[j+1]] is a segment on the tour
        std::unordered_map<vertex_t, number_t>::const_iterator it = segmentNumber.find(alternatingWalk[permutation[j]]);
        if (it != segmentNumber.end()) { // The segment on the new tour has successor direction
            segmentPermutation.emplace_back(it->second, true);
        } else { // The segment on the new tour has predecessor direction
            it = segmentNumber.find(alternatingWalk[permutation[(j + 1) % size]]);
            segmentPermutation.emplace_back(it->second, false);
        }
        // {alternatingWalk[permutation[j+1]], alternatingWalk[permutation[j+2]]} is an out-edge on the tour
    }

    SignedPermutation signedPermutation(segmentPermutation);
    dimension_t segmentsSize = segments.size();

    // Compute the reversal steps needed to transform segmentPermutation to the identity permutation and translate
    // these reversal steps to 2-opt exchanges

    // preStartVertex and postStartVertex are needed, because the reversal of the signed permutation does not have to
    // be in the successor direction on the tour
    std::pair<number_t, bool> preStartElement, startElement, endElement, postEndElement;
    std::pair<vertex_t, vertex_t> preStartSegment, startSegment, endSegment, postEndSegment;
    vertex_t preStartVertex, startVertex, endVertex, postEndVertex;
    while (!signedPermutation.isIdentityPermutation()) {
        std::pair<size_t, size_t> reversal = signedPermutation.nextReversal();
        preStartElement = signedPermutation.getElementAt((reversal.first + segmentsSize - 1) % segmentsSize);
        startElement = signedPermutation.getElementAt(reversal.first);
        endElement = signedPermutation.getElementAt(reversal.second);
        postEndElement = signedPermutation.getElementAt((reversal.second + 1) % segmentsSize);

        preStartSegment = segments[preStartElement.first];
        startSegment = segments[startElement.first];
        endSegment = segments[endElement.first];
        postEndSegment = segments[postEndElement.first];

        preStartVertex = preStartElement.second ? preStartSegment.second : preStartSegment.first;
        startVertex = startElement.second ? startSegment.first : startSegment.second;
        endVertex = endElement.second ? endSegment.second : endSegment.first;
        postEndVertex = postEndElement.second ? postEndSegment.first : postEndSegment.second;

        if (preStartVertex != predecessor(startVertex)) {
            std::swap(preStartVertex, startVertex);
        }
        if (postEndVertex != successor(endVertex)) {
            std::swap(postEndVertex, endVertex);
        }
        flip(startVertex, preStartVertex, endVertex, postEndVertex);

        signedPermutation.performReversal(reversal);
    }
}


// ================================================ ArrayTour class ====================================================

dimension_t ArrayTour::getDimension() const {
    return sequence.size();
}

dimension_t ArrayTour::distance(vertex_t vertex1, vertex_t vertex2) const {
    // Avoid problems with negative values and modulo
    ptrdiff_t diff = indices[vertex2] - indices[vertex1];
    return (diff + getDimension()) % getDimension();
}

void ArrayTour::setVertices(const std::vector<vertex_t> &vertexList) {
    sequence = vertexList;
    indices = inversePermutation(sequence);
}

ArrayTour::ArrayTour(const std::vector<vertex_t> &vertexList) {
    setVertices(vertexList);
}

vertex_t ArrayTour::predecessor(vertex_t vertex) const {
    return sequence[(indices[vertex] + getDimension() - 1) % getDimension()];
}

vertex_t ArrayTour::successor(vertex_t vertex) const {
    return sequence[(indices[vertex] + 1) % getDimension()];
}

bool ArrayTour::isBetween(vertex_t before, vertex_t vertex, vertex_t after) const {
    dimension_t distanceToVertex = distance(before, vertex);
    dimension_t distanceToAfter = distance(before, after);
    return distanceToVertex < distanceToAfter;
}

void ArrayTour::flip(vertex_t a, vertex_t b, vertex_t c, vertex_t d) {
    // Calculate the length of the two segments to decide which of them will be reversed
    dimension_t distanceSegment1 = distance(a, c);
    dimension_t distanceSegment2 = distance(d, b);
    dimension_t segmentStartIndex, segmentEndIndex;
    dimension_t segmentLength;
    if (distanceSegment1 <= distanceSegment2) {
        // The segment a-c will be reversed
        segmentStartIndex = indices[a];
        segmentEndIndex = indices[c];
        segmentLength = distanceSegment1;
    } else {
        // The segment d-b will be reversed
        segmentStartIndex = indices[d];
        segmentEndIndex = indices[b];
        segmentLength = distanceSegment2;
    }
    for (dimension_t i = 0; i < (segmentLength + 1) / 2; ++i) {
        vertex_t &vertex1 = sequence[(segmentStartIndex + i) % getDimension()];
        vertex_t &vertex2 = sequence[(segmentEndIndex + getDimension() - i) % getDimension()];
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
