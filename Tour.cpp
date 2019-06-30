//
// Created by Karl Welzel on 19.04.2019.
//

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <list>
#include <unordered_map>
#include <utility>
#include <vector>
#include "Tour.h"
#include "SignedPermutation.h"


// ============================================= AlternatingWalk class =================================================

std::size_t AlternatingWalk::size() const noexcept {
    return vertices.size();
}

vertex_t AlternatingWalk::operator[](std::size_t index) const {
    return vertices[index];
}

std::vector<vertex_t>::const_iterator AlternatingWalk::begin() const noexcept {
    return vertices.begin();
}

std::vector<vertex_t>::const_iterator AlternatingWalk::end() const noexcept {
    return vertices.end();
}

void AlternatingWalk::clear() noexcept {
    vertices.clear();
}

std::vector<vertex_t>::const_iterator AlternatingWalk::erase(std::vector<vertex_t>::const_iterator first,
                                                             std::vector<vertex_t>::const_iterator last) {
    return vertices.erase(first, last);
}

void AlternatingWalk::push_back(vertex_t vertex) {
    vertices.push_back(vertex);
}

AlternatingWalk AlternatingWalk::close() const {
    AlternatingWalk result(*this);
    result.push_back(vertices[0]); // Add the first vertex to the end to close the walk
    return result;
}

AlternatingWalk AlternatingWalk::appendAndClose(vertex_t vertex) const {
    AlternatingWalk result(*this);
    result.push_back(vertex); // Add vertex to the end
    result.push_back(vertices[0]); // Close the walk as in close
    return result;
}

bool AlternatingWalk::containsEdge(vertex_t v, vertex_t w) const {
    for (dimension_t i = 0; i < vertices.size() - 1; ++i) {
        // Check if the edge {vertices[i], vertices[i+1]} is the same as (v, w)
        if ((vertices[i] == v and vertices[i + 1] == w) or (vertices[i] == w and vertices[i + 1] == v)) {
            return true;
        }
    }
    return false;
}


// ================================================= BaseTour class ====================================================

std::vector<dimension_t> BaseTour::inversePermutation(const std::vector<dimension_t> &permutation) {
    std::vector<dimension_t> result(permutation.size());
    for (dimension_t i = 0; i < permutation.size(); ++i) {
        result[permutation[i]] = i;
    }
    return result;
}

std::vector<vertex_t> BaseTour::getNeighbors(vertex_t vertex) {
    return {predecessor(vertex), successor(vertex)};
}

bool BaseTour::containsEdge(vertex_t vertex1, vertex_t vertex2) const {
    return predecessor(vertex1) == vertex2 or successor(vertex1) == vertex2;
}

std::vector<dimension_t> BaseTour::cyclicPermutation(const AlternatingWalk &alternatingWalk) const {
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
    // The vertex 0 is chosen as an arbitrary starting point on the tour
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

bool BaseTour::isTourAfterExchange(const AlternatingWalk &alternatingWalk) const {
    // Compute the cyclic permutation and its inverse
    std::vector<dimension_t> permutation = cyclicPermutation(alternatingWalk);
    std::vector<dimension_t> indices = inversePermutation(permutation);

    dimension_t size = permutation.size();
    dimension_t i = permutation[1]; // index for alternatingWalk
    dimension_t j;                  // index for permutation
    dimension_t count = 0;

    // Traverse the in-edges of alternatingWalk and skip the segments of the tour in between the in-edges to simulate
    // a traversal of the tour after exchanging out-edges by in-edges
    // count is a counter for the number of in-edges that are visited. The exchange produces a tour if and only if
    // count is half as big as the number of edges in alternatingWalk.

    // {alternatingWalk[permutation[2i]], alternatingWalk[permutation[2i+1]]} are the out-edges in alternatingWalk
    do {
        j = ((indices[i] % 2 == 0) ? indices[i] + size - 1 : indices[i] + 1) % size;
        // alternatingWalk[i] -> alternatingWalk[permutation[j]] is a segment on the tour

        i = ((permutation[j] % 2 == 0) ? permutation[j] + size - 1 : permutation[j] + 1) % size;
        // {alternatingWalk[permutation[j]], alternatingWalk[i]} is an in-edge

        count++;
    } while (i != permutation[1]);

    return 2 * count == size;
}

void BaseTour::exchange(const AlternatingWalk &alternatingWalk) {
    // The "new" tour is the tour after exchanging all out-edges by in-edges, while the "old" tour is the current one
    // before this exchange.

    // Compute the cyclic permutation and its inverse
    std::vector<dimension_t> permutation = cyclicPermutation(alternatingWalk);
    std::vector<dimension_t> indices = inversePermutation(permutation);
    // {alternatingWalk[permutation[2i]], alternatingWalk[permutation[2i+1]]} are the out-edges in alternatingWalk

    // The completely unchanged segments v -> w of the tour between the out-edges as pairs (v, w)
    std::vector<std::pair<vertex_t, vertex_t>> segments;

    // A signed permutation corresponding to reordering and reversing of segments by the exchange
    // The new tour corresponds to the identity permutation and the current tour to segmentPermutation
    std::vector<std::pair<number_t, bool>> segmentPermutation;

    // Maps the first vertex of each such segment to its corresponding number in segmentPermutation
    std::unordered_map<vertex_t, number_t> segmentNumber;

    dimension_t size = permutation.size();
    dimension_t i = permutation[1]; // index for alternatingWalk
    dimension_t j;                  // index for permutation

    // Traverse the segments in the order they appear on the new tour (as in isTourAfterExchange) and number them
    do {
        j = ((indices[i] % 2 == 0) ? indices[i] + size - 1 : indices[i] + 1) % size;
        // alternatingWalk[i] -> alternatingWalk[permutation[j]] is a segment on the tour

        segments.emplace_back(alternatingWalk[i], alternatingWalk[permutation[j]]);
        segmentNumber[alternatingWalk[i]] = segments.size() - 1;
        // segmentNumber[segments[a].first] == a

        i = ((permutation[j] % 2 == 0) ? permutation[j] + size - 1 : permutation[j] + 1) % size;
        // {alternatingWalk[permutation[j]], alternatingWalk[i]} is an in-edge on the tour
    } while (i != permutation[1]);

    // Traverse the segments on the current tour and create segmentPermutation
    for (j = 1; j < size; j += 2) {
        // alternatingWalk[permutation[j]] -> alternatingWalk[permutation[j+1]] is a segment on the tour
        // Now we need to check the direction of the segment and append its number to segmentPermutation
        std::unordered_map<vertex_t, number_t>::const_iterator it = segmentNumber.find(alternatingWalk[permutation[j]]);
        if (it != segmentNumber.end()) { // The segment on the new tour has successor direction
            segmentPermutation.emplace_back(it->second, true);
        } else { // The segment on the new tour has predecessor direction
            it = segmentNumber.find(alternatingWalk[permutation[(j + 1) % size]]);
            segmentPermutation.emplace_back(it->second, false);
        }

        // {alternatingWalk[permutation[j+1]], alternatingWalk[permutation[j+2]]} is an out-edge on the tour
    }

    // Construct a SignedPermutation from segmentPermutation
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
        // Compute the next reversal step
        std::pair<std::size_t, std::size_t> reversal = signedPermutation.nextReversal();

        // The reversal reverses the elements with indices in [reversal.first, reversal.last]
        // Get the bounding elements of this reversal:
        //     preStartElement, [startElement, ..., endElement], postElement
        // where [ and ] indicate the reversal range
        preStartElement = signedPermutation.getElementAt((reversal.first + segmentsSize - 1) % segmentsSize);
        startElement = signedPermutation.getElementAt(reversal.first);
        endElement = signedPermutation.getElementAt(reversal.second);
        postEndElement = signedPermutation.getElementAt((reversal.second + 1) % segmentsSize);

        // Get the corresponding bounding segments of this reversal
        preStartSegment = segments[preStartElement.first];
        startSegment = segments[startElement.first];
        endSegment = segments[endElement.first];
        postEndSegment = segments[postEndElement.first];

        // Get the bounding vertices of this reversal which correspond to the start or end of a segment
        // The edges {preStartVertex, startVertex} and {endVertex, postEndVertex} should be replaced by two edges that
        // connect these vertices and form a tour. In this new tour the vertices between startVertex and endVertex
        // or the other half of the tour is reversed and this matches the reversal in the signed permutation.
        preStartVertex = preStartElement.second ? preStartSegment.second : preStartSegment.first;
        startVertex = startElement.second ? startSegment.first : startSegment.second;
        endVertex = endElement.second ? endSegment.second : endSegment.first;
        postEndVertex = postEndElement.second ? postEndSegment.first : postEndSegment.second;

        // Swap vertices if necessary to conform with the expectation of flip
        if (preStartVertex != predecessor(startVertex)) {
            std::swap(preStartVertex, startVertex);
        }
        if (postEndVertex != successor(endVertex)) {
            std::swap(postEndVertex, endVertex);
        }

        flip(startVertex, preStartVertex, endVertex, postEndVertex);

        // Perform the reversal on the signed permutation to reflect the flip.
        signedPermutation.performReversal(reversal);
    }
}


// ================================================ ArrayTour class ====================================================

dimension_t ArrayTour::getDimension() const {
    return sequence.size();
}

dimension_t ArrayTour::distance(vertex_t vertex1, vertex_t vertex2) const {
    // Avoid problems with negative values and modulo
    std::ptrdiff_t diff = indices[vertex2] - indices[vertex1];
    return (diff + getDimension()) % getDimension();
}

void ArrayTour::setVertices(const std::vector<vertex_t> &vertexList) {
    sequence = vertexList;
    indices = inversePermutation(sequence);
}

ArrayTour::ArrayTour(const std::vector<vertex_t> &tourSequence) {
    setVertices(tourSequence);
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
    dimension_t acDistance = distance(a, c);
    dimension_t dbDistance = distance(d, b);

    dimension_t segmentStartIndex, segmentEndIndex;
    dimension_t segmentLength;
    if (acDistance <= dbDistance) {
        // The segment a-c will be reversed
        segmentStartIndex = indices[a];
        segmentEndIndex = indices[c];
        segmentLength = acDistance;
    } else {
        // The segment d-b will be reversed
        segmentStartIndex = indices[d];
        segmentEndIndex = indices[b];
        segmentLength = dbDistance;
    }
    for (dimension_t i = 0; i < (segmentLength + 1) / 2; ++i) {
        vertex_t &vertex1 = sequence[(segmentStartIndex + i) % getDimension()];
        vertex_t &vertex2 = sequence[(segmentEndIndex + getDimension() - i) % getDimension()];
        std::swap(vertex1, vertex2);
        std::swap(indices[vertex1], indices[vertex2]);
    }
}


// ============================================ TwoLevelTreeTour class =================================================

//             ======================= TwoLevelTreeTour::SegmentParent class ===============================

bool operator==(const TwoLevelTreeTour::SegmentParent &parent, const TwoLevelTreeTour::SegmentParent &otherParent) {
    return parent.sequenceNumber == otherParent.sequenceNumber;
}

dimension_t TwoLevelTreeTour::SegmentParent::count() const {
    return vertices.size();
}

vertex_t TwoLevelTreeTour::SegmentParent::firstVertex() const {
    return reversed ? vertices.back().vertex : vertices.front().vertex;
}

vertex_t TwoLevelTreeTour::SegmentParent::lastVertex() const {
    return reversed ? vertices.front().vertex : vertices.back().vertex;
}

void TwoLevelTreeTour::SegmentParent::reverseVertices(std::list<SegmentVertex>::iterator first,
                                                      std::list<SegmentVertex>::iterator last) {
    auto insertIterator = std::next(last);

    // Move the elements to a separate list, reverse this list and move them back to vertices
    std::list<SegmentVertex> temporaryList{};
    temporaryList.splice(temporaryList.begin(), vertices, first, std::next(last));
    temporaryList.reverse();
    vertices.splice(insertIterator, temporaryList);

    // Reverse the sequenceNumbers, so that they stay consecutive within the segment
    std::swap(first, last);
    last++;
    while ((first != last) && (first != --last)) {
        std::swap((first++)->sequenceNumber, last->sequenceNumber);
    }
}

//             ============================== TwoLevelTreeTour class =======================================

TwoLevelTreeTour::TwoLevelTreeTour(const TwoLevelTreeTour &otherTour) : dimension(otherTour.dimension),
                                                                        groupSize(otherTour.groupSize),
                                                                        parents(otherTour.parents),
                                                                        iterators(otherTour.iterators) {

    // Update all iterators because they still point to elements in otherTour
    for (auto parentIterator = parents.begin(); parentIterator != parents.end(); ++parentIterator) {
        std::list<SegmentVertex> &vertices = parentIterator->vertices;
        for (auto vertexIterator = vertices.begin(); vertexIterator != vertices.end(); ++vertexIterator) {
            vertexIterator->parentIterator = parentIterator;
            iterators[vertexIterator->vertex] = vertexIterator;
        }
    }
}

TwoLevelTreeTour &TwoLevelTreeTour::operator=(const TwoLevelTreeTour &otherTour) {
    if (this != &otherTour) {
        dimension = otherTour.dimension;
        groupSize = otherTour.groupSize;
        parents = otherTour.parents;
        iterators = otherTour.iterators;

        // Update all iterators because they still point to elements in otherTour
        for (auto parentIterator = parents.begin(); parentIterator != parents.end(); ++parentIterator) {
            std::list<SegmentVertex> &vertices = parentIterator->vertices;
            for (auto vertexIterator = vertices.begin(); vertexIterator != vertices.end(); ++vertexIterator) {
                vertexIterator->parentIterator = parentIterator;
                iterators[vertexIterator->vertex] = vertexIterator;
            }
        }
    }

    return *this;
}


void TwoLevelTreeTour::setVertices(const std::vector<vertex_t> &tourSequence) {
    dimension = tourSequence.size();

    // For two segments the implementation does not work and for one segment (3/4)*groupSize may not be smaller
    // than the dimension, so for dimension <= 300 I set the groupSize to twice the dimension. This sets the number of
    // segments to one and the TwoLevelTreeTour essentially works like an ArrayTour

    // The choices of groupSize for all other cases is taken from the paper linked in Tour.h

    if (tourSequence.size() <= 300) {
        groupSize = 2 * dimension;
    } else if (tourSequence.size() <= 100000) {
        groupSize = 100;
    } else {
        groupSize = 200;
    }

    iterators.resize(tourSequence.size());
    dimension_t segmentLength = groupSize; // The first segments should have groupSize elements
    std::size_t parentIndex = 0; // The index for tourSequence
    std::size_t vertexIndex = 0; // The index for parents

    while (vertexIndex < dimension) {
        // The last segment should have a length between groupSize and 2*groupSize
        if (tourSequence.size() - vertexIndex < 2 * groupSize) {
            segmentLength = tourSequence.size() - vertexIndex;
        }

        // Create an empty SegmentParent and get a pointer to it
        parents.emplace_back();
        auto parentIterator = std::prev(parents.end());
        SegmentParent &parent = *parentIterator;

        // Initialize the parent values and fill parent.vertices
        parent.reversed = false;
        parent.sequenceNumber = parentIndex;
        for (std::size_t i = vertexIndex; i < vertexIndex + segmentLength; ++i) {
            // Create a new SegmentVertex and append it to parent.vertices
            parent.vertices.push_back(SegmentVertex{tourSequence[i], parentIterator, static_cast<long>(i)});
            iterators[tourSequence[i]] = std::prev(parent.vertices.end());
        }

        parentIndex++;
        vertexIndex += segmentLength;
    }
}

TwoLevelTreeTour::TwoLevelTreeTour(const std::vector<vertex_t> &tourSequence) {
    setVertices(tourSequence);
}


dimension_t TwoLevelTreeTour::getDimension() const {
    return dimension;
}

const TwoLevelTreeTour::SegmentParent &TwoLevelTreeTour::getPreviousParent(
        std::list<SegmentParent>::iterator parentIterator) const {
    return parentIterator != parents.begin() ? *std::prev(parentIterator) : parents.back();
}

const TwoLevelTreeTour::SegmentParent &TwoLevelTreeTour::getNextParent(
        std::list<SegmentParent>::iterator parentIterator) const {
    return std::next(parentIterator) != parents.end() ? *std::next(parentIterator) : parents.front();
}

vertex_t TwoLevelTreeTour::predecessor(vertex_t vertex) const {
    auto iterator = iterators[vertex];
    const SegmentParent &parent = *(iterator->parentIterator);
    if ((!parent.reversed and iterator == parent.vertices.begin()) or
        (parent.reversed and std::next(iterator) == parent.vertices.end())) {
        return getPreviousParent(iterator->parentIterator).lastVertex();
    } else {
        return parent.reversed ? std::next(iterator)->vertex : std::prev(iterator)->vertex;
    }
}

vertex_t TwoLevelTreeTour::successor(vertex_t vertex) const {
    auto iterator = iterators[vertex];
    const SegmentParent &parent = *(iterator->parentIterator);
    if ((!parent.reversed and std::next(iterator) == parent.vertices.end()) or
        (parent.reversed and iterator == parent.vertices.begin())) {
        return getNextParent(iterator->parentIterator).firstVertex();
    } else {
        return parent.reversed ? std::prev(iterator)->vertex : std::next(iterator)->vertex;
    }
}

bool TwoLevelTreeTour::isBetween(vertex_t before, vertex_t vertex, vertex_t after) const {
    auto beforeIterator = iterators[before];
    auto vertexIterator = iterators[vertex];
    auto afterIterator = iterators[after];
    const SegmentParent &beforeParent = *(beforeIterator->parentIterator);
    const SegmentParent &vertexParent = *(vertexIterator->parentIterator);
    const SegmentParent &afterParent = *(afterIterator->parentIterator);
    if (beforeParent == vertexParent and vertexParent == afterParent) {
        // All three vertices are in the same segment, so only the relative positions in the segment matter
        const SegmentParent &parent = vertexParent;
        long rawDistanceToVertex =
                (parent.reversed ? -1 : 1) * (vertexIterator->sequenceNumber - beforeIterator->sequenceNumber);
        long rawDistanceToAfter =
                (parent.reversed ? -1 : 1) * (afterIterator->sequenceNumber - beforeIterator->sequenceNumber);
        dimension_t segmentDistanceToVertex = (rawDistanceToVertex + parent.count()) % parent.count();
        dimension_t segmentDistanceToAfter = (rawDistanceToAfter + parent.count()) % parent.count();
        return segmentDistanceToVertex < segmentDistanceToAfter;
    } else if (beforeParent == vertexParent) {
        // before and vertex are in the same segment, while after is not, so we only need to compare before and vertex
        if (!vertexParent.reversed) {
            return vertexIterator->sequenceNumber >= beforeIterator->sequenceNumber;
        } else {
            return vertexIterator->sequenceNumber <= beforeIterator->sequenceNumber;
        }
    } else if (vertexParent == afterParent) {
        // vertex and after are in the same segment, while before is not, so we only need to compare vertex and after
        if (!vertexParent.reversed) {
            return vertexIterator->sequenceNumber < afterIterator->sequenceNumber;
        } else {
            return vertexIterator->sequenceNumber > afterIterator->sequenceNumber;
        }
    } else if (beforeParent == afterParent) {
        // before and after are in the same segment, while vertex is not, so we only need to compare before and after
        if (!beforeParent.reversed) {
            return afterIterator->sequenceNumber < beforeIterator->sequenceNumber;
        } else {
            return afterIterator->sequenceNumber > beforeIterator->sequenceNumber;
        }
    } else {
        // All three vertices are in different segments, so only the parent positions matter
        dimension_t parentDistanceToVertex =
                (vertexParent.sequenceNumber - beforeParent.sequenceNumber + parents.size()) % parents.size();
        dimension_t parentDistanceToAfter =
                (afterParent.sequenceNumber - beforeParent.sequenceNumber + parents.size()) % parents.size();
        return parentDistanceToVertex < parentDistanceToAfter;
    }
}

void TwoLevelTreeTour::mergeHalfSegment(vertex_t v, bool mergeToTheRight) {
    auto parentIterator = iterators[v]->parentIterator;
    SegmentParent &parent = *parentIterator;

    // Isolate the half segment we want to merge
    std::list<SegmentVertex> halfSegment;
    if (mergeToTheRight != parent.reversed) {
        halfSegment.splice(halfSegment.begin(), parent.vertices, iterators[v], parent.vertices.end());
    } else {
        halfSegment.splice(halfSegment.begin(), parent.vertices, parent.vertices.begin(), std::next(iterators[v]));
    }

    // Determine the other segment to merge with
    std::list<SegmentParent>::iterator otherParentIterator;
    if (mergeToTheRight) {
        otherParentIterator = std::next(parentIterator) != parents.end() ? std::next(parentIterator) : parents.begin();
    } else {
        otherParentIterator = parentIterator != parents.begin() ? std::prev(parentIterator) : std::prev(parents.end());
    }
    SegmentParent &otherParent = *otherParentIterator;

    // Reverse the half segment, if necessary, so that the direction of segment and halfSegment are the same
    if (parent.reversed != otherParent.reversed) {
        halfSegment.reverse();
    }

    // Determine where to insert the half-segment in the other segment and update sequenceNumber and parentIterator for
    // all elements in halfSegment
    std::list<SegmentVertex>::iterator insertIterator;
    long sequenceNumber;
    if (mergeToTheRight != otherParent.reversed) {
        insertIterator = otherParent.vertices.begin(); // Insert at the front of otherParent.vertices
        sequenceNumber = otherParent.vertices.front().sequenceNumber;
        for (auto it = halfSegment.rbegin(); it != halfSegment.rend(); ++it) {
            it->sequenceNumber = --sequenceNumber;
            it->parentIterator = otherParentIterator;
        }
    } else {
        insertIterator = otherParent.vertices.end(); // Insert at the back of otherParent.vertices
        sequenceNumber = otherParent.vertices.back().sequenceNumber;
        for (auto &it : halfSegment) {
            it.sequenceNumber = ++sequenceNumber;
            it.parentIterator = otherParentIterator;
        }
    }

    // Insert the halfSegment in otherParent.vertices
    otherParent.vertices.splice(insertIterator, halfSegment);
}

void TwoLevelTreeTour::reverseParents(std::list<SegmentParent>::iterator first,
                                      std::list<SegmentParent>::iterator last) {
    auto insertIterator = std::next(last);

    // Move the elements to a separate list, reverse this list and move them back to parents
    std::list<SegmentParent> temporaryList{};
    temporaryList.splice(temporaryList.begin(), parents, first, std::next(last));
    temporaryList.reverse();
    parents.splice(insertIterator, temporaryList);

    // Reverse the sequenceNumbers, so that they stay consecutive and swap the reversal bits for the affected parents
    std::swap(first, last);
    last++;
    while (first != last) {
        if (first == --last) {
            first->reversed = !first->reversed;
            break;
        }
        first->reversed = !first->reversed;
        last->reversed = !last->reversed;
        std::swap(first->sequenceNumber, last->sequenceNumber);
        first++;
    }
}


void TwoLevelTreeTour::flip(vertex_t a, vertex_t b, vertex_t c, vertex_t d) {
    // Get the segment vertices and parents corresponding to a, b, c and d
    SegmentVertex &aVertex = *iterators[a];
    SegmentVertex &bVertex = *iterators[b];
    SegmentVertex &cVertex = *iterators[c];
    SegmentVertex &dVertex = *iterators[d];

    SegmentParent &aParent = *(aVertex.parentIterator);
    SegmentParent &bParent = *(bVertex.parentIterator);
    SegmentParent &cParent = *(cVertex.parentIterator);
    SegmentParent &dParent = *(dVertex.parentIterator);

    // Case 1: The paths a-c and d-b are made up of segments
    // Reverse a range of consecutive parents that correspond to a-c or d-b
    if (aParent.firstVertex() == a and cParent.lastVertex() == c) {
        // Reverse one of these paths
        if (cParent.sequenceNumber >= aParent.sequenceNumber) {
            reverseParents(aVertex.parentIterator, cVertex.parentIterator);
        } else {
            reverseParents(dVertex.parentIterator, bVertex.parentIterator);
        }

        return;
    }

    // Case 2: One of the paths a-c and d-b is contained in a single segment
    // Reverse the paths inside the segment just as in ArrayTour when the path length is smaller than (3/4)*groupSize
    // otherwise split off and merge the ends of the segment with the adjacent segments, so that the path itself
    // becomes a segment and can be reversed by changing the reversal bit
    bool acInOneSegment = aParent == cParent and
                          ((!aParent.reversed and aVertex.sequenceNumber <= cVertex.sequenceNumber) or
                           (aParent.reversed and aVertex.sequenceNumber >= cVertex.sequenceNumber));
    bool dbInOneSegment = dParent == bParent and
                          ((!dParent.reversed and dVertex.sequenceNumber <= bVertex.sequenceNumber) or
                           (dParent.reversed and dVertex.sequenceNumber >= bVertex.sequenceNumber));
    if (acInOneSegment or dbInOneSegment) {
        // Determine which path will be reversed
        vertex_t start, end;
        if (acInOneSegment) {
            start = a;
            end = c;
        } else {
            start = d;
            end = b;
        }
        SegmentVertex &startVertex = *iterators[start];
        SegmentVertex &endVertex = *iterators[end];
        SegmentParent &parent = *(startVertex.parentIterator);

        dimension_t length = static_cast<dimension_t>(std::abs(endVertex.sequenceNumber - startVertex.sequenceNumber));
        if (length >= (groupSize * 3) / 4) {
            // This is a special case to reduce running time and achieve implicit rebalancing
            // Split off the two ends of the segments and merge them with their neighboring segments, so that the
            // complete segment only consists of the path that needs to be flipped. The actual flip (case 1) is
            // done with a new call to flip.
            if (start != parent.firstVertex()) {
                mergeHalfSegment(predecessor(start), false);
            }
            if (end != parent.lastVertex()) {
                mergeHalfSegment(successor(end), true);
            }
            flip(a, b, c, d); // This is now handled using case 1
        } else {
            // Flip the path inside of parent.vertices
            if (parent.reversed) {
                // iterators[start] needs to come before iterators[end] inside the segment
                std::swap(start, end);
            }
            parent.reverseVertices(iterators[start], iterators[end]);
        }

        return;
    }

    // Case 3: None of the above
    // Split up the segments so that case 1 or 2 is applicable
    if (aParent == bParent) {
        // Split the segment of aParent and merge one half into the neighboring segment
        long frontSequenceNumber = aParent.vertices.front().sequenceNumber;
        long backSequenceNumber = aParent.vertices.back().sequenceNumber;
        long aHalfSegmentLength = aParent.reversed ? aVertex.sequenceNumber - frontSequenceNumber :
                                  backSequenceNumber - aVertex.sequenceNumber;
        long bHalfSegmentLength = bParent.reversed ? backSequenceNumber - bVertex.sequenceNumber :
                                  bVertex.sequenceNumber - frontSequenceNumber;
        if (aHalfSegmentLength <= bHalfSegmentLength) {
            mergeHalfSegment(a, true);
        } else {
            mergeHalfSegment(b, false);
        }
    }

    if (cParent == dParent) {
        // Split the segment of cParent and merge one half into the neighboring segment
        long frontSequenceNumber = dParent.vertices.front().sequenceNumber;
        long backSequenceNumber = dParent.vertices.back().sequenceNumber;
        long dHalfSegmentLength = dParent.reversed ? dVertex.sequenceNumber - frontSequenceNumber :
                                  backSequenceNumber - dVertex.sequenceNumber;
        long cHalfSegmentLength = cParent.reversed ? backSequenceNumber - dVertex.sequenceNumber :
                                  cVertex.sequenceNumber - frontSequenceNumber;
        if (dHalfSegmentLength <= cHalfSegmentLength) {
            mergeHalfSegment(d, true);
        } else {
            mergeHalfSegment(c, false);
        }
    }

    // Normally, a-b and c-d now are the boundaries of segments, so it can be handled by case 1
    // If cParent (= dParent) was the neighboring segment of aParent (= bParent) it can happen that a half of aParent
    // gets merged to cParent and then the half that contains the newly merged vertices gets merged back to aParent. In
    // this case, the flip can now be handled by case 2.
    flip(a, b, c, d);
}
