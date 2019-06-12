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

// TODO: Replace swap with std::reverse

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
    dimension_t distanceAC = distance(a, c);
    dimension_t distanceDB = distance(d, b);
    dimension_t segmentStartIndex, segmentEndIndex;
    dimension_t segmentLength;
    if (distanceAC <= distanceDB) {
        // The segment a-c will be reversed
        segmentStartIndex = indices[a];
        segmentEndIndex = indices[c];
        segmentLength = distanceAC;
    } else {
        // The segment d-b will be reversed
        segmentStartIndex = indices[d];
        segmentEndIndex = indices[b];
        segmentLength = distanceDB;
    }
    for (dimension_t i = 0; i < (segmentLength + 1) / 2; ++i) {
        vertex_t &vertex1 = sequence[(segmentStartIndex + i) % getDimension()];
        vertex_t &vertex2 = sequence[(segmentEndIndex + getDimension() - i) % getDimension()];
        std::swap(vertex1, vertex2);
        std::swap(indices[vertex1], indices[vertex2]);
    }
}


// ============================================ TwoLevelTreeTour class =================================================

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

const TwoLevelTreeTour::SegmentParent &TwoLevelTreeTour::getPreviousParent(
        std::list<SegmentParent>::iterator parentIterator) const {
    return parentIterator != parents.begin() ? *std::prev(parentIterator) : parents.back();
}

const TwoLevelTreeTour::SegmentParent &TwoLevelTreeTour::getNextParent(
        std::list<SegmentParent>::iterator parentIterator) const {
    return std::next(parentIterator) != parents.end() ? *std::next(parentIterator) : parents.front();
}

void TwoLevelTreeTour::reverse(std::list<SegmentVertex> &list) {
    list.reverse();
    auto first = list.begin();
    auto last = list.end();
    while ((first != last) && (first != --last)) {
        std::swap((first++)->sequenceNumber, last->sequenceNumber);
    }
}

void TwoLevelTreeTour::reverse(std::list<SegmentVertex> &list, std::list<SegmentVertex>::iterator first,
                               std::list<SegmentVertex>::iterator last) {
    auto insertIterator = std::next(last);
    std::list<SegmentVertex> temporaryList;
    temporaryList.splice(temporaryList.begin(), list, first, last);
    reverse(temporaryList);
    list.splice(insertIterator, temporaryList);

}

void TwoLevelTreeTour::reverse(std::list<SegmentParent> &list) {
    list.reverse();
    auto first = list.begin();
    auto last = list.end();
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

void TwoLevelTreeTour::reverse(std::list<SegmentParent> &list, std::list<SegmentParent>::iterator first,
                               std::list<SegmentParent>::iterator last) {
    auto insertIterator = std::next(last);
    std::list<SegmentParent> temporaryList;
    temporaryList.splice(temporaryList.begin(), list, first, last);
    reverse(temporaryList);
    list.splice(insertIterator, temporaryList);
}


void TwoLevelTreeTour::setVertices(const std::vector<vertex_t> &vertexList) {
    dimension = vertexList.size();
    if (vertexList.size() <= 100000) {
        groupSize = 100;
    } else {
        groupSize = 200;
    }

    iterators.assign(vertexList.size(), std::list<SegmentVertex>::iterator());

    dimension_t segmentLength = groupSize;
    size_t parentIndex = 0;
    size_t vertexIndex = 0;
    while (vertexIndex < dimension) {
        // The last segment should have a length between groupSize and 2*groupSize
        if (vertexList.size() - vertexIndex < 2 * groupSize) {
            segmentLength = vertexList.size() - vertexIndex;
        }

        SegmentParent parent{};
        parent.reversed = false;
        parent.sequenceNumber = parentIndex;
        parents.push_back(parent);
        auto parentIterator = std::prev(parents.end());
        for (; vertexIndex < vertexIndex + segmentLength; vertexIndex++) {
            parentIterator->vertices.push_back(
                    SegmentVertex{vertexList[vertexIndex], parentIterator, static_cast<long>(vertexIndex)});
            iterators[vertexList[vertexIndex]] = std::prev(parent.vertices.end());
        }

        parentIndex++;
    }
}

TwoLevelTreeTour::TwoLevelTreeTour(const std::vector<vertex_t> &vertexList) {
    setVertices(vertexList);
}

dimension_t TwoLevelTreeTour::getDimension() const {
    return dimension;
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
    dimension_t parentDistanceToVertex =
            (vertexParent.sequenceNumber - beforeParent.sequenceNumber + parents.size()) % parents.size();
    dimension_t parentDistanceToAfter =
            (afterParent.sequenceNumber - beforeParent.sequenceNumber + parents.size()) % parents.size();
    if (parentDistanceToVertex < parentDistanceToAfter) {
        return true;
    } else {
        // vertex and after have the same parent
        const SegmentParent &parent = vertexParent;
        if (beforeParent == vertexParent) {
            // All three vertices have the same parent and are therefore in the same segment
            dimension_t segmentDistanceToVertex =
                    (vertexIterator->sequenceNumber - beforeIterator->sequenceNumber + parent.count()) % parent.count();
            dimension_t segmentDistanceToAfter =
                    (afterIterator->sequenceNumber - beforeIterator->sequenceNumber + parent.count()) % parent.count();
            if (!parent.reversed) {
                return segmentDistanceToVertex < segmentDistanceToAfter;
            } else {
                return segmentDistanceToVertex > segmentDistanceToAfter;
            }
        } else {
            // Since before is in a different segment, only the order inside the segment matters
            if (!parent.reversed) {
                return vertexIterator->sequenceNumber < afterIterator->sequenceNumber;
            } else {
                return vertexIterator->sequenceNumber > afterIterator->sequenceNumber;
            }
        }
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
        halfSegment.splice(halfSegment.begin(), parent.vertices, parent.vertices.begin(), iterators[v]);
    }

    // Determine the other segment to merge with
    std::list<SegmentParent>::iterator otherParentIterator;
    if (mergeToTheRight) {
        otherParentIterator = std::next(parentIterator) != parents.end() ? std::next(parentIterator) : parents.begin();
    } else {
        otherParentIterator = parentIterator != parents.begin() ? std::prev(parentIterator) : parents.end();
    }
    SegmentParent &otherParent = *otherParentIterator;

    // Reverse the half segment, if necessary, so that the direction of segment and halfSegment are the same
    if (parent.reversed != otherParent.reversed) {
        halfSegment.reverse();
    }

    // Determine where to insert the half-segment in the other segment and update sequence number and parent for all
    // elements in halfSegment
    std::list<SegmentVertex>::iterator insertIterator;
    std::list<SegmentVertex>::iterator countStartIterator;
    long sequenceNumber;
    if (mergeToTheRight != otherParent.reversed) {
        insertIterator = otherParent.vertices.begin();
        sequenceNumber = otherParent.vertices.front().sequenceNumber;
        for (auto it = halfSegment.rbegin(); it != halfSegment.rend(); ++it) {
            it->sequenceNumber = --sequenceNumber;
            it->parentIterator = otherParentIterator;
        }
    } else {
        insertIterator = otherParent.vertices.end();
        sequenceNumber = otherParent.vertices.back().sequenceNumber;
        for (auto it = halfSegment.begin(); it != halfSegment.end(); ++it) {
            it->sequenceNumber = ++sequenceNumber;
            it->parentIterator = otherParentIterator;
        }
    }

    // Merge
    otherParent.vertices.splice(insertIterator, halfSegment);
}


void TwoLevelTreeTour::flip(vertex_t a, vertex_t b, vertex_t c, vertex_t d) {
    SegmentVertex &aVertex = *iterators[a];
    SegmentVertex &bVertex = *iterators[b];
    SegmentVertex &cVertex = *iterators[c];
    SegmentVertex &dVertex = *iterators[d];

    SegmentParent &aParent = *(aVertex.parentIterator);
    SegmentParent &bParent = *(bVertex.parentIterator);
    SegmentParent &cParent = *(cVertex.parentIterator);
    SegmentParent &dParent = *(dVertex.parentIterator);

    // Case 1: The paths a-c and d-b are made up of segments
    if (aParent.firstVertex() == a and cParent.lastVertex() == c) {
        std::cout << "flip Case 1: " << a << ", " << b << ", " << c << ", " << d << std::endl;
        dimension_t numberOfParentsAC =
                (cParent.sequenceNumber - aParent.sequenceNumber + parents.size()) % parents.size();
        dimension_t numberOfParentsDB =
                (bParent.sequenceNumber - dParent.sequenceNumber + parents.size()) % parents.size();

        std::list<SegmentParent>::iterator parentStartIterator, parentEndIterator;
        if (numberOfParentsAC <= numberOfParentsDB) {
            parentStartIterator = aVertex.parentIterator;
            parentEndIterator = cVertex.parentIterator;
        } else {
            parentStartIterator = dVertex.parentIterator;
            parentEndIterator = bVertex.parentIterator;
        }

        reverse(parents, parentStartIterator, parentEndIterator);

        return;
    }

    // Case 2: The path a-c or d-b is contained in a single segment
    bool acInOneSegment = aParent == cParent and aVertex.sequenceNumber <= cVertex.sequenceNumber;
    bool dbInOneSegment = dParent == bParent and dVertex.sequenceNumber <= bVertex.sequenceNumber;
    if (acInOneSegment or dbInOneSegment) {
        std::cout << "flip Case 2: " << a << ", " << b << ", " << c << ", " << d << std::endl;
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
        dimension_t length = endVertex.sequenceNumber - startVertex.sequenceNumber;
        if (length >= (parent.count() * 3) / 4) {
            if (start != parent.firstVertex()) {
                mergeHalfSegment(predecessor(start), false);
            }
            if (end != parent.lastVertex()) {
                mergeHalfSegment(successor(end), true);
            }
            flip(a, b, c, d);
        } else {
            reverse(parent.vertices, iterators[start], iterators[end]);
        }
        return;
    }

    // Case 3:
    // Split up the segments so that case 1 is applicable
    std::cout << "flip Case 3: " << a << ", " << b << ", " << c << ", " << d << std::endl;
    if (aParent == bParent) {
        // Split the segment of aParent and merge it into the neighboring segment
        dimension_t aHalfSegmentLength = aParent.vertices.back().sequenceNumber - aVertex.sequenceNumber;
        dimension_t bHalfSegmentLength = bVertex.sequenceNumber - bParent.vertices.front().sequenceNumber;
        if (aHalfSegmentLength <= bHalfSegmentLength) {
            mergeHalfSegment(a, true);
        } else {
            mergeHalfSegment(b, false);
        }
    }

    if (cParent == dParent) {
        // Split the segment of cParent and merge it into the neighboring segment
        dimension_t dHalfSegmentLength = dParent.vertices.back().sequenceNumber - dVertex.sequenceNumber;
        dimension_t cHalfSegmentLength = cVertex.sequenceNumber - cParent.vertices.front().sequenceNumber;
        if (dHalfSegmentLength <= cHalfSegmentLength) {
            mergeHalfSegment(d, true);
        } else {
            mergeHalfSegment(c, false);
        }
    }
    flip(a, b, c, d);
}


// =============================================== TourWalker class ====================================================


// TODO: Use std::ostream_iterator to make std::ostream &operator<< more beautiful
/*
template <typename T>
void print(std::list<T> & listObj)
{
    std::copy(listObj.begin(), listObj.end(), std::ostream_iterator<T>(std::cout, " , "));
    std::cout<<std::endl;
}
*/


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

