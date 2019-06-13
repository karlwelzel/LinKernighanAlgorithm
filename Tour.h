//
// Created by Karl Welzel on 19/04/2019.
//

#ifndef LINKERNINGHANALGORITHM_TOUR_H
#define LINKERNINGHANALGORITHM_TOUR_H

#include <utility>
#include <limits>
#include <list>
#include <vector>
#include <list>
#include <iostream>

using vertex_t = size_t; // The type used for vertices
using dimension_t = vertex_t; // The type used for counting vertices
using distance_t = unsigned long; // The type used for distances and lengths
using signed_distance_t = long long; // The type used for distances and lengths that can be negative
// using Tour = ArrayTour; // The current implementation of a tour, that should be used


// ================================================= BaseTour class ====================================================

// This class represents a Tour. It is an abstract base class for different implementations of a tour.

class BaseTour {
public: // pure virtual functions that need to be implemented by subclasses

    // Returns the number of vertices in the tour
    virtual dimension_t getDimension() const = 0;

    // Initialize the tour with a sequence of vertices
    // This function expects a vector containing the numbers from 0 to tour.size()-1 and overrides all data
    virtual void setVertices(const std::vector<vertex_t> &vertexList) = 0;

    // Returns the predecessor of vertex in the current tour
    virtual vertex_t predecessor(vertex_t vertex) const = 0;

    // Returns the successor of vertex in the current tour
    virtual vertex_t successor(vertex_t vertex) const = 0;

    // Checks whether "vertex" is reached before "after" when walking in forward direction from "before"
    virtual bool isBetween(vertex_t before, vertex_t vertex, vertex_t after) const = 0;

    // Performs a 2-opt exchange: Replaces {a, b} and {c, d} by {b, c} and {d, a}
    // The new tour is a-c, b-d
    // Expects successor(b) = a and successor(c) = d
    virtual void flip(vertex_t a, vertex_t b, vertex_t c, vertex_t d) = 0;

public: // functions that all Tour classes have in common and that only depend on the functions above

    // Compute the inverse permutation to a permutation of the numbers 0 to n-1, i.e. a vector inv such that
    // permutation[inv[i]] == i for all 0 <= i < n
    // Expects that permutation contains every number 0 to n-1 exactly once
    static std::vector<dimension_t> inversePermutation(const std::vector<dimension_t> &permutation);

    // Returns the neighbors of vertex in the tour
    std::vector<vertex_t> getNeighbors(vertex_t vertex);

    // Checks whether the tour contains the edge {vertex1, vertex2}
    bool containsEdge(vertex_t vertex1, vertex_t vertex2) const;

    // Computes the cyclic permutation of the vertices in an alternating walk, i.e. a vector permutation such that
    // the elements in alternatingWalk appear on the tour in the order
    //   alternatingWalk[permutation[0]], alternatingWalk[permutation[1]], ..., alternatingWalk[permutation[n-1]]
    // Expects a closed alternating walk that starts with an edge on the tour
    std::vector<dimension_t> cyclicPermutation(const std::vector<vertex_t> &alternatingWalk) const;

    // Checks if the tour after exchanging all edges of alternatingWalk on the tour by edges not on the tour is still
    // a hamiltonian tour, does not change the tour itself
    // Expects a closed alternating walk that starts with an edge on the tour
    bool isTourAfterExchange(const std::vector<vertex_t> &alternatingWalk) const;

    // Exchanges all edges of alternatingWalk on the tour by edges not on the tour
    // Expects a closed alternating walk that starts with an edge on the tour
    // Expects that the exchange will lead to a hamiltonian tour, check with isTourAfterExchange beforehand
    void exchange(const std::vector<vertex_t> &alternatingWalk);
};

// Output the tour to the stream out, typically used with std::cout
std::ostream &operator<<(std::ostream &out, const BaseTour &tour);


// ================================================ ArrayTour class ====================================================

// This is the simplest implementation of a Tour using a vector to store the tour order and a second one to store the
// index of each vertex in the former vector. The running times for predecessor, successor and isBetween is O(1) while
// flip has running time O(n)

class ArrayTour : public BaseTour {
private:
    std::vector<vertex_t> sequence; // Stores the sequence of vertices in the tour
    std::vector<dimension_t> indices; // Stores the index of each vertex in sequence: sequence[indices[i]] = i

    // Compute the number of vertices between vertex1 and vertex2 in the successor direction (including vertex1 and
    // excluding vertex2)
    dimension_t distance(vertex_t vertex1, vertex_t vertex2) const;

public:
    ArrayTour() = default;

    // Initialize the tour with a sequence of vertices
    // This function expects a vector containing the numbers from 0 to tour.size()-1 and overrides all data
    void setVertices(const std::vector<vertex_t> &vertexList) override;

    // Initialize the tour with a sequence of vertices
    // This function expects a vector containing the numbers from 0 to tour.size()-1 and overrides all data
    explicit ArrayTour(const std::vector<vertex_t> &vertexList);

    // Returns the number of vertices in the tour
    dimension_t getDimension() const override;

    // Returns the predecessor of vertex in the current tour
    vertex_t predecessor(vertex_t vertex) const override;

    // Returns the successor of vertex in the current tour
    vertex_t successor(vertex_t vertex) const override;

    // Checks whether "vertex" is reached before "after" when walking in forward direction from "before"
    bool isBetween(vertex_t before, vertex_t vertex, vertex_t after) const override;

    // Performs a 2-opt exchange: Replaces {a, b} and {c, d} by {b, c} and {d, a}
    // Expects successor(b) = a and successor(c) = d
    void flip(vertex_t a, vertex_t b, vertex_t c, vertex_t d) override;
};


// ============================================ TwoLevelTreeTour class =================================================

// This tour implementation uses a tree with two levels to store the vertices. The first level is a list of parent
// nodes (here: SegmentParent) and the second level is a segment of the tour for every parent node. Each vertex in a
// segment (here: SegmentVertex) has a pointer to its parent and there is a map that stores a pointer to the
// SegmentVertex inside one of the segments corresponding to a vertex (of type vertex_t). The advantage over the
// ArrayTour implementation is that each parent node stores a reversal bit, which indicates if the complete segment
// should be traversed in the opposite direction. So, when the complete segment needs to be reversed, it takes only
// O(1) time instead of a linear running time in the segment length.

// The implementation was adapted from the paper
// Fredman, M.L., Johnson, D.S., McGeoch, L.A., & Ostheimer, G. (1993). Data Structures for Traveling Salesmen. SODA.
// which is available here: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.49.570&rep=rep1&type=pdf
// It sacrifices the tight O(sqrt(n)) worst-case bound for a flip operation for a simpler and faster implementation.

// groupSize is a parameter for how long a single segment should be. The choice in setVertices is similar to the one
// suggested in the paper.

class TwoLevelTreeTour : public BaseTour {
private:
    struct SegmentParent;

    struct SegmentVertex {
        // The represented vertex
        vertex_t vertex;

        // A pointer to its parent
        std::list<SegmentParent>::iterator parentIterator;

        // Used to order the elements in the segment (in isBetween). Inside a segment the numbers are consecutive,
        // but they might start and end with arbitrary numbers.
        long sequenceNumber;
    };

    struct SegmentParent {
        // A list of the vertices in this segment
        std::list<SegmentVertex> vertices;

        // The reversal bit, indicates whether this segments orientation is opposite to the tour
        bool reversed = false;

        // Used to order the parent nodes (in isBetween). The numbers are consecutive between 0 and the number of
        // segments
        dimension_t sequenceNumber = 0;

        // Counts the number of elements in this parent (= vertices.size())
        dimension_t count() const;

        // Returns the first vertex in this segment with respect to the tour direction
        vertex_t firstVertex() const;

        // Returns the last vertex in this segment with respect to the tour direction
        vertex_t lastVertex() const;

        // Reverses the segment vertices in vertices in the range [first, last] while correctly changing the sequence
        // numbers
        // Expects that first and last are iterators of vertices
        void reverseVertices(std::list<SegmentVertex>::iterator first, std::list<SegmentVertex>::iterator last);
    };

    // Stores the number of vertices
    dimension_t dimension = 0;

    // The average size the segments should have. The parameter is described in more detail in the paper above
    dimension_t groupSize = 0;

    // A list of the parents in the order their segments appear on the tour
    std::list<SegmentParent> parents;

    // A map of each vertex to the pointer that points to its SegmentVertex in one of the parents vertices list
    // Invariant: iterators[v]->vertex == v. These iterators do not have to be updated, because they are list iterators
    // and will not be invalidated even when the elements are moved between different lists.
    std::vector<std::list<SegmentVertex>::iterator> iterators;

    // Compare two SegmentParents by equality
    friend bool operator==(const SegmentParent &parent, const SegmentParent &otherParent);

    // Returns the SegmentParent that comes before *parentIterator
    const SegmentParent &getPreviousParent(std::list<SegmentParent>::iterator parentIterator) const;

    // Returns the SegmentParent that comes after *parentIterator
    const SegmentParent &getNextParent(std::list<SegmentParent>::iterator parentIterator) const;

    // Reverses the SegmentParent in parents in the range [first, last] while correctly changing the sequence numbers
    // and updating the reversal bits
    // Expects that first and last are iterators of parents
    void reverseParents(std::list<SegmentParent>::iterator first, std::list<SegmentParent>::iterator last);

    // Merge the half-segment to the right (or left) of vertex v including v with the right (or left) neighbor segment.
    // The direction is given with respect to the tour direction (right = successor-direction,
    // left = predecessor-direction)
    void mergeHalfSegment(vertex_t v, bool mergeToTheRight);

public:
    TwoLevelTreeTour() = default;

    // The copy-constructor: Copy the state of otherTour into this tour
    TwoLevelTreeTour(const TwoLevelTreeTour &otherTour);

    // Initialize the tour with a sequence of vertices
    // This function expects a vector containing the numbers from 0 to tour.size()-1 and overrides all data
    void setVertices(const std::vector<vertex_t> &vertexList) override;

    // Initialize the tour with a sequence of vertices
    // This function expects a vector containing the numbers from 0 to tour.size()-1 and overrides all data
    explicit TwoLevelTreeTour(const std::vector<vertex_t> &vertexList);

    // Returns the number of vertices in the tour
    dimension_t getDimension() const override;

    // Returns the predecessor of vertex in the current tour
    vertex_t predecessor(vertex_t vertex) const override;

    // Returns the successor of vertex in the current tour
    vertex_t successor(vertex_t vertex) const override;

    // Checks whether "vertex" is reached before "after" when walking in forward direction from "before"
    bool isBetween(vertex_t before, vertex_t vertex, vertex_t after) const override;

    // Performs a 2-opt exchange: Replaces {a, b} and {c, d} by {b, c} and {d, a}
    // Expects successor(b) = a and successor(c) = d
    void flip(vertex_t a, vertex_t b, vertex_t c, vertex_t d) override;
};


using Tour = TwoLevelTreeTour;

#endif //LINKERNINGHANALGORITHM_TOUR_H
