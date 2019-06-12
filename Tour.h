//
// Created by Karl Welzel on 19/04/2019.
//

#ifndef LINKERNINGHANALGORITHM_TOUR_H
#define LINKERNINGHANALGORITHM_TOUR_H

#include <utility>
#include <limits>
#include <list>
#include <vector>
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
    // This function expects a vector containing the numbers from 0 to tour.size()-1 and clears existing neighbors
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


// ================================================ ArrayTour class ====================================================

// This is the simplest implementation of a Tour using two vectors (as replacement for arrays)

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

using Tour = ArrayTour;


// ============================================ TwoLevelTreeTour class =================================================


class TwoLevelTreeTour : public BaseTour {
private:
    struct SegmentParent {
        std::vector<vertex_t> segmentVertices;
        bool reversed;
        dimension_t index;

        dimension_t count() {
            return segmentVertices.size();
        }
    };

    dimension_t dimension = 0;
    dimension_t groupSize = 0;
    std::vector<SegmentParent> parents;
    std::vector<std::pair<dimension_t, dimension_t>> indices;

public:
    TwoLevelTreeTour() = default;

    // Initialize the tour with a sequence of vertices
    // This function expects a vector containing the numbers from 0 to tour.size()-1 and overrides all data
    void setVertices(const std::vector<vertex_t> &vertexList) override;

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


// =============================================== TourWalker class ====================================================

// This class provides a way to walk through a Tour. It stores a reference to the tour, its current and next vertex
// (to store the direction) and can be advanced by "getNextVertex".


class TourWalker {
private:
    const BaseTour &tour;
    vertex_t current;

public:
    // Create a TourWalker that starts the walk at vertex first
    TourWalker(const BaseTour &tour, vertex_t first);

    // Advance the walk and get the next vertex
    vertex_t getNextVertex();
};

// Output the tour to the stream out, typically used with std::cout
std::ostream &operator<<(std::ostream &out, const BaseTour &tour);

#endif //LINKERNINGHANALGORITHM_TOUR_H
