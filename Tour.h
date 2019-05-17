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


// ================================================ VertexList class ===================================================

// This class maintains a map that stores a maximum of two neighbors of each vertex. There is no inherent direction
// encoded in these neighbors. This gives additional flexibility when compared to std::list while still having similar
// run times. It is possible to create tours that do not contain every vertex or multiple separate paths or even have a
// vertex x with neighbor y while y does not have x as its neighbor, so extra caution is necessary when manipulating
// these objects.
// This class is the base class for Tour.


class VertexList {
protected:
    std::vector<std::vector<vertex_t>> neighbors;

    // Set newNeighbor as one of the neighbors of vertex.
    // Throws a runtime_error if both of the neighbors of vertex are already set
    void addNeighbor(vertex_t vertex, vertex_t newNeighbor);

    // Deletes neighbor as a neighbor of vertex
    // Throws a runtime_error if neighbor is not a neighbor of vertex
    void removeNeighbor(vertex_t vertex, vertex_t neighbor);

public:
    VertexList();

    // Create dimension vertices, each with no neighbors
    explicit VertexList(dimension_t dimension);

    // Create VertexList from a given neighbors map
    explicit VertexList(std::vector<std::vector<vertex_t>> neighbors);

    // Returns the number of vertices
    dimension_t getDimension() const;

    // Returns the neighbors of vertex in the tour
    std::vector<vertex_t> getNeighbors(vertex_t vertex);

    // Which vertex comes after current, when previous comes before it?
    // previous is necessary to determine the direction
    // Throws a runtime_error if the vertex cannot be determined
    vertex_t next(vertex_t previous, vertex_t current) const;

    // Returns some vertex, that is a neighbor of current
    // Throws a runtime_error if current has no neighbor
    vertex_t next(vertex_t current) const;

    // Checks whether next is a neighbor of current
    bool isNeighbor(vertex_t vertex, vertex_t neighbor) const;

    // Checks whether the tour contains the edge {vertex1, vertex2}
    // This function is an alias for isNeighbor
    bool containsEdge(vertex_t vertex1, vertex_t vertex2) const;

    // Tries to make vertex1 a neighbor of vertex2 and vertex2 a neighbor of vertex1 and returns whether this was
    // successful. It only fails if one of them already has two neighbors.
    bool makeNeighbors(vertex_t vertex1, vertex_t vertex2);

};


// =================================================== Tour class ======================================================

// This class represents a Tour. It is a subclass of VertexList, but avoids the problems of VertexList by only allowing
// changes that don't destroy the tour property


class Tour : public VertexList {
protected:
    // Initialize the neighbors map with a std::list that represents a tour. For each vertex x the neighbors are the
    // adjacent entries in "vertexList" (start and end are also adjacent)
    // This function expects a list containing the numbers from 0 to tour.size()-1 and clears existing neighbors
    void setVertices(const std::vector<vertex_t> &vertexList);

public:
    Tour();

    // Create VertexList from a given neighbors map, checks if the map represents a tour
    explicit Tour(std::vector<std::vector<vertex_t>> neighbors);

    // see setVertices
    explicit Tour(const std::vector<vertex_t> &vertexList);

    // Checks if this Tour really is a hamiltonian tour
    bool isHamiltonianTour() const;

    // Checks if the tour after exchanging all edges of alternatingWalk on the tour by edges not on the tour is still
    // a hamiltonian tour
    // Expects a closed alternating walk that starts with an edge on the tour
    bool isTourAfterExchange(const std::vector<vertex_t> &alternatingWalk) const;

    // Computes the tour after exchanging all edges of alternatingWalk on the tour by edges not on the tour, but does
    // not modify the tour object itself
    // Expects a closed alternating walk that starts with an edge on the tour
    // Expects that the exchange will lead to a hamiltonian tour, check with isTourAfterExchange beforehand
    Tour exchange(const std::vector<vertex_t> &alternatingWalk) const;
};


// =============================================== TourWalker class ====================================================

// This class provides a way to walk through a Tour. It stores a copy of the tour, its current and next vertex (to store
// the direction) and can be advanced by "getNextVertex".


class TourWalker {
private:
    const Tour tour;
    vertex_t current;
    vertex_t next;

public:
    // Create a TourWalker that starts the walk at vertex first
    TourWalker(const Tour &tour, vertex_t first);

    // Create a TourWalker that starts the walk at vertex first and then walks in the direction of vertex second
    TourWalker(Tour tour, vertex_t first, vertex_t second);

    // Advance the walk and get the next vertex
    vertex_t getNextVertex();
};

// Output the tour to the stream out, typically used with std::cout
std::ostream &operator<<(std::ostream &out, const Tour &tour);

#endif //LINKERNINGHANALGORITHM_TOUR_H
