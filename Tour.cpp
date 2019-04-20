//
// Created by Karl Welzel on 19/04/2019.
//

#include <utility>
#include "Tour.h"


// ================================================ VertexList class ===================================================

// This class maintains a map that stores a maximum of two neighbors of each vertex. There is no inherent direction
// encoded in these neighbors. This gives additional flexibility when compared to std::list while still having similar
// run times. It is possible to create tours that do not contain every vertex or multiple separate paths or even have a
// vertex x with neighbor y while y does not have x as its neighbor, so extra caution is necessary when manipulating
// these objects.
// This class is the base class for Tour.

const vertex_t VertexList::NO_VERTEX = std::numeric_limits<unsigned int>::max();

VertexList::VertexList() = default;

// Create dimension vertices, each with no neighbors
VertexList::VertexList(const dimension_t dimension) {
    neighbors.assign(dimension, std::make_pair(NO_VERTEX, NO_VERTEX));
}

VertexList::VertexList(std::vector<std::pair<vertex_t, vertex_t>> neighbors) :
        neighbors(std::move(neighbors)) {}

dimension_t VertexList::getDimension() const {
    return neighbors.size();
}

// Which vertex comes after current, when previous comes before it?
// previous is necessary to determine the direction
vertex_t VertexList::next(const vertex_t previous, const vertex_t current) const {
    const std::pair<vertex_t, vertex_t> currentNeighbors = neighbors.at(current);
    if (previous != currentNeighbors.first and previous != currentNeighbors.second) {
        throw std::runtime_error(
                "Tour::next: previous (" + std::to_string(previous) + ") is not a neighbor of current (" +
                std::to_string(current) + ")");
    } else if (previous == currentNeighbors.first) {
        return currentNeighbors.second;
    } else { // previous == currentNeighbors.second
        return currentNeighbors.first;
    }
}

// Returns some vertex, that is a neighbor of current, if one exists, otherwise NO_VERTEX
vertex_t VertexList::next(const vertex_t current) const {
    std::pair<vertex_t, vertex_t> currentNeighbors = neighbors.at(current);
    if (currentNeighbors.first != NO_VERTEX) {
        return currentNeighbors.first;
    } else if (currentNeighbors.second != NO_VERTEX) {
        return currentNeighbors.second;
    } else {
        return NO_VERTEX;
    }
}

// Set the vertex that comes after current, when previous comes before it. The neighbors of next are not changed!
// previous is necessary to determine the direction
void VertexList::setNext(const vertex_t previous, const vertex_t current, const vertex_t next) {
    if (previous == next) {
        throw std::runtime_error("Tour::setNext: You cannot set both neighbors to the same vertex.");
    }
    const std::pair<vertex_t, vertex_t> currentNeighbors = neighbors.at(current);
    if (previous != currentNeighbors.first and previous != currentNeighbors.second) {
        throw std::runtime_error(
                "Tour::setNext: previous (" + std::to_string(previous) + ") is not a neighbor of current (" +
                std::to_string(current) + ")");
    }
    neighbors.at(current) = std::make_pair(previous, next);
}

// Tries to make vertex1 a neighbor of vertex2 and vertex2 a neighbor of vertex1 and returns whether this was
// successful. It only fails if one of them already has two neighbors.
bool VertexList::makeNeighbors(const vertex_t vertex1, const vertex_t vertex2) {
    const std::pair<vertex_t, vertex_t> neighbors1 = neighbors.at(vertex1);
    const std::pair<vertex_t, vertex_t> neighbors2 = neighbors.at(vertex2);
    if (neighbors1.first == NO_VERTEX and neighbors2.first == NO_VERTEX) {
        neighbors.at(vertex1).first = vertex2;
        neighbors.at(vertex2).first = vertex1;
        return true;
    } else if (neighbors1.first == NO_VERTEX and neighbors2.second == NO_VERTEX) {
        neighbors.at(vertex1).first = vertex2;
        neighbors.at(vertex2).second = vertex1;
        return true;
    } else if (neighbors1.second == NO_VERTEX and neighbors2.first == NO_VERTEX) {
        neighbors.at(vertex1).second = vertex2;
        neighbors.at(vertex2).first = vertex1;
        return true;
    } else if (neighbors1.second == NO_VERTEX and neighbors2.second == NO_VERTEX) {
        neighbors.at(vertex1).second = vertex2;
        neighbors.at(vertex2).second = vertex1;
        return true;
    } else {
        return false;
    }
}


// =================================================== Tour class ======================================================

// This class represents a Tour. It is a subclass of VertexList, but avoids the problems of VertexList by only allowing
// changes that don't destroy the tour property


Tour::Tour() = default;

Tour::Tour(std::vector<std::pair<vertex_t, vertex_t>> neighbors) : VertexList(std::move(neighbors)) {
    if (!isHamiltonianTour()) {
        throw std::runtime_error(
                "A Tour object cannot be initialized with neighbors that do not describe a hamiltonian tour");
    }
}

// Initialize the neighbors map with a std::list that represents a tour. For each vertex x the neighbors are the
// adjacent entries in "vertexList" (start and end are also adjacent)
// This expects a list containing the numbers from 0 to tour.size()-1 and clears neighbors
void Tour::setVertices(const std::list<vertex_t> &vertexList) {
    neighbors.assign(vertexList.size(), std::make_pair(NO_VERTEX, NO_VERTEX));
    vertex_t previous, current;
    auto it = vertexList.begin();
    current = *it;
    std::advance(it, 1);
    for (; it != vertexList.end(); ++it) {
        previous = current;
        current = *it;
        makeNeighbors(previous, current);
    }
    makeNeighbors(vertexList.back(), vertexList.front());
}

// Checks if this Tour really is a hamiltonian tour
bool Tour::isHamiltonianTour() const {
    dimension_t dimension = neighbors.size();
    std::vector<bool> visited(dimension, false);
    TourWalker tourWalker(*this, 0);
    vertex_t currentVertex = tourWalker.getNextVertex();
    do {
        if (currentVertex >= visited.size() or visited.at(currentVertex)) {
            return false;
        }
        visited.at(currentVertex) = true;
        currentVertex = tourWalker.getNextVertex();
    } while (currentVertex != 0);
    for (bool v : visited) {
        if (!v) return false;
    }
    return true;
}


// =============================================== TourWalker class ====================================================

// This class provides a way to walk through a Tour. It stores a copy of the tour, its current and next vertex (to store
// the direction) and can be advanced by "getNextVertex".


// Create a TourWalker that starts the walk at vertex first
TourWalker::TourWalker(const Tour &tour, vertex_t first) : TourWalker(tour, first, tour.next(first)) {}

// Create a TourWalker that starts the walk at vertex first and then walks in the direction of vertex second
TourWalker::TourWalker(Tour tour, vertex_t first, vertex_t second) : tour(std::move(tour)),
                                                                     current(first), next(second) {}

// Advance the walk and get the next vertex
vertex_t TourWalker::getNextVertex() {
    const vertex_t previous = current;
    current = next;
    next = tour.next(previous, current);
    return previous;
}

// Output the tour to the stream out, typically used with std::cout
std::ostream &operator<<(std::ostream &out, const Tour &tour) {
    std::string output;
    TourWalker tourWalker(tour, 0);
    vertex_t currentVertex = tourWalker.getNextVertex();
    do {
        output += std::to_string(currentVertex + 1) + ", ";
        currentVertex = tourWalker.getNextVertex();
    } while (currentVertex != 0);
    out << output.substr(0, output.length() - 2) << std::endl;
    return out;
}