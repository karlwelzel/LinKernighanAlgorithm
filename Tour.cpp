//
// Created by Karl Welzel on 19/04/2019.
//

#include <utility>
#include <limits>
#include <list>
#include <vector>
#include <iostream>
#include "Tour.h"


// ================================================ VertexList class ===================================================

const vertex_t VertexList::NO_VERTEX = std::numeric_limits<vertex_t>::max();

VertexList::VertexList() = default;

VertexList::VertexList(const dimension_t dimension) {
    neighbors.assign(dimension, std::make_pair(NO_VERTEX, NO_VERTEX));
}

VertexList::VertexList(std::vector<std::pair<vertex_t, vertex_t>> neighbors) :
        neighbors(std::move(neighbors)) {}

dimension_t VertexList::getDimension() const {
    return neighbors.size();
}

vertex_t VertexList::next(const vertex_t previous, const vertex_t current) const {
    const std::pair<vertex_t, vertex_t> &currentNeighbors = neighbors.at(current);
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

void VertexList::setNext(const vertex_t previous, const vertex_t current, const vertex_t next) {
    if (previous == next) {
        throw std::runtime_error("setNext: You cannot set both neighbors to the same vertex.");
    }
    const std::pair<vertex_t, vertex_t> &currentNeighbors = neighbors.at(current);
    if (previous != currentNeighbors.first and previous != currentNeighbors.second) {
        throw std::runtime_error(
                "setNext: previous (" + std::to_string(previous) + ") is not a neighbor of current (" +
                std::to_string(current) + ")");
    }
    neighbors.at(current) = std::make_pair(previous, next);
}

void VertexList::addNeighbor(vertex_t vertex, vertex_t newNeighbor) {
    std::pair<vertex_t, vertex_t> &neighborsOfVertex = neighbors.at(vertex);
    if (neighborsOfVertex.first == NO_VERTEX) {
        neighborsOfVertex.first = newNeighbor;
    } else if (neighborsOfVertex.second == NO_VERTEX) {
        neighborsOfVertex.second = newNeighbor;
    } else {
        throw std::runtime_error(
                "addNeighbor: Both neighbors of vertex (" + std::to_string(vertex) + ") are already set.");
    }
}

void VertexList::removeNeighbor(vertex_t vertex, vertex_t neighbor) {
    std::pair<vertex_t, vertex_t> &neighborsOfVertex = neighbors.at(vertex);
    if (neighborsOfVertex.first == neighbor) {
        neighborsOfVertex.first = NO_VERTEX;
    } else if (neighborsOfVertex.second == neighbor) {
        neighborsOfVertex.second = NO_VERTEX;
    } else {
        throw std::runtime_error(
                "removeNeighbor: neighbor (" + std::to_string(neighbor) + ") is not a neighbor of vertex (" +
                std::to_string(vertex) + ").");
    }
}

bool VertexList::makeNeighbors(const vertex_t vertex1, const vertex_t vertex2) {
    const std::pair<vertex_t, vertex_t> neighbors1 = neighbors.at(vertex1);
    const std::pair<vertex_t, vertex_t> neighbors2 = neighbors.at(vertex2);
    if ((neighbors1.first != NO_VERTEX and neighbors1.second != NO_VERTEX) or
        (neighbors2.first != NO_VERTEX and neighbors2.second != NO_VERTEX)) {
        return false;
    } else {
        addNeighbor(vertex1, vertex2);
        addNeighbor(vertex2, vertex1);
        return true;
    }
}


// =================================================== Tour class ======================================================

Tour::Tour() = default;

Tour::Tour(std::vector<std::pair<vertex_t, vertex_t>> neighbors) : VertexList(std::move(neighbors)) {
    if (!isHamiltonianTour()) {
        throw std::runtime_error(
                "A Tour object cannot be initialized with neighbors that do not describe a hamiltonian tour");
    }
}

Tour::Tour(const std::vector<vertex_t> &vertexList) {
    setVertices(vertexList);
}

void Tour::setVertices(const std::vector<vertex_t> &vertexList) {
    neighbors.assign(vertexList.size(), std::make_pair(NO_VERTEX, NO_VERTEX));
    vertex_t previous, current;
    auto it = vertexList.begin();
    current = *it;
    ++it;
    for (; it != vertexList.end(); ++it) {
        previous = current;
        current = *it;
        makeNeighbors(previous, current);
    }
    makeNeighbors(vertexList.back(), vertexList.front());
}

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

Tour Tour::exchange(std::vector<vertex_t> alternatingWalk) {
    Tour tour(*this);
    for (vertex_t i = 0; i < alternatingWalk.size(); ++i) {
        // The current edge in the alternating walk is from vertex1 to vertex2
        vertex_t vertex1 = alternatingWalk.at(i);
        vertex_t vertex2 = alternatingWalk.at((i + 1) % alternatingWalk.size());
        if (i % 2 == 0) { // {x_i, x_i+1} is part of the tour
            tour.removeNeighbor(vertex1, vertex2);
            tour.removeNeighbor(vertex2, vertex1);
        } else { // {x_i, x_i+1} is not part of the tour
            tour.addNeighbor(vertex1, vertex2);
            tour.addNeighbor(vertex2, vertex1);
        }
    }

    if (!tour.isHamiltonianTour()) {
        throw std::runtime_error("The exchange with the given alternating walk does not result in a hamiltonian tour");
    }
}


// =============================================== TourWalker class ====================================================


TourWalker::TourWalker(const Tour &tour, vertex_t first) : TourWalker(tour, first, tour.next(first)) {}

TourWalker::TourWalker(Tour tour, vertex_t first, vertex_t second) : tour(std::move(tour)),
                                                                     current(first), next(second) {}

vertex_t TourWalker::getNextVertex() {
    const vertex_t previous = current;
    current = next;
    next = tour.next(previous, current);
    return previous;
}

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