//
// Created by Karl Welzel on 29.03.2019.
//

// This library contains the code for the introduction assignment

#include <algorithm>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include "Tour.h"
#include "TsplibUtils.h"
#include "SimpleHeuristic.h"


// ================================================ VertexList class ===================================================

VertexList::VertexList() = default;

VertexList::VertexList(const dimension_t dimension) : neighbors(dimension) {}

VertexList::VertexList(std::vector<std::vector<vertex_t>> neighbors) :
        neighbors(std::move(neighbors)) {}

dimension_t VertexList::getDimension() const {
    return neighbors.size();
}

std::vector<vertex_t> VertexList::getNeighbors(vertex_t vertex) {
    return neighbors.at(vertex);
}

int VertexList::neighborIndex(vertex_t vertex, vertex_t neighbor) const {
    const std::vector<vertex_t> &vertexNeighbors = neighbors.at(vertex);
    if (!vertexNeighbors.empty() and vertexNeighbors[0] == neighbor) {
        return 0;
    } else if (vertexNeighbors.size() >= 2 and vertexNeighbors[1] == neighbor) {
        return 1;
    } else {
        return -1;
    }
}

vertex_t VertexList::next(const vertex_t previous, const vertex_t current) const {
    const std::vector<vertex_t> &currentNeighbors = neighbors.at(current);
    if (currentNeighbors.size() != 2) {
        throw std::runtime_error("current (" + std::to_string(current) + ") does not have two neighbors");
    } else {
        int index = neighborIndex(current, previous);
        if (index == -1) {
            throw std::runtime_error("previous (" + std::to_string(previous) + ") is not a neighbor of current (" +
                                     std::to_string(current) + ")");
        } else {
            return currentNeighbors[(index + 1) % 2];
        }
    }
}

vertex_t VertexList::next(const vertex_t current) const {
    const std::vector<vertex_t> &currentNeighbors = neighbors.at(current);
    if (currentNeighbors.empty()) {
        throw std::runtime_error("current (" + std::to_string(current) + ") does not have a neighbor");
    }
    return currentNeighbors[0];
}

bool VertexList::isNeighbor(vertex_t vertex, vertex_t neighbor) const {
    return neighborIndex(vertex, neighbor) != -1;
}

bool VertexList::containsEdge(vertex_t vertex1, vertex_t vertex2) const {
    return isNeighbor(vertex1, vertex2);
}

void VertexList::addNeighbor(vertex_t vertex, vertex_t newNeighbor) {
    std::vector<vertex_t> &vertexNeighbors = neighbors.at(vertex);
    if (vertexNeighbors.size() >= 2) {
        throw std::runtime_error(
                "Both neighbors of vertex (" + std::to_string(vertex) + ") are already set.");
    } else {
        vertexNeighbors.push_back(newNeighbor);
    }
}

void VertexList::removeNeighbor(vertex_t vertex, vertex_t neighbor) {
    std::vector<vertex_t> &vertexNeighbors = neighbors.at(vertex);
    int index = neighborIndex(vertex, neighbor);
    if (index == -1) {
        throw std::runtime_error(
                "removeNeighbor: neighbor (" + std::to_string(neighbor) + ") is not a neighbor of vertex (" +
                std::to_string(vertex) + ").");
    } else {
        vertexNeighbors.erase(vertexNeighbors.begin() + index);
    }
}

bool VertexList::makeNeighbors(const vertex_t vertex1, const vertex_t vertex2) {
    if (neighbors.at(vertex1).size() >= 2 or neighbors.at(vertex2).size() >= 2) {
        return false;
    } else {
        addNeighbor(vertex1, vertex2);
        addNeighbor(vertex2, vertex1);
        return true;
    }
}

std::vector<vertex_t> VertexList::toTourSequence() const {
    std::vector<vertex_t> result = std::vector<vertex_t>();
    vertex_t previous = 0;
    vertex_t current = next(0);
    vertex_t next_;
    do {
        result.push_back(previous);
        next_ = next(previous, current);
        previous = current;
        current = next_;
    } while (previous != 0);
    return result;
}

// ============================================== TourParts class ======================================================

TourParts::TourParts(const dimension_t dimension) : VertexList(dimension) {
    for (dimension_t i = 0; i < dimension; ++i) {
        parent.push_back(i); // parent[i] = i
        size.push_back(1);   // size[i] = 1
        pathEnds.emplace_back(i, i);
    }
}

vertex_t TourParts::find(const vertex_t x) {
    if (parent.at(x) != x) { // path compression
        parent.at(x) = find(parent.at(x));
    }
    return parent.at(x);
}

void TourParts::join(const vertex_t x, const vertex_t y) {
    vertex_t xr = find(x); // root of x
    vertex_t yr = find(y); // root of y
    if (xr == yr) {
        return; // A path cannot be joined with itself, this would lead to non-Hamiltonian tours
    }
    std::pair<vertex_t, vertex_t> newPathEnds;
    if (x == pathEnds.at(xr).first and y == pathEnds.at(yr).first) {
        newPathEnds = std::make_pair(pathEnds.at(xr).second, pathEnds.at(yr).second);
    } else if (x == pathEnds.at(xr).first and y == pathEnds.at(yr).second) {
        newPathEnds = std::make_pair(pathEnds.at(xr).second, pathEnds.at(yr).first);
    } else if (x == pathEnds.at(xr).second and y == pathEnds.at(yr).first) {
        newPathEnds = std::make_pair(pathEnds.at(xr).first, pathEnds.at(yr).second);
    } else if (x == pathEnds.at(xr).second and y == pathEnds.at(yr).second) {
        newPathEnds = std::make_pair(pathEnds.at(xr).first, pathEnds.at(yr).first);
    } else {
        return; // x and y cannot be joined, because one of them is not an end of its path
    }

    if (!makeNeighbors(x, y)) {
        throw std::runtime_error("You cannot join neighbors when one of them already has two neighbors.");
    }

    if (size.at(xr) > size.at(yr)) {
        std::swap(xr, yr);
    }
    parent.at(xr) = yr;
    size.at(yr) += size.at(xr);
    pathEnds.at(yr) = newPathEnds;
}

bool TourParts::isTourClosable() {
    return size.at(find(0)) == getDimension();
}

Tour TourParts::closeTour() {
    // if there is only one root, then its also the root of 0
    const vertex_t root = find(0);
    if (isTourClosable()) {
        // join the last ends to form a tour
        if (!makeNeighbors(pathEnds.at(root).first, pathEnds.at(root).second)) {
            throw std::runtime_error("You cannot join neighbors when one of them already has two neighbors.");
        }
        return Tour(toTourSequence());
    } else {
        throw std::runtime_error("The tour cannot be closed while there are multiple connected components");
    }
}


Tour simpleHeuristic(const TsplibProblem &tsplibProblem) {
    std::vector<std::tuple<distance_t, vertex_t, vertex_t>> edges;
    for (vertex_t i = 0; i < tsplibProblem.getDimension(); ++i) {
        for (vertex_t j = i + 1; j < tsplibProblem.getDimension(); ++j) {
            edges.emplace_back(tsplibProblem.dist(i, j), i, j);
        }
    }

    // Sort edges by edge cost/distance in ascending order
    std::sort(edges.begin(), edges.end());

    TourParts tourParts(tsplibProblem.getDimension());
    for (const auto edge : edges) {
        vertex_t vertex1, vertex2;
        std::tie(std::ignore, vertex1, vertex2) = edge;
        tourParts.join(vertex1, vertex2);
        if (tourParts.isTourClosable()) {
            return tourParts.closeTour();
        }
    }

    // This line will never be reached, because the graph is a complete graph
    throw std::runtime_error("SimpleHeuristic: Although every edge was checked, no tour has been found.");
}

Tour ascendingVerticesHeuristic(const TsplibProblem &tsplibProblem) {
    std::vector<vertex_t> ascendingVertices(tsplibProblem.getDimension());
    std::iota(ascendingVertices.begin(), ascendingVertices.end(), 0);
    return Tour(ascendingVertices);
}
