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

vertex_t VertexList::next(const vertex_t previous, const vertex_t current) const {
    const std::vector<vertex_t> &currentNeighbors = neighbors.at(current);
    if (currentNeighbors.size() != 2) {
        throw std::runtime_error("current (" + std::to_string(current) + ") does not have two neighbors");
    } else if (currentNeighbors[0] == previous) {
        return currentNeighbors[1];
    } else if (currentNeighbors[1] == previous) {
        return currentNeighbors[0];
    } else {
        throw std::runtime_error(
                "previous (" + std::to_string(previous) + ") is not a neighbor of current (" +
                std::to_string(current) + ")");
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
    const std::vector<vertex_t> &vertexNeighbors = neighbors.at(vertex);
    return std::find(vertexNeighbors.begin(), vertexNeighbors.end(), neighbor) != vertexNeighbors.end();
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
    const auto iterator = std::find(vertexNeighbors.begin(), vertexNeighbors.end(), neighbor);
    if (iterator == vertexNeighbors.end()) {
        throw std::runtime_error(
                "removeNeighbor: neighbor (" + std::to_string(neighbor) + ") is not a neighbor of vertex (" +
                std::to_string(vertex) + ").");
    } else {
        vertexNeighbors.erase(iterator);
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


// =================================================== Tour class ======================================================

Tour::Tour() = default;

Tour::Tour(std::vector<std::vector<vertex_t>> neighbors) : VertexList(std::move(neighbors)) {
    if (!isHamiltonianTour()) {
        throw std::runtime_error(
                "A Tour object cannot be initialized with neighbors that do not describe a hamiltonian tour");
    }
}

Tour::Tour(const std::vector<vertex_t> &vertexList) {
    setVertices(vertexList);
}

void Tour::setVertices(const std::vector<vertex_t> &vertexList) {
    neighbors.assign(vertexList.size(), std::vector<vertex_t>());
    for (auto it = vertexList.begin(); it != vertexList.end() - 1; ++it) {
        makeNeighbors(*it, *std::next(it));
    }
    makeNeighbors(vertexList.back(), vertexList.front());
}

bool Tour::isHamiltonianTour() const {
    dimension_t dimension = neighbors.size();
    std::vector<bool> visited(dimension, false);
    TourWalker tourWalker(*this, 0);
    vertex_t currentVertex = tourWalker.getNextVertex();
    do {
        if (visited[currentVertex]) {
            return false;
        }
        visited[currentVertex] = true;
        currentVertex = tourWalker.getNextVertex();
    } while (currentVertex != 0);
    // Check if all elements in visited are true
    return std::all_of(visited.begin(), visited.end(), [](bool v) { return v; });
}

Tour Tour::exchange(const std::vector<vertex_t> &alternatingWalk) const {
    Tour tour(*this);
    // Remove all edges in the alternatingWalk that are part of the tour (every edge with even i)
    for (vertex_t i = 0; i < alternatingWalk.size() - 1; i += 2) {
        // The current edge in the alternating walk is from vertex1 to vertex2
        vertex_t vertex1 = alternatingWalk.at(i);
        vertex_t vertex2 = alternatingWalk.at(i + 1);
        tour.removeNeighbor(vertex1, vertex2);
        tour.removeNeighbor(vertex2, vertex1);
    }
    // Add all edges in the alternatingWalk that are not part of the tour (every edge with odd i)
    for (vertex_t i = 1; i < alternatingWalk.size() - 1; i += 2) {
        // The current edge in the alternating walk is from vertex1 to vertex2
        vertex_t vertex1 = alternatingWalk.at(i);
        vertex_t vertex2 = alternatingWalk.at(i + 1);
        tour.addNeighbor(vertex1, vertex2);
        tour.addNeighbor(vertex2, vertex1);
    }
    return tour;
}

bool Tour::isTourAfterExchange(const std::vector<vertex_t> &alternatingWalk) const {
    return exchange(alternatingWalk).isHamiltonianTour();
}


// =============================================== TourWalker class ====================================================


TourWalker::TourWalker(const Tour &tour, vertex_t first) : TourWalker(tour, first, tour.next(first)) {}

TourWalker::TourWalker(const Tour &tour, vertex_t first, vertex_t second) : tour(tour), current(first), next(second) {}

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
    out << output.substr(0, output.length() - 2);
    return out;
}