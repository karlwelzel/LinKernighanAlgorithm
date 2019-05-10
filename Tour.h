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

class VertexList {
protected:
    std::vector<std::pair<vertex_t, vertex_t>> neighbors;

    void setNext(vertex_t previous, vertex_t current, vertex_t next);

public:
    static const vertex_t NO_VERTEX;

    VertexList();

    explicit VertexList(dimension_t dimension);

    explicit VertexList(std::vector<std::pair<vertex_t, vertex_t>> neighbors);

    dimension_t getDimension() const;

    vertex_t next(vertex_t previous, vertex_t current) const;

    vertex_t next(vertex_t current) const;

    bool makeNeighbors(vertex_t vertex1, vertex_t vertex2);
};

class Tour : public VertexList {
protected:
    void setVertices(const std::vector<vertex_t> &vertexList);

public:
    Tour();

    explicit Tour(std::vector<std::pair<vertex_t, vertex_t>> neighbors);

    explicit Tour(const std::vector<vertex_t> &vertexList);

    bool isHamiltonianTour() const;
};

class TourWalker {
private:
    const Tour tour;
    vertex_t current;
    vertex_t next;

public:
    TourWalker(const Tour &tour, vertex_t first);

    TourWalker(Tour tour, vertex_t first, vertex_t second);

    vertex_t getNextVertex();
};

std::ostream &operator<<(std::ostream &out, const Tour &tour);

#endif //LINKERNINGHANALGORITHM_TOUR_H
