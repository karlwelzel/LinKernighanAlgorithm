//
// Created by Karl Welzel on 19/04/2019.
//

#ifndef LINKERNINGHANALGORITHM_TOUR_H
#define LINKERNINGHANALGORITHM_TOUR_H

#include <limits>
#include <map>
#include <list>
#include <iostream>

using vertex_t = size_t; // The type used for vertices
using dimension_t = vertex_t; // The type used for counting vertices
using distance_t = unsigned long; // The type used for distances and lengths

class VertexList {
protected:
    std::map<vertex_t, std::pair<vertex_t, vertex_t>> neighbors;

public:
    vertex_t NO_VERTEX = std::numeric_limits<unsigned int>::max();

    VertexList();

    explicit VertexList(std::map<vertex_t, std::pair<vertex_t, vertex_t>> neighbors);

    vertex_t next(vertex_t previous, vertex_t current) const;

    vertex_t next(vertex_t current) const;

    void setNext(vertex_t previous, vertex_t current, vertex_t next);

    bool makeNeighbors(vertex_t vertex1, vertex_t vertex2);
};

class Tour : public VertexList {
protected:
    void setVertices(const std::list<vertex_t> &vertexList);

public:
    using VertexList::VertexList;

    bool isHamiltonianTour();
};

class TourWalker {
private:
    Tour tour;
    vertex_t current;
    vertex_t next;

public:
    TourWalker(const Tour &tour, vertex_t first);

    TourWalker(Tour tour, vertex_t first, vertex_t second);

    vertex_t getNextVertex();
};

std::ostream &operator<<(std::ostream &out, const Tour &tour);

#endif //LINKERNINGHANALGORITHM_TOUR_H
