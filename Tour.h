//
// Created by Karl Welzel on 19/04/2019.
//

#ifndef LINKERNINGHANALGORITHM_TOUR_H
#define LINKERNINGHANALGORITHM_TOUR_H

#include <limits>
#include <utility>
#include <regex>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <list>
#include <vector>
#include <cmath>
#include <iostream>
#include <utility>
#include "TsplibProblem.h"

class VertexList {
protected:
    std::map<unsigned int, std::pair<unsigned int, unsigned int>> neighbors;

public:
    unsigned int NO_VERTEX = std::numeric_limits<unsigned int>::max();

    VertexList();

    explicit VertexList(std::map<unsigned int, std::pair<unsigned int, unsigned int>> neighbors);

    unsigned int next(unsigned int previous, unsigned int current) const;

    unsigned int next(unsigned int current) const;

    void setNext(unsigned int previous, unsigned int current, unsigned int next);

    bool makeNeighbors(unsigned int vertex1, unsigned int vertex2);
};

class Tour : public VertexList {
protected:
    void setVertices(const std::list<unsigned int> &vertexList);

public:
    using VertexList::VertexList;

    unsigned int length(TsplibProblem &tsplibProblem);

    bool isHamiltonianTour();
};

class TourWalker {
private:
    Tour tour;
    unsigned int current;
    unsigned int next;

public:
    TourWalker(const Tour &tour, unsigned int first);

    TourWalker(Tour tour, unsigned int first, unsigned int second);

    unsigned int getNextVertex();
};

std::ostream &operator<<(std::ostream &out, const Tour &tour);

#endif //LINKERNINGHANALGORITHM_TOUR_H
