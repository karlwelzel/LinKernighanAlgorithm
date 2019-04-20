//
// Created by Karl Welzel on 29/03/2019.
//

#ifndef LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H
#define LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H

#include <utility>
#include <vector>
#include "Tour.h"
#include "TsplibUtils.h"

class TourParts : public VertexList {
private:
    std::vector<vertex_t> parent;
    std::vector<dimension_t> size;
    std::vector<std::pair<vertex_t, vertex_t>> pathEnds;

public:
    explicit TourParts(dimension_t dimension);

    vertex_t find(vertex_t x);

    void join(vertex_t x, vertex_t y);

    bool isTourClosable();

    Tour closeTour();
};

Tour simpleHeuristic(const TsplibProblem &tsplibProblem);

#endif //LINKERNINGHANALGORITHM_SIMPLEHEURISTIC_H
