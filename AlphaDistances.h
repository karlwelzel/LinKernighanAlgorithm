//
// Created by Karl Welzel on 26.06.19.
//

#ifndef LINKERNINGHANALGORITHM_ALPHADISTANCES_H
#define LINKERNINGHANALGORITHM_ALPHADISTANCES_H


#include <functional>
#include <vector>
#include "Tour.h"

vertex_t secondNearestNeighbor(vertex_t v, std::vector<vertex_t> &vertices,
                               const std::function<distance_t(vertex_t, vertex_t)> &dist);

std::vector<std::vector<distance_t>>
betaValues(dimension_t dimension, const std::function<distance_t(vertex_t, vertex_t)> &dist);

std::vector<std::vector<distance_t>>
alphaDistances(dimension_t dimension, const std::function<distance_t(vertex_t, vertex_t)> &dist);

#endif //LINKERNINGHANALGORITHM_ALPHADISTANCES_H
