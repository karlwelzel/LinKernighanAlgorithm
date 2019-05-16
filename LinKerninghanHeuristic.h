//
// Created by karl on 14.05.19.
//

#ifndef LINKERNINGHANALGORITHM_LINKERNINGHANHEURISTIC_H
#define LINKERNINGHANALGORITHM_LINKERNINGHANHEURISTIC_H

#include "Tour.h"
#include "TsplibUtils.h"

bool containsEdge(const std::vector<vertex_t> &walk, vertex_t vertex1, vertex_t vertex2);

Tour linKerninghanHeuristic(const TsplibProblem &tsplibProblem, const Tour &startTour);

#endif //LINKERNINGHANALGORITHM_LINKERNINGHANHEURISTIC_H
