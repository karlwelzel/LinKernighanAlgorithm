//
// Created by Karl Welzel on 13.06.19.
//

#ifndef LINKERNINGHANALGORITHM_PRIMSALGORITHM_H
#define LINKERNINGHANALGORITHM_PRIMSALGORITHM_H

#include <vector>
#include <tuple>
#include <functional>

std::tuple<std::vector<std::vector<vertex_t>>, std::vector<vertex_t>>
primsAlgorithm(dimension_t dimension, std::function<distance_t(vertex_t, vertex_t)> dist, vertex_t root);

#endif //LINKERNINGHANALGORITHM_PRIMSALGORITHM_H
