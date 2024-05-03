//
// Created by jonas on 03/05/2024.
//

#ifndef DA2324_PRJ2_G11_3_TSPGRAPH_H
#define DA2324_PRJ2_G11_3_TSPGRAPH_H

#include "Graph.h"
#include "TSPVertex.h"

using namespace std;

class TSPGraph : public Graph<int> {
public:
    TSPGraph() = default;
    TSPGraph(const TSPGraph& other);
    bool addVertex(int id, double longitude, double latitude);
};


#endif //DA2324_PRJ2_G11_3_TSPGRAPH_H
