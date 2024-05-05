#ifndef DA2324_PRJ2_G11_3_FUNCTIONALITY_H
#define DA2324_PRJ2_G11_3_FUNCTIONALITY_H

#include <iostream>
#include "vector"
#include "../classes/Graph.h"
#include "../classes/TSPGraph.h"

class Functionality {
public:
    static void backtracking(TSPGraph &graph);
    /*
    static void triangularInequality(WaterSupply& graph);
    static void otherHeuristic(WaterSupply& graph);
    static void notFullyConnected(WaterSupply& graph);*/

private:
    static void
    tspBacktracking(TSPGraph &graph, vector<bool> &visitedVector, Vertex<int> *currNode, vector<int> *currPath,
                    int distance, double cost);
};


#endif //DA2324_PRJ2_G11_3_FUNCTIONALITY_H
