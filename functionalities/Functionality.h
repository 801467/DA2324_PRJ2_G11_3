#ifndef DA2324_PRJ2_G11_3_FUNCTIONALITY_H
#define DA2324_PRJ2_G11_3_FUNCTIONALITY_H

#include <iostream>
#include "vector"
#include "../utils/haversine.h"
#include "../classes/Graph.h"
#include "../classes/TSPGraph.h"
#include "../utils/MutablePriorityQueue.h"

class Functionality {
public:
    static void backtracking(TSPGraph &graph);
    static void triangularInequality(TSPGraph& graph);
    /*
    static void otherHeuristic(WaterSupply& graph);
    static void notFullyConnected(WaterSupply& graph);*/

private:
    static void
    tspBacktracking(TSPGraph &graph, vector<bool> &visitedVector, Vertex<int>* currNode, vector<int> *currPath,
                    int distance, double cost);
    static vector<Vertex<int>*> prim(TSPGraph& graph);
    static double tspTour(vector<Vertex<int>*> visitOrder);
};


#endif //DA2324_PRJ2_G11_3_FUNCTIONALITY_H
