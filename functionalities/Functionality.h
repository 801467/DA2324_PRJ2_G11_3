#ifndef DA2324_PRJ2_G11_3_FUNCTIONALITY_H
#define DA2324_PRJ2_G11_3_FUNCTIONALITY_H

#include <iostream>
#include "vector"
#include <numeric>
#include "../utils/haversine.h"
#include "../classes/Graph.h"
#include "../classes/TSPGraph.h"
#include "../utils/MutablePriorityQueue.h"

class Functionality {
public:
    static void backtracking(TSPGraph &graph);
    static void triangularInequality(TSPGraph &graph);
    static void nearestNeighbour(TSPGraph &graph);
    static void backtrackedNearestNeighbour(TSPGraph &graph, int originId, int destinationId);
    static bool checkHamiltonianFeasibility(TSPGraph &graph, int originId);
    static void runClusteredGraphs(TSPGraph &graph, int originId);

private:
    static void
    tspBacktracking(TSPGraph &graph, vector<bool> &visitedVector, Vertex<int> *currNode, vector<int> *currPath,
                    int distance, double cost);
    static void prim(TSPGraph &graph);
    static vector<Vertex<int> *> preOrderWalk(TSPGraph &graph);
    static void tspDfsVisit(Vertex<int>* v, vector<Vertex<int>*> &orderedPath);
    static double tspTour(vector<Vertex<int> *> visitOrder);
    static void
    tspBacktrackingNearestNeighbour(TSPGraph &graph, vector<bool> &visitedVector, Vertex<int> *currNode,
                                    vector<int> *currPath, int distance,
                                    double cost, Vertex<int> *destinationNode);
    static vector<TSPGraph> generateClusters(TSPGraph &graph);
    static TSPGraph generateSubgraph(TSPGraph &originalGraph, int fromPosition, int toPosition);
    static bool isValidSubgraph(TSPGraph &originalGraph, TSPGraph &subgraph);
};


#endif //DA2324_PRJ2_G11_3_FUNCTIONALITY_H
