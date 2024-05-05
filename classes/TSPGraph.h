#ifndef DA2324_PRJ2_G11_3_TSPGRAPH_H
#define DA2324_PRJ2_G11_3_TSPGRAPH_H

#include "Graph.h"
#include "TSPVertex.h"

using namespace std;

class TSPGraph : public Graph<int> {
    Vertex<int> *origin{};
    double minCost{};
    vector<int> minPath;
public:
    TSPGraph() = default;
    bool addVertex(int id, double longitude, double latitude);
    Vertex<int> *getOrigin() { return origin; };
    void setOrigin(int id) { origin = findVertex(id); };
    double getMinCost() const { return minCost; };
    void setMinCost(double cost) { minCost = cost; };
    vector<int> getMinPath() { return minPath; };
    void setMinPath(vector<int> path) { minPath = std::move(path); };
};


#endif //DA2324_PRJ2_G11_3_TSPGRAPH_H
