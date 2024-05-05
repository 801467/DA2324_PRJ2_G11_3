#include "TSPGraph.h"

bool TSPGraph::addVertex(int id, double longitude, double latitude) {
    if (findVertex(id) != nullptr)
        return false;
    vertexSet.push_back(new TSPVertex(id,longitude,latitude));
    return true;
}