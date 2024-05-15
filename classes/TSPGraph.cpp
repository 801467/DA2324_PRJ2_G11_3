#include "TSPGraph.h"

bool TSPGraph::addVertex(int id, double longitude, double latitude) {
    if (findVertex(id) != nullptr)
        return false;
    vertexSet.push_back(new Vertex<int>(id,longitude,latitude));
    vertexIndex[id] = vertexSet.size() - 1;
    return true;
}
