//
// Created by jonas on 03/05/2024.
//

#include "TSPGraph.h"

TSPGraph::TSPGraph(const TSPGraph& other)  : Graph(other) {
    if (this != &other) {
        vertexSet.clear();
        for (Vertex<int>* vertex : other.vertexSet) {
            auto newVertexPtr = new Vertex<int>(vertex->getInfo());
            vertexSet.push_back(newVertexPtr);
        }
        for (Vertex<int>* vertex : other.vertexSet){
            for (Edge<int>* e : vertex->getAdj()){
                addEdge(vertex->getInfo(),e->getDest()->getInfo(),e->getWeight());
            }
        }
    }
}

bool TSPGraph::addVertex(int id, double longitude, double latitude) {
    if (findVertex(id) != nullptr)
        return false;
    vertexSet.push_back(new TSPVertex(id,longitude,latitude));
    return true;
}