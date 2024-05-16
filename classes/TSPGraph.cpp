#include "TSPGraph.h"

bool TSPGraph::addVertex(int id, double longitude, double latitude) {
    if (findVertex(id) != nullptr)
        return false;
    vertexSet.push_back(new Vertex<int>(id,longitude,latitude));
    vertexIndex[id] = vertexSet.size() - 1;
    return true;
}

vector<vector<double>> TSPGraph::getAdjacencyMatrix() const {
    // Get the number of vertices
    auto vertices = getVertexSet();
    size_t numVertices = vertices.size();

    // Initialize the adjacency matrix with infinity for unreachable vertices
    vector<vector<double>> adjacencyMatrix(numVertices, vector<double>(numVertices, INF));

    // Iterate through the vertices
    for (size_t i = 0; i < numVertices; ++i) {
        // Get the current vertex
        auto& currentVertex = vertices[i];
        // Set the diagonal element to zero (distance to itself)
        adjacencyMatrix[i][i] = 0;

        // Iterate through the edges of the current vertex
        for (const auto& edge : currentVertex->getAdj()) {
            // Get the destination vertex index
            size_t destIndex = distance(vertices.begin(), find(vertices.begin(), vertices.end(), edge->getDest()));

            // Set the distance between the current vertex and its destination
            // You need to implement a method in the Edge class to get the weight
            adjacencyMatrix[i][destIndex] = edge->getWeight();
        }
    }

    return adjacencyMatrix;
}

void TSPGraph::clearState() {
    origin = nullptr;
    minCost = INF;
    minPath.clear();

    // vertex
    for(auto v : vertexSet){
        v->setVisited(false);
        v->setProcessing(false);
        v->setDist(0);
        v->setPath(nullptr);
        v->setIncFlow(0);
    }
}
