#include "TSPGraph.h"
#include "haversine.h"

bool TSPGraph::addVertex(int id, double longitude, double latitude) {
    if (findVertex(id) != nullptr)
        return false;
    vertexSet.push_back(new Vertex<int>(id, longitude, latitude));
    vertexIndex[id] = vertexSet.size() - 1;
    return true;
}

/**
 * @brief Generates the adjacency matrix for the TSPGraph.
 * 
 * This function generates the adjacency matrix for the graph. The adjacency matrix is a 2D vector where the element at the i-th row and j-th column represents the weight of the edge between the i-th and j-th vertices.
 * If there's no direct edge between the i-th and j-th vertices, the element is set to INF (representing "infinity").
 * The diagonal elements of the matrix (representing the distance of a vertex to itself) are set to 0.
 * 
 * @return A 2D vector representing the adjacency matrix of the graph.
 */
vector<vector<double>> TSPGraph::getAdjacencyMatrix() const {
    // Get the number of vertices
    auto vertices = getVertexSet();
    size_t numVertices = vertices.size();

    // Initialize the adjacency matrix with infinity for unreachable vertices
    vector<vector<double>> adjacencyMatrix(numVertices, vector<double>(numVertices, INF));

    // Iterate through the vertices
    for (size_t i = 0; i < numVertices; i++) {
        // Get the current vertex
        auto &currentVertex = vertices[i];
        // Set the diagonal element to zero (distance to itself)
        adjacencyMatrix[i][i] = 0;

        // Iterate through the edges of the current vertex
        for (const auto &edge: currentVertex->getAdj()) {
            // Get the destination vertex index
            size_t destIndex = vertexIndex.at(edge->getDest()->getInfo());

            // Set the distance between the current vertex and its destination
            adjacencyMatrix[i][destIndex] = edge->getWeight();
        }

        // Assume all other connections are possible at a Haversine distance
        for (size_t n = 0; n < numVertices; n++) {
            if (adjacencyMatrix[i][n] != INF)
                continue;

            auto &nextVertex = vertices[n];
            adjacencyMatrix[i][n] =
                    1000 * calculate_distance(currentVertex->getLatitude(), currentVertex->getLongitude(),
                                              nextVertex->getLatitude(), nextVertex->getLongitude());

        }
    }

    return adjacencyMatrix;
}

void TSPGraph::clearState() {
    origin = nullptr;
    minCost = INF;
    minPath.clear();
    forceStop = false;

    // vertex
    for (auto v: vertexSet) {
        v->setVisited(false);
        v->setProcessing(false);
        v->setDist(0);
        v->setPath(nullptr);
        v->setIncFlow(0);
    }
}

void TSPGraph::reorderByGeographicalPosition() {
    // sort vertexSet
    sort(vertexSet.begin(), vertexSet.end(), [](Vertex<int> *a, Vertex<int> *b) {
        // distinguish opposite coordinates, by giving more weight to x-axis.
        return (a->getLatitude() * 2) + a->getLongitude() < (b->getLatitude() * 2) + b->getLongitude();
    });

    // recreate map
    vertexIndex.clear();
    for(int i = 0; i < vertexSet.size(); i++){
        auto v = vertexSet[i];
        vertexIndex[v->getInfo()] = i;
    }
}
