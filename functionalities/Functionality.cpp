#include <unordered_set>
#include "Functionality.h"

void Functionality::backtracking(TSPGraph &graph) {
    cout << "Running Backtracking..." << endl;
    cout << endl;

    graph.clearState();

    unsigned int nodesSize = graph.getNumVertex();
    vector<bool> visitedVector(nodesSize);
    for (int i = 0; i < nodesSize; i++)
        visitedVector[i] = false;

    // start from "origin"
    graph.setOrigin(0);
    visitedVector[0] = true;

    // set initial cost to max number possible
    graph.setMinCost(std::numeric_limits<double>::max());
    vector<int> currPath;

    // recursively branch out,
    // keep track of visited nodes, distance and min cost found
    tspBacktracking(graph, visitedVector, graph.findVertex(0), &currPath, 1, 0);

    cout << "Min Cost: " << graph.getMinCost() << endl;
    cout << "Path: " << graph.getOrigin()->getInfo() << " ";
    unsigned int i = 0;
    for (auto element: graph.getMinPath()) {
        i++;
        (i % 20 == 0) ? cout << endl : cout << element << " ";
    }
    cout << graph.getOrigin()->getInfo() << endl << endl;
}

void Functionality::tspBacktracking(TSPGraph &graph, vector<bool> &visitedVector, Vertex<int> *currNode,
                                    vector<int> *currPath, int distance,
                                    double cost) {
    if (distance == graph.getNumVertex()) {
        for (auto adj: currNode->getAdj()) {
            // If "last node" is same as source
            if (adj->getDest() == graph.getOrigin()) {
                // check if cost is less than previously found
                double currCost = cost + adj->getWeight();
                if (graph.getMinCost() > currCost) {
                    graph.setMinCost(currCost);
                    graph.setMinPath(*currPath);
                }
            }
        }
        return;
    }

    // Backtracking Part
    // Loop equally through adjacent nodes
    for (auto adj: currNode->getAdj()) {
        // branch filtering based on cost
        if (cost + adj->getWeight() >= graph.getMinCost()) {
            continue;
        }

        int destId = adj->getDest()->getInfo();
        int destPosition = graph.getVertexIndex()[destId];

        // if not visited
        if (!visitedVector[destPosition]) {
            // mark as visited
            visitedVector[destPosition] = true;
            currPath->push_back(destId);

            // Then recursively run this method on the adjacent node,
            // until all paths have either failed or found the source
            tspBacktracking(graph, visitedVector, adj->getDest(), currPath, distance + 1,
                            cost + adj->getWeight());

            // mark as unvisited when backtracking
            visitedVector[destPosition] = false;
            currPath->pop_back();
        }
    }
}


void Functionality::triangularInequality(TSPGraph &graph) {
    cout << "Running Triangular Approximation Heuristic..." << endl;
    cout << endl;

    graph.clearState();
    double cost;

    // select root vertex
    graph.setOrigin(0);

    //compute MST from origin
    prim(graph);

    //perform a pre-order walk of the MST
    vector<Vertex<int> *> path = preOrderWalk(graph);

    //perform the travel
    cost = tspTour(path);

    graph.setMinCost(cost);

    cout << "Min Cost: " << graph.getMinCost() << endl;
    cout << "Path: ";
    unsigned int i = 0;
    for (auto element: path) {
        i++;
        (i % 20 == 0) ? cout << endl : cout << element->getInfo() << " ";
    }

    cout << path.front()->getInfo() << endl << endl;

}


void Functionality::prim(TSPGraph &graph) {
    vector<Vertex<int> *> path;
    MutablePriorityQueue<Vertex<int>> q;

    for (Vertex<int> *v: graph.getVertexSet()) {
        v->setVisited(false);
        v->setPath(nullptr);
        v->setDist(INF);
    }
    Vertex<int> *origin = graph.getOrigin();
    origin->setDist(0);
    q.insert(origin);
    while (!q.empty()) {
        Vertex<int> *u = q.extractMin();
        u->setVisited(true);
        path.push_back(u);
        for (Edge<int> *e: u->getAdj()) {
            Vertex<int> *w = e->getDest();
            if (!w->isVisited()) {
                double oldDist = w->getDist();
                if (e->getWeight() < oldDist) {
                    w->setDist(e->getWeight());
                    w->setPath(e);
                    if (oldDist == INF) {
                        q.insert(w);
                    } else {
                        q.decreaseKey(w);
                    }
                }
            }
        }
    }
}

vector<Vertex<int> *> Functionality::preOrderWalk(TSPGraph &graph) {
    vector<Vertex<int> *> orderedPath;
    Vertex<int> *origin = graph.getOrigin();

    for (auto v: graph.getVertexSet()) {
        v->setVisited(false);
    }

    tspDfsVisit(origin, orderedPath);


    return orderedPath;
}

void Functionality::tspDfsVisit(Vertex<int> *v, vector<Vertex<int> *> &orderedPath) {
    v->setVisited(true);
    orderedPath.push_back(v);
    for (auto e: v->getAdj()) {
        auto w = e->getDest();
        if (!w->isVisited() && (w->getPath() == e)) tspDfsVisit(w, orderedPath);
    }
}

double Functionality::tspTour(vector<Vertex<int> *> path) {
    double cost = 0;
    bool connected = false;
    for (int i = 0; i < path.size() - 1; i++) {
        Vertex<int> *nextVertex = path[i + 1];
        Vertex<int> *currVertex = path[i];
        for (auto e: currVertex->getAdj()) {
            if (e->getDest() == nextVertex) { // if there is an edge connecting both nodes, simply add the edge cost
                connected = true;
                cost += e->getWeight();
                break;
            }
        }
        if (!connected) {    // if not, calculate the distance between them and add to the cost
            cost += 1000 * calculate_distance(currVertex->getLatitude(), currVertex->getLongitude(),
                                              nextVertex->getLatitude(), nextVertex->getLongitude());
        }
        connected = false;
    }

    //  go from the last node to the origin
    Vertex<int> *nextVertex = path.front();
    Vertex<int> *currVertex = path.back();
    for (auto e: currVertex->getAdj()) {
        if (e->getDest() == nextVertex) {
            connected = true;
            cost += e->getWeight();
        }
    }
    if (!connected)
        cost += 1000 * calculate_distance(currVertex->getLatitude(), currVertex->getLongitude(),
                                          nextVertex->getLatitude(), nextVertex->getLongitude());

    return cost;
}


void Functionality::nearestNeighbour(TSPGraph &graph) {
    cout << "Running Nearest Neighbour Heuristic..." << endl;
    cout << endl;

    graph.clearState();

    // select root vertex
    int origin = 0;
    graph.setOrigin(origin);

    double cost = 0;
    vector<int> path;

    auto nodes = graph.getVertexSet();
    auto distance = graph.getAdjacencyMatrix();
    auto neighbor = graph.getOrigin();
    auto start_node_index = graph.getVertexIndex()[origin];

    int neighbor_index;
    double MIN;

    while (true) {
        neighbor->setVisited(true);
        path.push_back(neighbor->getInfo());
        neighbor_index = graph.getVertexIndex()[neighbor->getInfo()];
        MIN = INF;

        for (int i = 0; i < nodes.size(); i++) {
            // look for unvisited nodes
            if (!nodes[i]->isVisited()) {
                if (distance[neighbor_index][i] < MIN) {
                    MIN = min(distance[neighbor_index][i], MIN);
                    neighbor = nodes[i];
                }
            }
        }
        // if all nodes were already visited
        if (MIN == INF)
            break;

        cost += MIN;
    }
    // add cost of last node to origin
    cost += distance[path.back()][start_node_index];

    cout << "Min Cost: " << cost << endl;
    cout << "Path: ";
    unsigned int i = 0;
    for (auto element: path) {
        i++;
        (i % 20 == 0) ? cout << endl : cout << element << " ";
    }
    cout << path.front() << endl << endl;
}

void Functionality::runClusteredGraphs(TSPGraph &graph, int originId) {
    cout << "Running Clustered, Backtracked with Nearest Neighbour Heuristic..." << endl;
    cout << endl;

    graph.clearState();
    graph.setOrigin(originId);

    // before clustering, sort by geographical position
    graph.reorderByGeographicalPosition();

    auto clusters = generateClusters(graph);
    cout << "Generated a total of " << clusters.size() << " clusters. Targeted less than 100 nodes." << endl << endl;

    double cost = 0;
    vector<int> path;
    for (auto cluster: clusters) {
        int origin = cluster.getOrigin()->getInfo();
        int destination = cluster.getVertexSet().back()->getInfo();
        backtrackedNearestNeighbour(cluster, origin, destination);

        // TODO add transition costs
        cost += cluster.getMinCost();
        auto clusterPath = cluster.getMinPath();
        path.insert(path.end(), clusterPath.begin(), clusterPath.end());
    }

    cout << "Min Cost: " << cost << " (+ transitions)" << endl;
    cout << "Path: " << originId << " ";
    unsigned int i = 0;
    for (auto element: path) {
        i++;
        (i % 30 == 0) ? cout << endl : cout << element << " ";
    }
    cout << originId << endl << endl;
}

vector<TSPGraph> Functionality::generateClusters(TSPGraph &graph) {
    vector<TSPGraph> clusters;

    // Avoid unnecessary smaller clustering
    if (graph.getVertexSet().size() < 101) {
        clusters.push_back(graph);
        return clusters;
    }

    // split in exact half
    int median = graph.getNumVertex() / 2;
    TSPGraph graph1 = generateSubgraph(graph, 0, median - 1);
    TSPGraph graph2 = generateSubgraph(graph, median, graph.getNumVertex() - 1);

    // TODO: check frontier
    if (!isValidSubgraph(graph, graph1) || !isValidSubgraph(graph, graph2)) {
        int firstId = graph.getVertexSet().front()->getInfo();
        int lastId = graph.getVertexSet().back()->getInfo();
        cout << "Rolled back to previous graph from " << firstId << " to " << lastId << " (" << graph.getVertexSet().size() << " elements)" << endl;

        clusters.push_back(graph);
        return clusters;
    }

    // attempt to subdivide even deeper any of the subgraphs
    auto newClusters1 = generateClusters(graph1);
    auto newClusters2 = generateClusters(graph2);

    clusters.insert(clusters.end(), newClusters1.begin(), newClusters1.end());
    clusters.insert(clusters.end(), newClusters2.begin(), newClusters2.end());

    return clusters;
}

TSPGraph Functionality::generateSubgraph(TSPGraph &originalGraph, int fromPosition, int toPosition) {
    TSPGraph newGraph;

    auto vertexSet = originalGraph.getVertexSet();
    unordered_set<int> validKeys;

    // add all vertexes
    for (int i = fromPosition; i < toPosition; i++) {
        newGraph.addVertex(vertexSet[i]->getInfo(), vertexSet[i]->getLongitude(), vertexSet[i]->getLatitude());
        validKeys.insert(vertexSet[i]->getInfo());
    }

    // add all valid edges
    for (auto v: newGraph.getVertexSet()) {
        auto originalVertex = originalGraph.findVertex(v->getInfo());
        for (auto e: originalVertex->getAdj()) {
            if (validKeys.count(e->getDest()->getInfo()) == 0) continue;

            // ensure you pass a vertex of the subgraph and not origin graph
            v->addEdge(newGraph.findVertex(e->getDest()->getInfo()), e->getWeight());
        }
    }

    // needed for subgraph validation
    newGraph.setOrigin(vertexSet[fromPosition]->getInfo());

    return newGraph;
}

bool Functionality::isValidSubgraph(TSPGraph &originalGraph, TSPGraph &subgraph) {
    int originId = originalGraph.getOrigin()->getInfo();
    int subgraphOrigin = subgraph.findVertex(originId) ? originId : subgraph.getVertexSet().front()->getInfo();

    bool validSubgraph = checkHamiltonianFeasibility(subgraph, subgraphOrigin);
    // compensate state clearance of previous check
    subgraph.setOrigin(subgraphOrigin);

    if (!validSubgraph) {
        int firstId = subgraph.getVertexSet().front()->getInfo();
        int lastId = subgraph.getVertexSet().back()->getInfo();
        cout << "Unable to divide into smaller subgraph from " << firstId << " to " << lastId << " (" << subgraph.getVertexSet().size() << " elements)" << endl;
    }

    return validSubgraph;
}

void Functionality::backtrackedNearestNeighbour(TSPGraph &graph, int originId, int destinationId) {
    graph.clearState();

    // backtrack system
    unsigned int nodesSize = graph.getNumVertex();
    vector<bool> visitedVector(nodesSize);
    for (int i = 0; i < nodesSize; i++)
        visitedVector[i] = false;

    // start from "origin"
    graph.setOrigin(originId);
    visitedVector[graph.getVertexIndex()[originId]] = true;

    // set initial cost to max number possible
    graph.setMinCost(std::numeric_limits<double>::max());
    vector<int> currPath;

    // recursively branch out,
    // keep track of visited nodes, distance and min cost found
    tspBacktrackingNearestNeighbour(graph, visitedVector, graph.findVertex(originId), &currPath, 1, 0,
                                    graph.findVertex(destinationId));
}

void Functionality::tspBacktrackingNearestNeighbour(TSPGraph &graph, vector<bool> &visitedVector, Vertex<int> *currNode,
                                                    vector<int> *currPath, int distance,
                                                    double cost, Vertex<int> *destinationNode) {
    if (graph.getForceStop()) return;

    if (distance == graph.getNumVertex()) {
        // If "last node" is same as source
        if (currNode == destinationNode) {
            // check if cost is less than previously found
            double currCost = cost;
            graph.setMinCost(currCost);
            graph.setMinPath(*currPath);

            // IMPORTANT:
            // script finishes on first successful route back to origin
            graph.setForceStop(true);
        }

        return;
    } else if (currNode == destinationNode) return;

    // Backtracking Part
    // Pick "nearest neighbours" before other adjacent nodes
    std::vector<Edge<int> *> sortedEdges = currNode->getAdj();
    sort(sortedEdges.begin(), sortedEdges.end(),
         [](Edge<int> *a, Edge<int> *b) {
             return a->getWeight() < b->getWeight();
         });

    for (auto adj: sortedEdges) {
        int destId = adj->getDest()->getInfo();
        int destPosition = graph.getVertexIndex()[destId];
        // if not visited
        if (!visitedVector[destPosition]) {
            // mark as visited
            visitedVector[destPosition] = true;
            currPath->push_back(destId);

            // Then recursively run this method on the adjacent node,
            // until all paths have either failed or found the source
            tspBacktrackingNearestNeighbour(graph, visitedVector, adj->getDest(), currPath, distance + 1,
                                            cost + adj->getWeight(), destinationNode);

            // mark as unvisited when backtracking
            visitedVector[destPosition] = false;
            currPath->pop_back();
        }
    }
}

bool Functionality::checkHamiltonianFeasibility(TSPGraph &graph, int originId) {
    graph.clearState();

    bool feasible = true;
    vector<int> dfsRes = graph.dfs(originId);

    // number of visitable nodes is smaller than number of total nodes
    if (dfsRes.size() != graph.getNumVertex()) feasible = false;

    // number of edges of each node must be >= 2
    for (auto v: graph.getVertexSet()) {
        if (v->getAdj().size() <= 1) feasible = false;
    }
    return feasible;
}
