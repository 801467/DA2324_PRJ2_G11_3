#include <unordered_set>
#include "Functionality.h"

/**
 * @brief Runs the backtracking algorithm on a TSPGraph.
 * 
 * This function runs a backtracking algorithm to solve the Traveling Salesman Problem (TSP) on a given graph.
 * It starts from the "origin" node, and recursively explores all possible paths, keeping track of the minimum cost found.
 * 
 * @param graph The TSPGraph to run the algorithm on.
 */
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

/**
 * @brief Recursive helper function for the backtracking algorithm.
 * 
 * This function is called by the backtracking function to recursively explore all possible paths in the TSPGraph.
 * It updates the visitedVector, currPath, distance, and cost as it explores the graph.
 * 
 * @param graph The TSPGraph to run the algorithm on.
 * @param visitedVector A vector keeping track of which nodes have been visited.
 * @param currNode The current node being visited.
 * @param currPath A vector representing the current path being explored.
 * @param distance The number of nodes visited so far.
 * @param cost The total cost of the current path.
 */
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


/**
 * @brief Runs the Triangular Approximation Heuristic on a TSPGraph.
 * 
 * This function runs the Triangular Approximation Heuristic on the given graph.
 * The heuristic works by computing the Minimum Spanning Tree (MST) of the graph, performing a pre-order walk of the MST,
 * and then performing the travel of the TSP tour based on the pre-order walk. 
 * If it cannot find an edge connecting two consecutive nodes of the path,
 * it calculates the distance between them using the Haversine formula.
 * The cost of the tour is then calculated and displayed.
 * 
 * @param graph The TSPGraph to run the heuristic on.
 */
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


/**
 * @brief Computes the minimum spanning tree (MST) from the origin using Prim's algorithm.
 * 
 * This function uses Prim's algorithm to compute the minimum spanning tree of the given graph.
 * The MST is a subset of the edges of the graph that connects all the vertices together,
 *  without any cycles and with the minimum possible total edge weight.
 * 
 * @param graph The TSPGraph to compute the MST on.
 */
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

/**
 * @brief Performs a pre-order walk of the minimum spanning tree (MST).
 * 
 * This function performs a pre-order walk (root, left, right) of the MST of the given graph.
 * The result is a vector of vertices representing the order in which they are visited during the walk.
 * 
 * @param graph The TSPGraph to perform the pre-order walk on.
 * @return A vector of vertices representing the order of the pre-order walk.
 */
vector<Vertex<int> *> Functionality::preOrderWalk(TSPGraph &graph) {
    vector<Vertex<int> *> orderedPath;
    Vertex<int> *origin = graph.getOrigin();

    for (auto v: graph.getVertexSet()) {
        v->setVisited(false);
    }

    tspDfsVisit(origin, orderedPath);


    return orderedPath;
}


/**
 * @brief Helper function for the pre-order walk of the MST.
 * 
 * This function is called by the preOrderWalk function to perform the pre-order walk of the MST.
 * It visits the given vertex and recursively checks its adjacent vertices to see if they have been visited
 * and if the edge connecting them is the path to follow. If so, it visits the adjacent vertex.
 * 
 * @param v The vertex to visit.
 * @param orderedPath A vector of vertices representing the order of the pre-order walk.
 */
void Functionality::tspDfsVisit(Vertex<int> *v, vector<Vertex<int> *> &orderedPath) {
    v->setVisited(true);
    orderedPath.push_back(v);
    for (auto e: v->getAdj()) {
        auto w = e->getDest();
        if (!w->isVisited() && (w->getPath() == e)) tspDfsVisit(w, orderedPath);
    }
}

/**
 * @brief Performs the travel of the TSP tour.
 * 
 * This function performs the travel of the TSP tour based on the given path.
 * If it cannot find an edge connecting two nodes, it calculates the distance between them using the Haversine formula,
 * and adds it to the cost of the tour. If there is an edge connecting the nodes, the cost of the edge is added to the cost of the tour.
 * The cost of the tour is then returned.
 * 
 * @param path A vector of vertices representing the order of the TSP tour.
 * @return The cost of the TSP tour.
 */
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


/**
 * @brief Runs the Nearest Neighbour Heuristic on a TSPGraph.
 * 
 * This function runs the Nearest Neighbour Heuristic to solve the Traveling Salesman Problem (TSP) on a given graph.
 * First it uses the Adjacency Matrix to store the distance between each pair of nodes.
 * Then starts from the "origin" node, and at each step, it moves to the unvisited node that is closest to the current node.
 * The cost of the tour is calculated as the sum of the distances of each step.
 * 
 * @param graph The TSPGraph to run the algorithm on.
 */
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

/**
 * @brief Runs the Clustered, Backtracked with Nearest Neighbour Heuristic on a TSPGraph.
 * 
 * This function runs the Clustered, Backtracked with Nearest Neighbour Heuristic on the given graph.
 * It generates clusters of nodes with less than 150 nodes in each, and then runs the backtracked nearest neighbour heuristic on each cluster.
 * The cost of the tour is then calculated and displayed.
 * 
 * @param graph The TSPGraph to run the heuristic on.
 * @param originId The ID of the origin node.
 */
void Functionality::runClusteredGraphs(TSPGraph &graph, int originId) {
    cout << "Running Clustered, Backtracked with Nearest Neighbour Heuristic..." << endl;
    cout << endl;

    graph.clearState();
    graph.setOrigin(originId);

    auto clusters = generateClusters(graph);
    cout << "Generated a total of " << clusters.size() << " clusters. Less than 150 nodes in each." << endl << endl;

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

/**
 * @brief Generates clusters of nodes with less than 150 nodes in each.
 * 
 * This function generates clusters of nodes with less than 150 nodes in each from the given graph.
 * It splits the graph into two subgraphs, and then attempts to subdivide even deeper any of the subgraphs.
 * The function then returns a vector of TSPGraphs representing the clusters.
 * 
 * @param graph The TSPGraph to generate the clusters from.
 * @return A vector of TSPGraphs representing the clusters.
 */
vector<TSPGraph> Functionality::generateClusters(TSPGraph &graph) {
    vector<TSPGraph> clusters;

    // Avoid unnecessary smaller clustering
    if (graph.getVertexSet().size() < 150) {
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
        cout << "Rolled back to previous graph from " << firstId << " to " << lastId << endl;

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

/**
 * @brief Generates a subgraph from the original graph.
 * 
 * This function generates a subgraph from the original graph, starting from the "fromPosition" node and ending at the "toPosition" node.
 * It adds all the vertexes from the original graph to the new graph, and then adds all the valid edges from the original graph to the new graph.
 * The function then returns the new graph.
 * 
 * @param originalGraph The original graph to generate the subgraph from.
 * @param fromPosition The position of the starting node.
 * @param toPosition The position of the ending node.
 * @return The new subgraph.
 */
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

/**
 * @brief Checks if a subgraph is valid.
 * 
 * This function checks if a subgraph is valid by checking if it is Hamiltonian feasible.
 * If the subgraph is not Hamiltonian feasible, the function displays an error message.
 * 
 * @param originalGraph The original graph.
 * @param subgraph The subgraph to check.
 * @return True if the subgraph is valid, false otherwise.
 */
bool Functionality::isValidSubgraph(TSPGraph &originalGraph, TSPGraph &subgraph) {
    int originId = originalGraph.getOrigin()->getInfo();
    int subgraphOrigin = subgraph.findVertex(originId) ? originId : subgraph.getVertexSet().front()->getInfo();

    bool validSubgraph = checkHamiltonianFeasibility(subgraph, subgraphOrigin);
    // compensate state clearance of previous check
    subgraph.setOrigin(subgraphOrigin);

    if (!validSubgraph) {
        int firstId = subgraph.getVertexSet().front()->getInfo();
        int lastId = subgraph.getVertexSet().back()->getInfo();
        cout << "Unable to divide into smaller subgraph from " << firstId << " to " << lastId << endl;
    }

    return validSubgraph;
}

/**
 * @brief Runs the backtracked nearest neighbour heuristic on a TSPGraph.
 * 
 * This function runs the backtracked nearest neighbour heuristic on the given graph, starting from the "originId" node and ending at the "destinationId" node.
 * It clears the state of the graph, and then runs the backtracked nearest neighbour heuristic on the graph.
 * 
 * @param graph The TSPGraph to run the heuristic on.
 * @param originId The ID of the origin node.
 * @param destinationId The ID of the destination node.
 */
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

/**
 * @brief Recursive helper function for the backtracked nearest neighbour heuristic.
 * 
 * This function is called by the backtrackedNearestNeighbour function to recursively explore all possible paths in the TSPGraph.
 * It updates the visitedVector, currPath, distance, and cost as it explores the graph.
 * 
 * @param graph The TSPGraph to run the heuristic on.
 * @param visitedVector A vector keeping track of which nodes have been visited.
 * @param currNode The current node being visited.
 * @param currPath A vector representing the current path being explored.
 * @param distance The number of nodes visited so far.
 * @param cost The total cost of the current path.
 * @param destinationNode The destination node.
 */
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

/**
 * @brief Checks if a graph is Hamiltonian feasible.
 * 
 * This function checks if a graph is Hamiltonian feasible by checking if the number of visitable nodes is smaller than the number of total nodes,
 * and if the number of edges of each node is greater than or equal to 2.
 * If the graph is not Hamiltonian feasible, the function returns false.
 * 
 * @param graph The TSPGraph to check.
 * @param originId The ID of the origin node.
 * @return True if the graph is Hamiltonian feasible, false otherwise.
 */
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
