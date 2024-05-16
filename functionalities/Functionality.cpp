#include "Functionality.h"

void Functionality::backtracking(TSPGraph &graph) {
    cout << "Running Backtracking..." << endl;
    cout << endl;

    unsigned int nodesSize = graph.getNumVertex();
    vector<bool> visitedVector(nodesSize);
    for (int i = 0; i < nodesSize; i++)
        visitedVector[i] = false;

    // start from "origin"
    graph.setOrigin(0);
    visitedVector[0] = true;

    // set initial cost to max number possible
    graph.setMinCost(std::numeric_limits<double>::max());

    // reset min path
    vector<int> currPath;
    graph.setMinPath(currPath);

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
        // if not visited
        if (!visitedVector[destId]) {
            // mark as visited
            visitedVector[destId] = true;
            currPath->push_back(destId);

            // Then recursively run this method on the adjacent node,
            // until all paths have either failed or found the source
            tspBacktracking(graph, visitedVector, adj->getDest(), currPath, distance + 1,
                            cost + adj->getWeight());

            // mark as unvisited when backtracking
            visitedVector[destId] = false;
            currPath->pop_back();
        }
    }
}


void Functionality::triangularInequality(TSPGraph &graph) {
    cout << "Running Triangular Approximation Heuristic..." << endl;
    cout << endl;

    double cost;

    // select root vertex
    graph.setOrigin(0);

    //compute MST from origin
    vector<Vertex<int> *> path = prim(graph);

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

vector<Vertex<int> *> Functionality::prim(TSPGraph &graph) {
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
    return path;
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

    double cost;
    vector<Vertex<int> *> visited;
    vector<double> minimum_distance_traveled;

    auto nodes = graph.getVertexSet();
    auto distance = graph.getAdjacencyMatrix();
    auto neighbor = graph.getOrigin();
    // HERE
    auto start_node_index = std::distance(nodes.begin(), find(nodes.begin(), nodes.end(), neighbor));

    std::vector<Vertex<int> *>::size_type no_nodes = nodes.size();
    int noN = 0;
    while (noN < no_nodes && find(visited.begin(), visited.end(), neighbor) == visited.end()) {

        visited.push_back(neighbor);
        auto neighbor_index = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), neighbor));
        int noNeighbour = 0;
        double MIN = INF;

        while (noNeighbour < distance[neighbor_index].size()) {

            if (std::find(visited.begin(), visited.end(), nodes[noNeighbour]) ==
                visited.end()) { //look for unvisitied nodes
                if (MIN == INF) {
                    MIN = distance[neighbor_index][noNeighbour];
                    neighbor = nodes[noNeighbour];
                } else {
                    double min_distance = min(distance[neighbor_index][noNeighbour], MIN);
                    if (distance[neighbor_index][noNeighbour] < MIN) {
                        MIN = min_distance;
                        neighbor = nodes[noNeighbour];
                    }
                }
            }
            noNeighbour += 1;
        }
        minimum_distance_traveled.push_back(MIN);
        noN += 1;
    }
    auto last_node_index = std::distance(nodes.begin(), find(nodes.begin(), nodes.end(), visited.back()));
    minimum_distance_traveled.back() = distance[last_node_index][start_node_index];

    cost = std::accumulate(minimum_distance_traveled.begin(), minimum_distance_traveled.end(), 0.0);
    graph.setMinCost(cost);

    cout << "Min Cost: " << graph.getMinCost() << endl;
    cout << "Path: ";
    unsigned int i = 0;
    for (auto element: visited) {
        i++;
        (i % 20 == 0) ? cout << endl : cout << element->getInfo() << " ";
    }
    cout << visited.front()->getInfo() << endl << endl;
}

/*
void notFullyConnected(WaterSupply& graph){
    
}*/