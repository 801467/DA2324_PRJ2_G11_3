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

    graph.clearState();
    double cost;

    // select root vertex
    graph.setOrigin(0);

    //compute MST from origin
    prim(graph);

    //perform a pre-order walk of the MST
    vector<Vertex<int>*> path = preOrderWalk(graph);

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

vector<Vertex<int> *> Functionality::preOrderWalk(TSPGraph &graph){
    vector<Vertex<int>*> orderedPath;
    Vertex<int>* origin = graph.getOrigin();

    for (auto v : graph.getVertexSet()){
        v->setVisited(false);
    }

    tspDfsVisit(origin,orderedPath);


    return orderedPath;
}

void Functionality::tspDfsVisit(Vertex<int>* v, vector<Vertex<int>*> &orderedPath) {
    v->setVisited(true);
    orderedPath.push_back(v);
    for (auto e: v->getAdj()){
        auto w = e->getDest();
        if (!w->isVisited() && (w->getPath() == e)) tspDfsVisit(w,orderedPath);
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
        if(MIN == INF)
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

void Functionality::backtrackedNearestNeighbour(TSPGraph &graph, int origin) {
    cout << "Running Backtracked Nearest Neighbour Heuristic..." << endl;
    cout << endl;

    graph.clearState();

    // TODO
}

void Functionality::tspBacktrackingNearestNeighbour(TSPGraph &graph, vector<bool> &visitedVector, Vertex<int> *currNode,
                                                    vector<int> *currPath, int distance,
                                                    double cost) {
    // TODO
}

void Functionality::checkHamiltonianFeasibility(TSPGraph &graph, int originId) {
    graph.clearState();

    bool feasible = true;
    vector<int> dfsRes = graph.dfs(originId);

    // number of visitable nodes is smaller than number of total nodes
    if (dfsRes.size() != graph.getNumVertex()) feasible = false;

    // number of edges of each node must be >= 2
    for(auto v : graph.getVertexSet()){
        if (v->getAdj().size() <= 1) feasible = false;
    }
    feasible ? cout << "True" << endl : cout << "False" << endl;

}
