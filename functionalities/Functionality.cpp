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
    cout << "Min Path: ";
    for (auto element : graph.getMinPath()) {
        cout << element << " ";
    }
    cout << endl << endl;
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
                if(graph.getMinCost() > currCost){
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
        double currCost = cost + adj->getWeight();
        if(currCost >= graph.getMinCost()){
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

/*
void triangularInequality(WaterSupply& graph){

}

void otherHeuristic(WaterSupply& graph){

}

void notFullyConnected(WaterSupply& graph){
    
}*/