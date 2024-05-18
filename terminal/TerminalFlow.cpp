#include <iostream>
#include "TerminalFlow.h"

using namespace std;

void TerminalFlow::welcome(TSPGraph &graph) {
    cout << endl;
    cout << "Welcome to our Routing Management System" << endl;
    cout << "---------------------------------------------" << endl;
    cout << endl;
    call(graph);
}

void TerminalFlow::call(TSPGraph &graph) {
    chooseGraph(graph);
    runFunctionality(graph);
}

void TerminalFlow::chooseGraph(TSPGraph &graph) {
    cout << "Please choose a graph type of the following:" << endl;
    cout << "1. Toy graphs" << endl;
    cout << "2. Extra fully connected graphs" << endl;
    cout << "3. Real world graphs" << endl;
    int selected;
    cin >> selected;
    switch (selected) {
        case 1 :
            chooseToyGraph(graph);
            break;
        case 2 :
            chooseFullyConnectedGraph(graph);
            break;
        case 3 :
            chooseRealWorldGraph(graph);
            break;
        default:
            cout << "Invalid choice. Please try again." << endl;
            chooseGraph(graph);
    }
    cout << "Graph has been loaded correctly" << endl;
}

void TerminalFlow::runFunctionality(TSPGraph &graph) {
    cout << endl;
    cout << "Please choose a functionality of the following:" << endl;
    cout << "---------------------------------------------" << endl;
    cout << "1. Backtracking algorithm for TSP" << endl;
    cout << "2. Approximation algorithm for TSP with triangular inequality" << endl;
    cout << "3. Approximation algorithm for TSP with nearest neighbour heuristic" << endl;
    cout << "4. Approximation algorithm for TSP when not fully connected" << endl;
    cout << "5. Choose a different graph." << endl;
    cout << "6. Exit" << endl;

    int selected;
    cin >> selected;
    cin.ignore(9999999, '\n'); // always clean cin (to avoid issue with getline)
    cout << endl;

    chrono::time_point<chrono::system_clock> start;

    switch (selected) {
        case 1 :
            start = chrono::system_clock::now();
            Functionality::backtracking(graph);
            printTimeLapsed(start);
            runFunctionality(graph);
            break;
        case 2 :
            start = chrono::system_clock::now();
            Functionality::triangularInequality(graph);
            printTimeLapsed(start);
            runFunctionality(graph);
            break;
        case 3 :
            start = chrono::system_clock::now();
            Functionality::nearestNeighbour(graph);
            printTimeLapsed(start);
            runFunctionality(graph);
            break;
        case 4 :
            start = chrono::system_clock::now();
            Functionality::backtrackedNearestNeighbour(graph, chooseOriginNode(graph));
            printTimeLapsed(start);
            runFunctionality(graph);
            break;
        case 5 :
            resetGraph(graph);
            call(graph);
            break;
        case 6 :
            exit(0);
        default:
            cout << "Invalid choice. Please try again." << endl;
            runFunctionality(graph);
    }
}

// private

void TerminalFlow::chooseToyGraph(TSPGraph &graph) {
    cout << "Please choose a toy graph of the following:" << endl;
    cout << "1. Shipping" << endl;
    cout << "2. Stadiums" << endl;
    cout << "3. Tourism" << endl;

    int selected;
    cin >> selected;
    switch (selected) {
        case 1 :
            loadGraph("TGraphShipping", graph);
            break;
        case 2 :
            loadGraph("TGraphStadiums", graph);
            break;
        case 3 :
            loadGraph("TGraphTourism", graph);
            break;
        default:
            cout << "Invalid choice. Please try again." << endl;
            chooseGraph(graph);
    }
}

void TerminalFlow::chooseFullyConnectedGraph(TSPGraph &graph) {
    cout << "Please choose a fully connected graph of the following:" << endl;
    cout << "1. 25 Node Graph " << endl;
    cout << "2. 50 Node Graph " << endl;
    cout << "3. 75 Node Graph " << endl;
    cout << "4. 100 Node Graph " << endl;
    cout << "5. 200 Node Graph " << endl;
    cout << "6. 300 Node Graph " << endl;
    cout << "7. 400 Node Graph " << endl;
    cout << "8. 500 Node Graph " << endl;
    cout << "9. 600 Node Graph " << endl;
    cout << "10. 700 Node Graph " << endl;
    cout << "11. 800 Node Graph " << endl;
    cout << "12. 900 Node Graph " << endl;

    int selected;
    cin >> selected;
    switch (selected) {
        case 1 :
            loadGraph("EFCGraph25", graph);
            break;
        case 2 :
            loadGraph("EFCGraph50", graph);
            break;
        case 3 :
            loadGraph("EFCGraph75", graph);
            break;
        case 4 :
            loadGraph("EFCGraph100", graph);
            break;
        case 5 :
            loadGraph("EFCGraph200", graph);
            break;
        case 6 :
            loadGraph("EFCGraph300", graph);
            break;
        case 7 :
            loadGraph("EFCGraph400", graph);
            break;
        case 8 :
            loadGraph("EFCGraph500", graph);
            break;
        case 9 :
            loadGraph("EFCGraph600", graph);
            break;
        case 10 :
            loadGraph("EFCGraph700", graph);
            break;
        case 11 :
            loadGraph("EFCGraph800", graph);
            break;
        case 12 :
            loadGraph("EFCGraph900", graph);
            break;
        default:
            cout << "Invalid choice. Please try again." << endl;
            chooseGraph(graph);
    }
}

void TerminalFlow::chooseRealWorldGraph(TSPGraph &graph) {
    cout << "Please choose a real world graph of the following:" << endl;
    cout << "1. Graph 1" << endl;
    cout << "2. Graph 2" << endl;
    cout << "3. Graph 3" << endl;

    int selected;
    cin >> selected;
    switch (selected) {
        case 1 :
            loadGraph("RWGraph1", graph);
            break;
        case 2 :
            loadGraph("RWGraph2", graph);
            break;
        case 3 :
            loadGraph("RWGraph3", graph);
            break;
        default:
            cout << "Invalid choice. Please try again." << endl;
            chooseGraph(graph);
    }
}

void TerminalFlow::loadGraph(const std::string &chosenGraph, TSPGraph &graph) {
    if (chosenGraph == "TGraphShipping")
        FileReader::loadGraph("../files/Toy-Graphs/shipping.csv", graph);
    else if (chosenGraph == "TGraphStadiums")
        FileReader::loadGraph("../files/Toy-Graphs/stadiums.csv", graph);
    else if (chosenGraph == "TGraphTourism")
        FileReader::loadGraph("../files/Toy-Graphs/tourism.csv", graph);
    else if (chosenGraph == "RWGraph1") {
        cout << "Takes ~10 sec to load. Please wait..." << endl;
        FileReader::loadGraph("../files/Real-world-Graphs/graph1/edges.csv",
                              "../files/Real-world-Graphs/graph1/nodes.csv",
                              graph);
    } else if (chosenGraph == "RWGraph2") {
        cout << "Takes ~1 min to load. Please wait..." << endl;
        FileReader::loadGraph("../files/Real-world-Graphs/graph2/edges.csv",
                              "../files/Real-world-Graphs/graph2/nodes.csv", graph);
    } else if (chosenGraph == "RWGraph3") {
        cout << "Takes ~2 min to load. Please wait..." << endl;
        FileReader::loadGraph("../files/Real-world-Graphs/graph3/edges.csv",
                              "../files/Real-world-Graphs/graph3/nodes.csv", graph);
    } else if (chosenGraph == "EFCGraph25")
        FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_25.csv",
                              "../files/Extra-Fully-Connected-Graphs/nodes.csv", graph);
    else if (chosenGraph == "EFCGraph50")
        FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_50.csv",
                              "../files/Extra-Fully-Connected-Graphs/nodes.csv", graph);
    else if (chosenGraph == "EFCGraph75")
        FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_75.csv",
                              "../files/Extra-Fully-Connected-Graphs/nodes.csv", graph);
    else if (chosenGraph == "EFCGraph100")
        FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_100.csv",
                              "../files/Extra-Fully-Connected-Graphs/nodes.csv", graph);
    else if (chosenGraph == "EFCGraph200")
        FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_200.csv",
                              "../files/Extra-Fully-Connected-Graphs/nodes.csv", graph);
    else if (chosenGraph == "EFCGraph300")
        FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_300.csv",
                              "../files/Extra-Fully-Connected-Graphs/nodes.csv", graph);
    else if (chosenGraph == "EFCGraph400")
        FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_400.csv",
                              "../files/Extra-Fully-Connected-Graphs/nodes.csv", graph);
    else if (chosenGraph == "EFCGraph500")
        FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_500.csv",
                              "../files/Extra-Fully-Connected-Graphs/nodes.csv", graph);
    else if (chosenGraph == "EFCGraph600")
        FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_600.csv",
                              "../files/Extra-Fully-Connected-Graphs/nodes.csv", graph);
    else if (chosenGraph == "EFCGraph700")
        FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_700.csv",
                              "../files/Extra-Fully-Connected-Graphs/nodes.csv", graph);
    else if (chosenGraph == "EFCGraph800")
        FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_800.csv",
                              "../files/Extra-Fully-Connected-Graphs/nodes.csv", graph);
    else if (chosenGraph == "EFCGraph900")
        FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_900.csv",
                              "../files/Extra-Fully-Connected-Graphs/nodes.csv", graph);
    else
        cout << "Unknown Graph. Check for a possible typo." << endl;
}


void TerminalFlow::printTimeLapsed(const chrono::time_point<chrono::system_clock> &start) {
    auto end = chrono::system_clock::now();

    chrono::duration<double> elapsed_seconds = end - start;

    std::time_t start_time = chrono::system_clock::to_time_t(start);
    std::time_t end_time = chrono::system_clock::to_time_t(end);

    cout << "Started algorithm at " << ctime(&start_time)
         << "Finished algorithm at " << ctime(&end_time)
         << "Elapsed time: " << elapsed_seconds.count() << "s"
         << endl;
}

void TerminalFlow::resetGraph(TSPGraph &graph) {
    TSPGraph emptyGraph;
    graph = emptyGraph;
}

int TerminalFlow::chooseOriginNode(TSPGraph &graph) {
    int selected = -1;
    while (selected == -1) {
        cout << "Please choose an Origin Vertex ID, from 0 to " << graph.getNumVertex() - 1 << ":" << endl;
        cin >> selected;
        if (selected >= graph.getNumVertex() || selected < 0)
            selected = -1;
    }
    return selected;
}
