#include <iostream>
#include "TerminalFlow.h"


static TSPGraph TGraphTourism, TGraphStadiums, TGraphShipping, RWGraph1, RWGraph2, RWGraph3, EFCGraph25, EFCGraph50, EFCGraph75, EFCGraph100, EFCGraph200, EFCGraph300,
        EFCGraph400, EFCGraph500, EFCGraph600, EFCGraph700, EFCGraph800, EFCGraph900;


using namespace std;

void TerminalFlow::call(TSPGraph &graph) {
    cout << endl;
    cout << "Welcome to our Routing Management System" << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Parsing files... Please wait.";
    loadGraphs();
    cout << '\r';   // this should delete the loading line but doesn't...
    chooseGraph(graph);
    mainMenu(graph);
}

void TerminalFlow::chooseGraph(TSPGraph &graph) {
    cout << "Please choose a graph type of the following:" << endl;
    cout << "1. Toy graphs" << endl;
    cout << "2. Extra fully connected graphs" << endl;
    cout << "3. Real world graphs" << endl;
    int selected, selected2;
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
    cout << "Graph has been chosen correctly" << endl;
}

void TerminalFlow::chooseToyGraph(TSPGraph &graph) {
    cout << "Please choose a toy graph of the following:" << endl;
    cout << "1. Shipping" << endl;
    cout << "2. Stadiums" << endl;
    cout << "3. Tourism" << endl;

    int selected;
    cin >> selected;
    switch (selected) {
        case 1 :
            graph = TGraphShipping;
            break;
        case 2 :
            graph = TGraphStadiums;
            break;
        case 3 :
            graph = TGraphTourism;
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
            graph = EFCGraph25;
            break;
        case 2 :
            graph = EFCGraph50;
            break;
        case 3 :
            graph = EFCGraph75;
            break;
        case 4 :
            graph = EFCGraph100;
            break;
        case 5 :
            graph = EFCGraph200;
            break;
        case 6 :
            graph = EFCGraph300;
            break;
        case 7 :
            graph = EFCGraph400;
            break;
        case 8 :
            graph = EFCGraph500;
            break;
        case 9 :
            graph = EFCGraph600;
            break;
        case 10 :
            graph = EFCGraph700;
            break;
        case 11 :
            graph = EFCGraph800;
            break;
        case 12 :
            graph = EFCGraph900;
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
            graph = RWGraph1;
            break;
        case 2 :
            graph = RWGraph2;
            break;
        case 3 :
            graph = RWGraph3;
            break;
        default:
            cout << "Invalid choice. Please try again." << endl;
            chooseGraph(graph);
    }
}

void TerminalFlow::mainMenu(TSPGraph &graph) {
    cout << endl;
    cout << "Please choose a functionality of the following:" << endl;
    cout << "---------------------------------------------" << endl;
    cout << "1. Backtracking algorithm for TSP" << endl;
    cout << "2. Approximation algorithm for TSP with triangular inequality" << endl;
    cout << "3. Approximation algorithm for TSP with another heuristic" << endl;
    cout << "4. Approximation algorithm for TSP when not fully connected" << endl;
    cout << "5. Exit" << endl;

    int selected;
    cin >> selected;
    cin.ignore(9999999, '\n'); // always clean cin (to avoid issue with getline)
    cout << endl;

    switch (selected) {
        case 1 :
            cout << "Functionality::backtracking(graph)" << endl;
            mainMenu(graph);
            break;
        case 2 :
            cout << "Functionality::triangularInequality(graph)" << endl;
            mainMenu(graph);
            break;
        case 3 :
            cout << "Functionality::otherHeuristic(graph)" << endl;
            mainMenu(graph);
            break;
        case 4 :
            cout << "Functionality::notFullyConnected(graph)" << endl;
            mainMenu(graph);
            break;
        case 5 :
            exit(0);
        default:
            cout << "Invalid choice. Please try again." << endl;
            mainMenu(graph);
    }
}

void TerminalFlow::loadGraphs() {
    FileReader::loadGraph("../files/Toy-Graphs/shipping.csv", TGraphShipping);
    FileReader::loadGraph("../files/Toy-Graphs/stadiums.csv", TGraphStadiums);
    FileReader::loadGraph("../files/Toy-Graphs/tourism.csv", TGraphTourism);
    FileReader::loadGraph("../files/Real-world-Graphs/graph1/edges.csv", "../files/Real-world-Graphs/graph1/nodes.csv",
                          RWGraph1);
    //FileReader::loadGraph("../files/Real-world-Graphs/graph2/edges.csv","../files/Real-world-Graphs/graph2/nodes.csv",RWGraph2); takes 1~min to load
    //FileReader::loadGraph("../files/Real-world-Graphs/graph3/edges.csv","../files/Real-world-Graphs/graph3/nodes.csv",RWGraph3); takes 10+min to load
    FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_25.csv",
                          "../files/Extra-Fully-Connected-Graphs/nodes.csv", EFCGraph25);
    FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_50.csv",
                          "../files/Extra-Fully-Connected-Graphs/nodes.csv", EFCGraph50);
    FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_75.csv",
                          "../files/Extra-Fully-Connected-Graphs/nodes.csv", EFCGraph75);
    FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_100.csv",
                          "../files/Extra-Fully-Connected-Graphs/nodes.csv", EFCGraph100);
    FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_200.csv",
                          "../files/Extra-Fully-Connected-Graphs/nodes.csv", EFCGraph200);
    FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_300.csv",
                          "../files/Extra-Fully-Connected-Graphs/nodes.csv", EFCGraph300);
    FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_400.csv",
                          "../files/Extra-Fully-Connected-Graphs/nodes.csv", EFCGraph400);
    FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_500.csv",
                          "../files/Extra-Fully-Connected-Graphs/nodes.csv", EFCGraph500);
    FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_600.csv",
                          "../files/Extra-Fully-Connected-Graphs/nodes.csv", EFCGraph600);
    FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_700.csv",
                          "../files/Extra-Fully-Connected-Graphs/nodes.csv", EFCGraph700);
    FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_800.csv",
                          "../files/Extra-Fully-Connected-Graphs/nodes.csv", EFCGraph800);
    FileReader::loadGraph("../files/Extra-Fully-Connected-Graphs/edges_900.csv",
                          "../files/Extra-Fully-Connected-Graphs/nodes.csv", EFCGraph900);
}




