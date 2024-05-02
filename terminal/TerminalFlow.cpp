#include <iostream>
#include "TerminalFlow.h"

using namespace std;

void TerminalFlow::call(WaterSupply &ws) {
    cout << endl;
    cout << "Welcome to our Routing Management System" << endl;
    cout << "---------------------------------------------" << endl;
    getReadDataMenu(ws);
    mainMenu(ws);
}

void TerminalFlow::mainMenu(WaterSupply &ws) {
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
            cout << Functionality::backtracking(ws) << endl;
            mainMenu(ws);
            break;
        case 2 :
            cout << Functionality::triangularInequality(ws) << endl;
            mainMenu(ws);
            break;
        case 3 :
            cout << Functionality::otherHeuristic(ws) << endl;
            mainMenu(ws);
            break;
        case 4 :
            cout << Functionality::notFullyConnected(ws) << endl;
            mainMenu(ws);
            break;
        case 5 :
            exit(0);
        default:
            cout << "Invalid choice. Please try again." << endl;
            mainMenu(ws);
    }
}

/*

void TerminalFlow::printPipeLoad(WaterSupply& ws) {
    vector<double> diffVector;
    cout << "Current Network / Pipe Load:" << endl;
    for(auto node : ws.getNodeSet()){
        if(node->getInfo() == "master_source" || node->getInfo() == "master_sink")
            continue;

        if(node->getAdj().empty())
            continue;

        if(node->getAdj().size() == 1 && node->getAdj()[0]->getDest()->getInfo() == "master_sink")
            continue;

        cout << "From node " << node->getInfo() << ":" << endl;
        for(auto pipe : node->getAdj()){
            if(pipe->getDest()->getInfo() == "master_sink")
                continue;

            diffVector.push_back(pipe->getWeight()-pipe->getFlow());
            cout << "- to " << pipe->getDest()->getInfo() << " (" << to_string(int(pipe->getFlow())) << "/" << to_string(int(pipe->getWeight())) <<")" << endl;
        }
    }
    cout << endl;
    vector<double> balanceStats = ws.getNetworkBalanceStats();
    cout << "Average flow: " << balanceStats[0] << endl;
    cout << "Std Dev of flow: " << balanceStats[1] << endl;
    cout << "Total flow: " << balanceStats[2] << endl;
    cout << endl;
}

string TerminalFlow::getValidCityCode(WaterSupply &ws) {
    string cityCode;
    while (true) {
        cout << "Insert a city code: (eg: C_1)" << endl;
        getline(cin, cityCode);
        if (ws.findNode(cityCode) != nullptr)
            break;

        cout << "You wrote an invalid input or city not found." << endl;
        cout << "Example of a valid input: C_2" << endl;
        cout << endl;
    }

    return cityCode;
}

string TerminalFlow::getValidReservoirCode(WaterSupply &ws) {
    string reservoirCode;
    while (true) {
        cout << "Insert a reservoir code: (eg: R_1)" << endl;
        getline(cin, reservoirCode);
        if (ws.findNode(reservoirCode) != nullptr)
            break;

        cout << "You wrote an invalid input or reservoir not found." << endl;
        cout << "Example of a valid input: R_2" << endl;
        cout << endl;
    }

    return reservoirCode;
}

void TerminalFlow::printVector(const vector<string> &resultVector) {
    for (const string &s: resultVector) {
        cout << s << endl;
    }
}

*/

void TerminalFlow::getReadDataMenu(WaterSupply &ws) {
    cout << "Please choose a graph type of the following:" << endl;
    cout << "1. Toy graphs" << endl;
    cout << "2. Extra fully connected graphs" << endl;
    cout << "3. Real world graphs" << endl;
    int selected, selected2;
    cin >> selected;
    string edges_file_path = "";
    switch (selected) {
        case 1 :
            cout << "Please choose a toy graph of the following:" << endl;
            cout << "1. Shipping" << endl;
            cout << "2. Stadiums" << endl;
            cout << "3. Tourism" << endl;
            cin >> selected2;
            switch (selected2) {
                case 1 :
                    edges_file_path = "../files/Toy-Graphs/shipping.csv";
                    break;
                case 2 :
                    edges_file_path = "../files/Toy-Graphs/stadiums.csv";
                    break;
                case 3 :
                    edges_file_path = "../files/Toy-Graphs/tourism.csv";
                    break;
                default:
                    cout << "Invalid choice. Please try again." << endl;
                    getReadDataMenu(ws);
            }
        case 2 :
            edges_file_path = "../files/Extra-Fully-Connected-Graphs/edges_25.csv";
            break;
        case 3 :
            cout << "Please choose a real world graph of the following:" << endl;
            cout << "1. Graph 1" << endl;
            cout << "2. Graph 2" << endl;
            cout << "3. Graph 3" << endl;
            cin >> selected2;
            switch (selected2) {
                case 1 :
                    edges_file_path = "../files/Real-world-Graphs/graph1/edges.csv";
                    break;
                case 2 :
                    edges_file_path = "../files/Real-world-Graphs/graph2/edges.csv";;
                    break;
                case 3 :
                    edges_file_path = "../files/Real-world-Graphs/graph3/edges.csv";
                    break;
                default:
                    cout << "Invalid choice. Please try again." << endl;
                    getReadDataMenu(ws);
            }
        default:
            cout << "Invalid choice. Please try again." << endl;
            getReadDataMenu(ws);
    }
    FileReader::addPipes(edges_file_path, ws);
    cout << "Files have been read correctly" << endl;
}



