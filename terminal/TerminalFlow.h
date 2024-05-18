#ifndef DA2324_PRJ2_G11_3_TERMINALFLOW_H
#define DA2324_PRJ2_G11_3_TERMINALFLOW_H

#include <iostream>
#include <chrono>
#include <ctime>
#include "../classes/TSPGraph.h"
#include "../utils/FileReader.h"
#include "../functionalities/Functionality.h"

class TerminalFlow {
public:
    static void welcome(TSPGraph &graph);
    static void call(TSPGraph &graph);
    static void chooseGraph(TSPGraph &graph);
    static void runFunctionality(TSPGraph &graph);

private:
    static void loadGraph(const std::string &chosenGraph, TSPGraph &graph);
    static void chooseToyGraph(TSPGraph &graph);
    static void chooseFullyConnectedGraph(TSPGraph &graph);
    static void chooseRealWorldGraph(TSPGraph &graph);
    static void printTimeLapsed(const chrono::time_point<chrono::system_clock> &start);
    static void resetGraph(TSPGraph &graph);
    static int chooseOriginNode(TSPGraph &graph);
    static void runNotConnectedFunctionalities(TSPGraph &graph);
};


#endif //DA2324_PRJ2_G11_3_TERMINALFLOW_H
