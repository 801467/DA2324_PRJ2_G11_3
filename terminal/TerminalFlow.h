#ifndef DA2324_PRJ2_G11_3_TERMINALFLOW_H
#define DA2324_PRJ2_G11_3_TERMINALFLOW_H

#include <iostream>
#include "../classes/TSPGraph.h"
#include "../utils/FileReader.h"
#include "../functionalities/Functionality.h"

class TerminalFlow {
public:
    static void call(TSPGraph& graph);
    static void mainMenu(TSPGraph& graph);
    static void chooseGraph(TSPGraph& graph);
private:
    static void loadGraphs();
    static void chooseToyGraph(TSPGraph& graph);
    static void chooseFullyConnectedGraph(TSPGraph& graph);
    static void chooseRealWorldGraph(TSPGraph& graph);

};


#endif //DA2324_PRJ2_G11_3_TERMINALFLOW_H
