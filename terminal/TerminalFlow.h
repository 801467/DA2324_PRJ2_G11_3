#ifndef DA2324_PRJ2_G11_3_TERMINALFLOW_H
#define DA2324_PRJ2_G11_3_TERMINALFLOW_H

#include <iostream>
#include "../classes/WaterSupply.h"
#include "../utils/FileReader.h"
#include "../functionalities/Functionality.h"

class TerminalFlow {
public:
    static void call(WaterSupply& ws);
    static void getReadDataMenu(WaterSupply& ws);
    static void mainMenu(WaterSupply& ws);
private:
    static void printVector(const vector<string>& resultVector);
    static string getValidCityCode(WaterSupply& ws);
    static string getValidReservoirCode(WaterSupply& ws);
    static void printPipeLoad(WaterSupply& ws);
};


#endif //DA2324_PRJ2_G11_3_TERMINALFLOW_H
