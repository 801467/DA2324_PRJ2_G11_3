#ifndef DA2324_PRJ2_G11_3_FUNCTIONALITY_H
#define DA2324_PRJ2_G11_3_FUNCTIONALITY_H

#include "../classes/WaterSupply.h"
#include <iostream>

class Functionality {
public:
    static void backtracking(WaterSupply& graph);
    static void triangularInequality(WaterSupply& graph);
    static void otherHeuristic(WaterSupply& graph);
    static void notFullyConnected(WaterSupply& graph);
};


#endif //DA2324_PRJ2_G11_3_FUNCTIONALITY_H
