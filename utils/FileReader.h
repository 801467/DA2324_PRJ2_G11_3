#ifndef FILEREADER_H
#define FILEREADER_H

#include "../classes/TSPGraph.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>


class FileReader {
public:
    static void loadNodes(const std::string &filename);
    static void loadGraph(const std::string &filename, TSPGraph &graph);
    static void loadGraph(const std::string &edgesfilename, const std::string &nodesfilename,TSPGraph &graph);
};

#endif //FILEREADER_H