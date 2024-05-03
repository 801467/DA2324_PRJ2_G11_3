#include "FileReader.h"

unordered_map<int,pair<double,double>> nodes;

using namespace std;

void FileReader::loadNodes(const std::string &filename) {
    nodes.clear();
    ifstream inputFile(filename);

    if(inputFile.is_open()) {
        string line;
        getline(inputFile, line);

        while (getline(inputFile, line)) {
            istringstream iss(line);
            string id, longitude, latitude;

            if (getline(iss, id, ',') &&
                getline(iss, longitude, ',') &&
                getline(iss, latitude, '\r')) {
                nodes[stoi(id)] = make_pair(stod(longitude),stod(latitude));
            }
        }
    }
    inputFile.close();
}

void FileReader::loadGraph(const std::string &filename, TSPGraph &graph) {
    ifstream inputFile(filename);

    if (inputFile.is_open()) {
        string line;
        getline(inputFile, line);

        while (getline(inputFile, line)) {
            istringstream iss(line);
            string origin, destination, distance;

            if (filename != "../files/Toy-Graphs/tourism.csv") {
                if (getline(iss, origin, ',') &&
                    getline(iss, destination, ',') &&
                    getline(iss, distance, '\r')) {
                    int originID = stoi(origin);
                    int destinationID = stoi(destination);
                    graph.Graph::addVertex(originID);
                    graph.Graph::addVertex(destinationID);
                    graph.addBidirectionalEdge(originID, destinationID, stod(distance));
                }
            } else {
                if (getline(iss, origin, ',') &&
                    getline(iss, destination, ',') &&
                    getline(iss, distance, ',') &&
                    getline(iss, line, ',')) {
                    int originID = stoi(origin);
                    int destinationID = stoi(destination);
                    graph.Graph::addVertex(originID);
                    graph.Graph::addVertex(destinationID);
                    graph.addBidirectionalEdge(originID, destinationID, stod(distance));
                }
            }
        }
    }

    inputFile.close();
}

void FileReader::loadGraph(const std::string &edgesfilename,const std::string &nodesfilename,TSPGraph &graph) {
    loadNodes(nodesfilename);
    ifstream inputFile(edgesfilename);

    if(inputFile.is_open())
    {
        string line;
        getline(inputFile, line);

        while(getline(inputFile,line))
        {
            istringstream iss(line);
            string origin,destination,distance;

            if( getline(iss, origin,',') &&
                getline(iss, destination,',') &&
                getline(iss, distance,'\r'))
            {
                int originID = stoi(origin);
                int destinationID = stoi(destination);
                pair<double, double> originNode = nodes[originID];
                pair<double, double> destinationNode = nodes[destinationID];
                graph.addVertex(originID, originNode.first, originNode.second);
                graph.addVertex(destinationID, destinationNode.first, destinationNode.second);
                graph.addBidirectionalEdge(originID,destinationID,stod(distance));
            }
        }
    }

    inputFile.close();
}