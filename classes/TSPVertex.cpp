#include "TSPVertex.h"

TSPVertex::TSPVertex() : Vertex<int>() {
    this->longitude = 0.0;
    this->latitude = 0.0;
}

TSPVertex::TSPVertex(int id, double longitude, double latitude) : Vertex<int>(id) {
    this->longitude = longitude;
    this->latitude = latitude;
}


double TSPVertex::getLongitude() const {
    return this->longitude;
}

double TSPVertex::getLatitude() const{
    return this->latitude;
}

void TSPVertex::setLongitude(double longitude) {
    this->longitude = longitude;
}

void TSPVertex::setLatitude(double latitude) {
    this->latitude = latitude;
}

