//
// Created by jonas on 03/05/2024.
//

#ifndef DA2324_PRJ2_G11_3_TSPVERTEX_H
#define DA2324_PRJ2_G11_3_TSPVERTEX_H

#include "Graph.h"

class TSPVertex : public Vertex<int>{
public:
    TSPVertex();
    TSPVertex(int id,double longitude, double latitude);
    double getLongitude() const;
    double getLatitude() const;

    void setLongitude(double longitude);
    void setLatitude(double latitude);
protected:
    double longitude;
    double latitude;
};


#endif //DA2324_PRJ2_G11_3_TSPVERTEX_H
