#pragma once

#include "initialization.h"

void cubeTest(){
    double rho = 1;

    Points p;

    for (int i = 0; i < 361; i += 10){
        for (int j = 0; j < 181; j += 10){
            double theta = i*PI/180.0;
            double phi = j*PI/180.0;

            double x = rho*cos(theta)*sin(phi);
            double y = rho*sin(theta)*sin(phi);
            double z = rho*cos(phi);

            p.addPoint(Point(x, y, z));
        }
    }

    p.saveMatrix("file.txt");
    cout << endl;






}