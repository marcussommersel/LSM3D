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

void isInsideSphereTest(){
    Point c0 = Point(0,0,0);
    Point c1 = Point(1,2,1);
    Point p0 = Point(1,0,0);
    Point p1 = Point(3,2,4);
    double r0 = 5;
    double r1 = 0.5;
    double r2 = 1;
    cout << isInsideSphere(r0, c0, p0) << endl;
    cout << isInsideSphere(r1, c0, p0) << endl;
    cout << isInsideSphere(r0, c1, p0) << endl;
    cout << isInsideSphere(r0, c1, p1) << endl;
}