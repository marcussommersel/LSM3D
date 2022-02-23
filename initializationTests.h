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

void signedDistanceSphereTest(){
    Point p0 = Point(5, 5, 5);
    Point p1 = Point(5, 0, 5);
    Point p2 = Point(5, 3, 5);
    Point p3 = Point(5, 4, 5);
    cout << length(p1, p0) << endl;
    cout << signedDistanceSphere(2, p0, p1) << endl;
    cout << signedDistanceSphere(2, p0, p2) << endl;
    cout << signedDistanceSphere(2, p0, p3) << endl;
}

void signedDistanceFieldTest(){
    const int m = 10;
    const int n = 10; 
    const int p = 10;
    const double xStart = 0;
    const double xEnd = 10;
    const double yStart = 0;
    const double yEnd = 10;
    const double zStart = 0;
    const double zEnd = 10;
    vector<double> x = linspace(xStart, xEnd, m);
    vector<double> y = linspace(yStart, yEnd, n);
    vector<double> z = linspace(zStart, zEnd, p);
    vector<double> phi;
    Point c = Point(5,5,5);
    double r = 3;

    signedDistanceField(phi, x, y, z, r, c, m, n, p);
    for (int k = 0; k < p; ++k){
        for (int j = 0; j < n; ++j){
            for (int i = 0; i < m; ++i){
                cout << phi[i + j*n + k*p*p] << " ";
            }
            cout << endl;
        }
        cout << endl << endl;
    }
    saveScalarField("scaleField.txt", phi, x, y, z, m, n, p);
}
