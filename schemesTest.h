#pragma once

#include "schemes.h"
#include "initialization.h"
#include "testCases.h"

using namespace std;

void euler_upwindTest(){
    const int m = 32;
    const int n = 32; 
    const int p = 32;
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
    auto [ax, ay, az] = vortexVelocity(m, n, p, x, y, z);

    double dt = 1;
    double dx = 1;
    double dy = 1;
    double dz = 1;

    for (double t = 0; t < 5; t += dt){
        saveScalarField(to_string(t) + ".txt", phi, x, y, z, m, n, p);
        euler_weno(phi, ax, ay, az, m, n, p, dx, dy, dz, dt);
    } 
}