#pragma once

#include <vector>
#include <tuple>
#include <cmath>

#include "initialization.h"
#include "schemes.h"

using namespace std;

#define PI 3.14159265

// velocity vector with one value for each direction in 3D
struct Velocity {
    vector<double> x;
    vector<double> y;
    vector<double> z;
};

// velocity field for deformation in 3D. Taken from LeVeque (1996)
Velocity vortexVelocity(int M, int N, int P, vector<double> X, vector<double> Y, vector<double> Z, double t, double T);

// velocity field for deformation in 2D. Taken from Morgan and Waltz (2017)
Velocity shearedSphereVelocity(int M, int N, int P, vector<double> X, vector<double> Y, vector<double> Z, double t, double T);

// simple velocity field with u = v = w = 1 for all grid nodes
Velocity simpleVelocity(int M, int N, int P);

// returns the volume of the domain bounded by the zero contour in a signed distance field
double volume(vector<double> &phi, double dx, double dy, double dz);

// returns the surface area of the domain bounded by the zero contour in a signed distance field
double surfaceArea(vector<double> &phi, double dx, double dy, double dz, double M, double N, double P);

// error measure of the interface error
double interfaceError(vector<double> &phi0, vector<double> &phi, double dx, double dy, double dz, double M, double N, double P);

// error measure of the average mass error
double massError(vector<double> &phi, double dx, double dy, double dz, double M, double N, double P);
