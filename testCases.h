#pragma once

#include <vector>
#include <tuple>
#include <cmath>

#include "initialization.h"

using namespace std;

#define PI 3.14159265

struct Velocity {
    vector<double> x;
    vector<double> y;
    vector<double> z;
};

Velocity vortexVelocity(int M, int N, int P, vector<double> X, vector<double> Y, vector<double> Z, double t, double T);

Velocity shearedSphereVelocity(int M, int N, int P, vector<double> X, vector<double> Y, vector<double> Z, double t, double T);

Velocity simpleVelocity(int M, int N, int P);

double volume(vector<double> &phi, double dx, double dy, double dz);

double surfaceArea(vector<double> &phi, double dx, double dy, double dz, double M, double N, double P);

double interfaceError(vector<double> &phi0, vector<double> &phi, double dx, double dy, double dz, double M, double N, double P);

double massError(vector<double> &phi, double dx, double dy, double dz, double M, double N, double P);