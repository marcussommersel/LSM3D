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

Velocity vortexVelocity(int M, int N, int P, vector<double> X, vector<double> Y, vector<double> Z);

Velocity simpleVelocity(int M, int N, int P);