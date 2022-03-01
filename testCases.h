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

Velocity vortexVelocity(int M, int N, int P, vector<double> U, vector<double> V, vector<double> W);

Velocity simpleVelocity(int M, int N, int P);