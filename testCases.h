#pragma once

#include <vector>
#include <tuple>

#include "initialization.h"

using namespace std;

// auto vortexVelocity(int M, int N, int P);
struct Velocity {
    vector<double> x;
    vector<double> y;
    vector<double> z;
};

Velocity simpleVelocity(int M, int N, int P);