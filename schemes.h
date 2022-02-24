#pragma once

#include <vector>
#include <tuple>
#include <functional>
#include "vectorUtilities.h"

using namespace std;

struct Derivative {
    vector<double> x;
    vector<double> y;
    vector<double> z;
};

Derivative upwind(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz);

void euler_upwind(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt);
