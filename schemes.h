#pragma once

#include <vector>
#include <array>
#include <tuple>
#include <functional>
#include <cmath>
#include "vectorUtilities.h"

using namespace std;

struct Derivative {
    vector<double> x;
    vector<double> y;
    vector<double> z;
};

Derivative upwind(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz);

Derivative weno(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int const N, int const P, double dx, double dy, double dz);

Derivative godunov(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int const N, int const P, double dx, double dy, double dz);

void euler_upwind(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt);

void TVDRK3_upwind(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt);

void TVDRK3_weno(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt);

void euler_weno(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt);

void TVDRK3_godunov_reinit(vector<double> &arr, vector<double> X, vector<double> Y, vector<double> Z, int M, int N, int P, double dx, double dy, double dz, double dt, vector<double> phi0);

void euler_upwind_reinit(vector<double> &arr, int M, int N, int P, double dx, double dy, double dz, double dt, const vector<double> &phi0);

void second_Order_Reinit(vector<double> &arr, vector<double> X, vector<double> Y, vector<double> Z, int M, int N, int P, double dx, double dy, double dz, double dt, const vector<double> &phi0);

int sign(double num);

double minAbs(double a, double b);

double minMod(double a, double b);
