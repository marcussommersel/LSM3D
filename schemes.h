#pragma once

#include <vector>
#include <array>
#include <tuple>
#include <functional>
#include <cmath>
#include "vectorUtilities.h"

using namespace std;

// a deriative vector with a value in each direction in 3D
struct Derivative {
    vector<double> x;
    vector<double> y;
    vector<double> z;
};

// first-order upwind scheme
Derivative upwind(vector<double> &phi, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz);

// WENO scheme. Third-order accurate and fifth-order accurate in smooth regions
Derivative weno(vector<double> &phi, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int const N, int const P, double dx, double dy, double dz);

// Godunov scheme used for the reinitialization equation
Derivative godunov(vector<double> &phi, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int const N, int const P, double dx, double dy, double dz);

// first-order explicit Euler scheme used with the upwind scheme
void euler_upwind(vector<double> &phi, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt);

// third-order TVDRK scheme used with the upwind scheme
void TVDRK3_upwind(vector<double> &phi, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt);

// third-order TVDRK scheme used with the WENO scheme
void TVDRK3_weno(vector<double> &phi, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt);

// first-order explicit Euler scheme used with the weno scheme
void euler_weno(vector<double> &phi, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt);

// third-order TVDRK scheme used with the Godunov scheme to solve the reinitialization equation
void TVDRK3_godunov_reinit(vector<double> &phi, int M, int N, int P, double dx, double dy, double dz, double dt, vector<double> phi0);

// first-order explicit Euler scheme used with the Godunov scheme to solve the reinitialization equation
void euler_godunov_reinit(vector<double> &phi, int M, int N, int P, double dx, double dy, double dz, double dt, const vector<double> &phi0);

// sign function that returns 1 for a positive value, -1 for a negative value, and 0 for a value of 0
int sign(double num);
