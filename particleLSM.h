#pragma once
#include <vector>
#include <random>
#include "initialization.h"
#include "schemes.h"

// particle class with coordinates, radius and a bool, where true is a positive particle, and false is a negative particle
class Particle 
{
public:
    double x;
    double y;
    double z;
    double r;
    bool positive;
    Particle(double x1, double y1, double z1){x = x1; y = y1; z = z1;}
    Particle(){x = 0; y = 0; z = 0;}
    Particle(double x1, double y1, double z1, double r1, bool pos){x = x1; y = y1; z = z1; r = r1; positive = pos;}
};

// initializes particles in a cell where (x0, y0, z0) are the coordinates of the cell closest ot origo.
// dx, dy, dz are the grid spacing. X, Y, Z are all grid nodes in the computational domain.
// phi is the signed distance field. normal is the normal-vector. M, N, P is the number of grid nodes in each direction.
// numParticles are the number of particles of each type to be initialized in the cell.
vector<Particle> initializeParticles(double x0, double y0, double z0, double dx, double dy, double dz,
    vector<double> &X, vector<double> &Y, vector<double> &Z, vector<double> &phi, Derivative &normal,
    int M, int N, int P, int numParticles);

// returns the corrected values of the signed distance field for each corner of the cells.
// the xi, yj, zk values are the coordinates of the corners of the cell, where x0 < x1, y0 < y1, z0 < z1.
// the phiijk are the values of the signed distance fiield at (i,j,k)
vector<double> correctInterface(Particle p, double x0, double x1, double y0, double y1, double z0, double z1,
    double phi000, double phi100, double phi110, double phi010, double phi001, double phi101, double phi111, double phi011, double phip);

// returns the interpolated value at (x, y, z)
// the xi, yj, zk values are the coordinates of the corners of the cell, where x0 < x1, y0 < y1, z0 < z1.
// the phiijk are the values of the signed distance fiield at (i,j,k)
double trilinearInterpolation(double x, double y, double z, double x0, double x1, double y0, double y1, double z0, double z1,
    double c000, double c100, double c110, double c010, double c001, double c101, double c111, double c011);

// returns the normal vectors Nx, Ny, Nz for the signed distance field
Derivative normal(vector<double> &arr, double dx, double dy, double dz, double M, double N, double P);

// saves the coordinates of all the particles to a .txt-files
void plotParticles(string filename, vector<Particle> particles);
