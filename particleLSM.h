#pragma once
#include <vector>
#include <random>
#include "initialization.h"
#include "schemes.h"

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
    Particle(double x1, double y1, double z1, double r1, bool pos){x = x1; y = y1; z = y1; r = r1; positive = pos;}
};

vector<Particle> initializeParticles(double x0, double y0, double z0, double dx, double dy, double dz, vector<double> &phi, Derivative &normal,
    int i, int j, int k, int M, int N, int P);

vector<double> correctInterface(Particle p, double x0, double x1, double y0, double y1, double z0, double z1,
    double phi000, double phi100, double phi110, double phi010, double phi001, double phi101, double phi111, double phi011);

double trilinearInterpolation(double x, double y, double z, double x0, double x1, double y0, double y1, double z0, double z1,
    double c000, double c100, double c110, double c010, double c001, double c101, double c111, double c011);

Derivative normal(vector<double> &arr, double dx, double dy, double dz, double M, double N, double P);