#pragma once

#include <vector>
#include <iostream>
#include <array>
#include <cmath>
#include <fstream>
using namespace std;

#define PI 3.14159265

// class for a 3D point
class Point 
{
public:
    double x;
    double y;
    double z;
    Point(double x1, double y1, double z1){x = x1; y = y1; z = z1;}
    Point(){x = 0; y = 0; z = 0;}
    Point operator+(Point const &p);
    void operator=(Point const &p);
};

// returns length between two points
double length(Point const &p0, Point const &p1);

// check if a point is within a sphere of center c and radius r
bool isInsideSphere(double r, Point c, Point p);

// returns the signed distance from a point to the surface of a sphere of center c and radius r
double signedDistanceSphere(double r, Point c, Point p);

// generates a signed distance field for all points in [xmin, xmax] * [ymin, ymax] * [zmin, zmax] with reference to a sphere of center c and radius r
void signedDistanceField(vector<double> &arr, vector<double> x, vector<double> y, vector<double> z, double r, Point c, int M, int N, int P);

// returns a vector of n indexes with equally spaced values from start to end
vector<double> linspace(double start, double end, int n);

// saves a signed distance field to .txt-file
void saveScalarField(string filename, vector<double> const &arr, vector<double> x, vector<double> y, vector<double> z, int M, int N, int P);
