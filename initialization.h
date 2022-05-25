#pragma once

#include <vector>
#include <iostream>
#include <array>
#include <cmath>
#include <fstream>
using namespace std;

#define PI 3.14159265

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

double length(Point const &p0, Point const &p1);

bool isInsideSphere(double r, Point c, Point p);

double signedDistanceSphere(double r, Point c, Point p);

void signedDistanceField(vector<double> &arr, vector<double> x, vector<double> y, vector<double> z, double r, Point c, int M, int N, int P);

vector<double> linspace(double start, double end, int n);

void saveScalarField(string filename, vector<double> const &arr, vector<double> x, vector<double> y, vector<double> z, int M, int N, int P);
