#pragma once

#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

#define PI 3.14159265

class Point 
{
public:
    double x;
    double y;
    double z;
    Point(double x1, double y1, double z1){x = x1; y = y1; z = z1;}
};

class Points
{
public:
    vector<Point> points;
    void addPoint(Point point){this->points.push_back(point);}
    void saveMatrix(string filename);
    // double minDist(Point point);
    // bool isInside(Point p0, Point p1);
    // double signedDistance(Point p);
};



// bool planeIntercepts(Point p0, Point p1, Point p2, Point p3);