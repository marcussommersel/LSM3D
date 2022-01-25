#pragma once

#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

class Point 
{
public:
    double x;
    double y;
    Point(double x1, double y1){x = x1; y = y1;}
};

class Points
{
public:
    vector<Point> points;
    double minDist(Point point);
    void addPoint(Point point){this->points.push_back(point);}
};