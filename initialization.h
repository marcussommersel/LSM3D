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
    void addPoint(Point point){this->points.push_back(point);}
    double minDist(Point point);
    bool isInside(Point p0, Point p1);
    double signedDistance(Point p);
};



bool intercepts(Point p0, Point p1, Point p2, Point p3);