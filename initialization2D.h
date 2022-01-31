#pragma once

#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

class Point2D 
{
public:
    double x;
    double y;
    Point2D(double x1, double y1){x = x1; y = y1;}
};

class Points2D
{
public:
    vector<Point2D> points;
    void addPoint(Point2D point){this->points.push_back(point);}
    double minDist(Point2D point);
    bool isInside(Point2D p0, Point2D p1);
    double signedDistance(Point2D p);
};



bool lineIntercepts(Point2D p0, Point2D p1, Point2D p2, Point2D p3);