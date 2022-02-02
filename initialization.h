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
    Point(){x = 0; y = 0; z = 0;}
    Point operator+(Point const &p);
    void operator=(Point const &p);
    double length(Point const &p);
};

class Points
{
public:
    vector<Point> points;
    void addPoint(Point point){this->points.push_back(point);}
    void saveMatrix(string filename);
    Point operator+(Point const &p);
    
    // double minDist(Point point);
    // bool isInside(Point p0, Point p1);
    // double signedDistance(Point p);
};

class Planes
{
private:
    void addPlane(Points points){this->planes.push_back(points);}
public:
    Planes(Point p0, Point p1, Point p2, Point p3);
    vector<Points> planes;
    void findPlanes(Points pts);
};

bool isInsideSphere(double r, Point c, Point p);
// bool planeIntercepts(Point p0, Point p1, Point p2, Point p3);

double signedDistanceSphere(double r, Point c, Point p);