#include "initialization.h"

double Points::minDist(Point point)
{
    double min = INFINITY;
    for (auto const& p: this->points){
        double dist = sqrt((point.x - p.x)*(point.x - p.x) + (point.y - p.y)*(point.y - p.y));
        if (dist < min){
            min = dist;
        }
    }
    return min;
}

bool intercepts(Point p0, Point p1, Point p2, Point p3){
    double a1 = p1.y - p0.y;
    double b1 = p0.x - p1.x;
    double c1 = a1*p0.x + b1*p0.y;
    double a2 = p3.y - p2.y;
    double b2 = p2.x - p3.x;
    double c2 = a2*p2.x + b2*p2.y;
    double denominator  = a1*b2 - a2*b1;

    if (denominator == 0)
    {
        return false;
    } 
    else
    {
        return true;
    }

    // double x = (b2*c1 - b1*c2)/denominator;
    // double y = (a1*c2 - a1*c1)/denominator;
}

// int numIntercept(Point p0, Point p1, Point p2, Point p3){
//     int num = 0;

//     return num;
// }

bool Points::isInside(Point p0, Point p1){
    int num = 0;
    for (int i = 0; i < this->points.size()-2; ++i) {
        if (intercepts(p0, p1, this->points[i], this->points[i+1])){
            ++num;
        }
    }
    if (intercepts(p0, p1, this->points[points.size()-1], this->points[0])){
        ++num;
    }
    if (num%2 == 0) {
        return false;
    } else {
        return true;
    }
}

double Points::signedDistance(Point p){
    double dist = this->minDist(p);
    if (isInside(p, Point(INFINITY, p.y))){
        dist = -dist;
    }
    return dist;
}