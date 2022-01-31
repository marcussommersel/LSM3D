#include "initialization2D.h"

double Points2D::minDist(Point2D point)
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

bool intercepts(Point2D p0, Point2D p1, Point2D p2, Point2D p3){
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

    double x = (b2*c1 - b1*c2)/denominator;
    double y = (a1*c2 - a2*c1)/denominator;

    double x0 = (x - p0.x)/(p1.x - p0.x);
    double y0 = (y - p0.y)/(p1.y - p0.y);
    double x1 = (x - p2.x)/(p3.x - p2.x);
    double y1 = (y - p2.y)/(p3.y - p2.y);

    if (((x0 >= 0 && x0 <= 1) || (y0 >= 0 && y0 <= 1)) && ((x1 >= 0 && x1 <= 1) || (y1 >= 0 && y1 <= 1))){
        return true;
    } else {
        return false;
    }
}

bool Points2D::isInside(Point2D p0, Point2D p1){
    int num = 0;
    for (int i = 0; i < this->points.size()-1; ++i) {
        if (intercepts(p0, p1, this->points[i], this->points[i+1])){
            ++num;
            if (p0.y == this->points[i+1].y){
                ++i;
            }
        }
    }
    if (num%2 == 0) {
        return false;
    } else {
        return true;
    }
}

double Points2D::signedDistance(Point2D p){
    double dist = this->minDist(p);
    if (isInside(p, Point2D(1000, p.y))){
        dist = -dist;
    }
    return dist;
}