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
