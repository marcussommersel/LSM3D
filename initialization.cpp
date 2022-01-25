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

// bool Points::isInside(Point point){
//     double dx = 0;
//     for (auto const& p: this->points){
//         if (this->points[0].x != p.x){
//             dx = abs(p.x - this->points[0].x);
//             break;
//         }
//     }

// }