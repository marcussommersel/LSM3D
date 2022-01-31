#include "initialization.h"
#include <iostream>
#include <fstream>

void Points::saveMatrix(string filename){
    ofstream file;
    file.open(filename);
    if (file.is_open()) {
        file << "x, y, z" << endl;
        for (int i = 0; i < this->points.size(); ++i){
            file << points[i].x << ", " << points[i].y << ", " << points[i].z << endl;
        }
    } else {cout << "could not open file." << endl;}

    file.close();
}

// double Points::minDist(Point point)
// {
//     double min = INFINITY;
//     for (auto const& p: this->points){
//         double dist = sqrt((point.x - p.x)*(point.x - p.x) + (point.y - p.y)*(point.y - p.y));
//         if (dist < min){
//             min = dist;
//         }
//     }
//     return min;
// }

bool planeIntercepts(Point p0, Point p1, Point p2, Point p3){
    return 0;
}

// bool Points::isInside(Point p0, Point p1){
//     int num = 0;
//     for (int i = 0; i < this->points.size()-1; ++i) {
//         if (planeIntercepts(p0, p1, this->points[i], this->points[i+1])){
//             ++num;
//             if (p0.y == this->points[i+1].y){
//                 ++i;
//             }
//         }
//     }
//     if (num%2 == 0) {
//         return false;
//     } else {
//         return true;
//     }
// }

// double Points::signedDistance(Point p){
//     double dist = this->minDist(p);
//     if (isInside(p, Point(1000, p.y))){
//         dist = -dist;
//     }
//     return dist;
// }