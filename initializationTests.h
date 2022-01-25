#pragma once

#include "initialization.h"

void minDistTest(){
    Points p;
    p.addPoint(Point(0, 0));
    p.addPoint(Point(1, 0));
    p.addPoint(Point(2, 0));
    p.addPoint(Point(0, 1));
    p.addPoint(Point(1, 1));
    p.addPoint(Point(2, 1));
    p.addPoint(Point(3, 1));
    cout << p.minDist(Point(1.1, 0.9));
}