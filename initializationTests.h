#pragma once

#include "initialization.h"

void minDistTest(){
    cout << "minDistTest: Should return 0.141421" << endl;
    Points p;
    p.addPoint(Point(0, 0));
    p.addPoint(Point(1, 0));
    p.addPoint(Point(2, 0));
    p.addPoint(Point(0, 1));
    p.addPoint(Point(1, 1));
    p.addPoint(Point(2, 1));
    p.addPoint(Point(3, 1));
    cout << p.minDist(Point(1.1, 0.9)) << endl;
}

void interceptsTest(){
    cout << "interceptsTest: Should return 1, 0, 1, 1, 0" << endl;
    cout << intercepts(Point(-1, 0), Point(1, 0), Point(0, 1), Point(0, -1)) << endl;
    cout << intercepts(Point(-1, 0), Point(1, 0), Point(-1, 1), Point(1, 1)) << endl;
    cout << intercepts(Point(-1, 0), Point(0, 0), Point(0, 1), Point(0, -1)) << endl;
    cout << intercepts(Point(-1, 0), Point(1, 0), Point(0.0, 0.5), Point(1, -0.5)) << endl;
    cout << intercepts(Point(-1, 0), Point(1, 0), Point(0, 1), Point(0, 2)) << endl;
}

void isInsideTest(){
    cout << "isInsideTest: Should return 1, 0" << endl;
    Points p;
    p.addPoint(Point(0, 0));
    p.addPoint(Point(1, 0));
    p.addPoint(Point(1, 1));
    p.addPoint(Point(0, 1));
    cout << p.isInside(Point(0.5, 0.5), Point(2, 0.5)) << endl;
    cout << p.isInside(Point(-1, 0.5), Point(2, 0.5)) << endl;
}

void signedDistanceTest(){
    cout << "signedDistanceTest: " << endl;
    Points p;
    p.addPoint(Point(0, 0));
    p.addPoint(Point(4, 0));
    p.addPoint(Point(4, 4));
    p.addPoint(Point(0, 4));
    cout.precision(2);
    for (int i = -4; i < 9; ++i){
        for (int j = -4; j < 9; ++j){
            cout << fixed << p.signedDistance(Point(i,j)) << " ";
        }
        cout << endl;
    }
}