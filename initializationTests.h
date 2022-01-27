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
    p.addPoint(Point(0, 0));
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
    p.addPoint(Point(0, 0));
    cout.precision(2);
    for (int j = 8; j > -5; --j){
        for (int i = 8; i > -5; --i){
            cout << fixed << p.signedDistance(Point(i,j)) << " ";
        }
        cout << endl;
    }
}

void simple2DAdvectionTest(){
    cout << "simple2DAdvectionTest: " << endl;
    Points p;
    for (int i = 2; i < 8; ++i){
        p.addPoint(Point(i, 2));
    }
    for (int i = 3; i < 8; ++i){
        p.addPoint(Point(7, i));
    }
    for (int i = 6; i > 1; --i){
        p.addPoint(Point(i, 7));
    }
    for (int i = 6; i > 1; --i){
        p.addPoint(Point(2, i));
    }
    cout.precision(2);
    int n = 10;
    int m = 10;

    for (int j = m; j > -1; --j){
        for (int i = 0; i < n+1; ++i){
            if (p.signedDistance(Point(i,j)) > 0){
                cout << " ";
            }
            cout << fixed << p.signedDistance(Point(i,j)) << " ";
        }
        cout << endl;
    }
    cout << endl;

    double phi[n][m];
    for (int i = 0; i < n+1; ++i){
        for (int j = 0; j < m+1; ++j){
            phi[i][j] = (p.signedDistance(Point(i,j)));
        }
    }
    for (int j = m; j > -1; --j){
        for (int i = 0; i < n+1; ++i){
            if (phi[i][j] > 0){
                cout << " ";
            }
            cout << fixed << phi[i][j] << " ";
        }
        cout << endl;
    }

    cout << endl;
    int it = 5;
    double phiNew[n][m];
    for (int k = 0; k < it ; ++k){
        for (int i = 0; i < n; ++i){
            for (int j = 0; j < m; ++j){
                    phiNew[i][j] = phi[i][j] - (phi[i+1][j] - phi[i][j]);
            }
        }
        memcpy(phi, phiNew, sizeof(phi));
        for (int j = m; j > -1; --j){
            for (int i = 0; i < n+1; ++i){
                if (phi[i][j] > 0){
                    cout << " ";
                }
                cout << fixed << phi[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}