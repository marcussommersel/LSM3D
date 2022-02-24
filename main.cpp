#include <vector>
#include <iostream>
#include "initialization2D.h"
#include "initialization2DTests.h"
#include "initialization.h"
#include "initializationTests.h"
#include "schemes.h"
#include "schemesTest.h"
#include "vectorUtilities.h"
using namespace std;

int main(){ 

    bool runTests = true;

    if (runTests)
    {
        // minDistTest();
        // interceptsTest();
        // isInsideTest();
        // signedDistanceTest();
        // simple2DAdvectionTest();
        // cubeTest();
        // isInsideSphereTest();
        // signedDistanceSphereTest();
        // signedDistanceFieldTest();
        euler_upwindTest();
    }

    // int n = 500;
    // int* array = new int[n*n*n];
    // for (int k = 0; k < n; ++k){
    //     for (int j = 0; j < n; ++j){
    //         for (int i = 0; i < n; ++i){
    //             array[i + j*n + k*n*n] = i + j + k;
    //             cout << array[i + j*n + k*n*n] << endl;
    //             cout << i << " " << j << " " << k << endl;
    //         }
    //     }
    // }



    // const int n = 10;
    // double x[n];
    // double y[n];


    // for (int i = 0; i < n; ++i){
    //     x[i] = i/(double)n;
    // }

    // for (int j = 0; j < n; ++j){
    //     y[j] = j/(double)n;
    // }

    // double u[n][n];
    // double v[n][n];

    // double phi[n][n];
    // vector<double> init[2];
    // for (int i = 0; i < n; ++i){
    //     if (x[i] > 0.1 && x[i] < 0.8){
    //         init[0].push_back(x[i] + 1.0/n);
    //     } 
    //     else {
    //         init[0].push_back(0);
    //     }
    //     cout << init[0][i] << " ";
    // }

    // for (int j = 0; j < n; ++j){
    //     if (y[j] > 0.2 && y[j] < 0.8){
    //         init[1].push_back(y[j] + 1.0/n);
    //     } 
    //     else {
    //         init[1].push_back(0);
    //     }
    // }

    // // double phi0[n][n] = sqrt(x^2 + init[0]);

    // for (int i = 0; i < n; ++i){
    //     for (int j = 0; j < n; ++i){
    //         // phi = 
    //     }
    // }

    // for (int i = 0; i < n; ++i){
    //     for (int j = 0; j < n; ++j){
    //         u[i][j] = 1;
    //         v[i][j] = 1;
    //         // if ((x[i] == 0.2 || x[i] == 0.8) && (y[j] > 0.2 && y[j] < 0.80)){
    //         // init[i][j] = x[];
    //         // }
    //     }
    // }
    return 0;
}