#include <vector>
#include <iostream>
#include <ctime>
#include <chrono>
#include <fstream>
#include "initialization2D.h"
#include "initialization2DTests.h"
#include "initialization.h"
#include "initializationTests.h"
#include "schemes.h"
#include "schemesTest.h"
#include "vectorUtilities.h"
using namespace std;

int main(){ 

    chrono::steady_clock::time_point startTime = chrono::steady_clock::now();

    bool runTests = false;

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

    const int m = 128;
    const int n = 128; 
    const int p = 128;

    const double xStart = 0;
    const double xEnd = 1;
    const double yStart = 0;
    const double yEnd = 1;
    const double zStart = 0;
    const double zEnd = 1;

    vector<double> x = linspace(xStart, xEnd, m);
    vector<double> y = linspace(yStart, yEnd, n);
    vector<double> z = linspace(zStart, zEnd, p);

    double dx = x[1] - x[0];
    double dy = y[1] - y[0];
    double dz = z[1] - z[0];

    vector<double> phi;
    Point c = Point(0.35,0.35,0.35);
    double r = 0.15;

    signedDistanceField(phi, x, y, z, r, c, m, n, p);
    auto [ax, ay, az] = vortexVelocity(m, n, p, x, y, z);

    double CFL = 0.9;
    double T = 3.0;

    double dt = CFL/(vectorMax(vectorAbs(ax)/dx + vectorAbs(ay)/dy + vectorAbs(az)/dz));
    
    double dtau = 0.5*dx;
    int reinitSteps = 10;
    bool doReinit = true;
    int reinitFreq = 50;

    bool halfplot = true;
    int it = 0;

    // clock_t startTime;

    vector<string> plotTimes;

    saveScalarField(to_string(0.000000) + ".txt", phi, x, y, z, m, n, p);
    plotTimes.push_back(to_string(0.000000));

    for (double t = 0; t < T; t += dt){

        ++it;

        if (t + dt > T){
            dt = T - t;
        }

        if (t > 1.5 && halfplot){
            saveScalarField(to_string(t) + ".txt", phi, x, y, z, m, n, p);
            plotTimes.push_back(to_string(t));
            halfplot = false;
        }

        TVDRK3_weno(phi, ax, ay, az, m, n, p, dx, dy, dz, dt);

        if (doReinit && (reinitFreq%it == 0)){
            for (int i = 0; i < reinitSteps - 1; ++i){
                TVDRK3_godunov_reinit(phi, m, n, p, dx, dy, dz, dtau);
            }
        }

        if (it%10 == 0){
            cout << "Iteration: " << it << endl;
            chrono::steady_clock::time_point currentTime = chrono::steady_clock::now();
            cout << "Elapsed time: " << (chrono::duration_cast<chrono::seconds>(currentTime - startTime).count())  << " s." << endl;
            // cout << "Iteration time: " << ((double)(clock() - startTime))/CLOCKS_PER_SEC << endl;
            cout << "t = " << t << endl;


            


        }
        if (it%50 == 0){
            saveScalarField(to_string(t) + ".txt", phi, x, y, z, m, n, p);
            plotTimes.push_back(to_string(t));
        }
    }

    saveScalarField(to_string(T) + ".txt", phi, x, y, z, m, n, p);
    plotTimes.push_back(to_string(T));

    ofstream file;
    file.open("plotTimes.txt");
    if (!file.is_open()){cerr << "could not open file." << endl;}
    for (int i = 0; i < plotTimes.size(); ++i){
        file << plotTimes[i] << endl;
    }
    file.close();
    return 0;
}