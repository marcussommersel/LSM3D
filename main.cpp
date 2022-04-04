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
#include "particleLSM.h"
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
        // euler_upwindTest();
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
    int reinitFreq = 5;

    bool halfplot = true;
    int it = 0;

    vector<string> plotTimes;

    saveScalarField(to_string(0.000000) + ".txt", phi, x, y, z, m, n, p);
    plotTimes.push_back(to_string(0.000000));

    vector<Particle> particles;

    // Initializing particles
    Derivative norm = normal(phi, dx, dy, dz, m, n, p); 
    for (int k = 0; k < p; ++k){
        for (int j = 0; j < n; ++j){
            for (int i = 0; i < m; ++i){
                if (abs(phi[i + j*n + k*p*p] < 3*max(dx, max(dy,dz)))){
                    vector<Particle> newParticles = initializeParticles(x[i], y[j], z[k], dx, dy, dz, phi, norm, i, j, k, m, n, p);
                    particles.insert(particles.end(), newParticles.begin(), newParticles.end()); // fix not double cells
                }
            }
        }
    }

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

        // Advect particles
        for (int a = 0; a < particles.size(); ++a){
            double currentI;
            double currentJ;
            double currentK;
            for (int b = 0; b < m; ++b){
                if (particles[a].x >= x[b] && particles[a].x <= x[b+1]){
                    currentI = b;
                }
                if (particles[a].y >= y[b] && particles[a].y <= y[b+1]){
                    currentJ = b;
                }
                if (particles[a].z >= z[b] && particles[a].z <= z[b+1]){
                    currentK = b;
                }
            }

            double Up = trilinearInterpolation(particles[a].x, particles[a].y, particles[a].z, x[currentI], x[currentI+1], y[currentJ], y[currentJ+1],
                z[currentK], z[currentK+1], ax[currentI+currentJ*n+currentK*p*p], ax[(currentI+1)+currentJ*n+currentK*p*p], ax[(currentI+1)+(currentJ+1)*n+currentK*p*p], 
                ax[currentI+(currentJ+1)*n+currentK*p*p], ax[currentI+currentJ*n+(currentK+1)*p*p], ax[(currentI+1)+currentJ*n+(currentK+1)*p*p], 
                ax[(currentI+1)+(currentJ+1)*n+(currentK+1)*p*p], ax[currentI+(currentJ+1)*n+(currentK+1)*p*p]);

            double Vp = trilinearInterpolation(particles[a].x, particles[a].y, particles[a].z, x[currentI], x[currentI+1], y[currentJ], y[currentJ+1],
                z[currentK], z[currentK+1], ay[currentI+currentJ*n+currentK*p*p], ay[(currentI+1)+currentJ*n+currentK*p*p], ay[(currentI+1)+(currentJ+1)*n+currentK*p*p], 
                ay[currentI+(currentJ+1)*n+currentK*p*p], ay[currentI+currentJ*n+(currentK+1)*p*p], ay[(currentI+1)+currentJ*n+(currentK+1)*p*p], 
                ay[(currentI+1)+(currentJ+1)*n+(currentK+1)*p*p], ay[currentI+(currentJ+1)*n+(currentK+1)*p*p]);

            double Wp = trilinearInterpolation(particles[a].x, particles[a].y, particles[a].z, x[currentI], x[currentI+1], y[currentJ], y[currentJ+1],
                z[currentK], z[currentK+1], az[currentI+currentJ*n+currentK*p*p], az[(currentI+1)+currentJ*n+currentK*p*p], az[(currentI+1)+(currentJ+1)*n+currentK*p*p], 
                az[currentI+(currentJ+1)*n+currentK*p*p], az[currentI+currentJ*n+(currentK+1)*p*p], az[(currentI+1)+currentJ*n+(currentK+1)*p*p], 
                az[(currentI+1)+(currentJ+1)*n+(currentK+1)*p*p], az[currentI+(currentJ+1)*n+(currentK+1)*p*p]);

            particles[a].x = particles[a].x + dt*Up;
            particles[a].y = particles[a].y + dt*Vp;
            particles[a].z = particles[a].z + dt*Wp;

        }

        if (doReinit && (it%reinitFreq == 0)){ // it!=0?
            vector<double> phi0 = phi;
            for (int i = 0; i < reinitSteps - 1; ++i){
                euler_upwind_reinit(phi, m, n, p, dx, dy, dz, dtau, phi0);
            }
        }

        for (int a = 0; a < particles.size(); ++a){
            double currentI;
            double currentJ;
            double currentK;
            for (int b = 0; b < m; ++b){
                if (particles[a].x >= x[b] && particles[a].x <= x[b+1]){
                    currentI = b;
                }
                if (particles[a].y >= y[b] && particles[a].y <= y[b+1]){
                    currentJ = b;
                }
                if (particles[a].z >= z[b] && particles[a].z <= z[b+1]){
                    currentK = b;
                }
            }
            double phip = trilinearInterpolation(particles[a].x, particles[a].y, particles[a].z, x[currentI], x[currentI+1], y[currentJ], y[currentJ+1],
                z[currentK], z[currentK+1], phi[currentI+currentJ*n+currentK*p*p], phi[(currentI+1)+currentJ*n+currentK*p*p], phi[(currentI+1)+(currentJ+1)*n+currentK*p*p], 
                phi[currentI+(currentJ+1)*n+currentK*p*p], phi[currentI+currentJ*n+(currentK+1)*p*p], phi[(currentI+1)+currentJ*n+(currentK+1)*p*p], 
                phi[(currentI+1)+(currentJ+1)*n+(currentK+1)*p*p], phi[currentI+(currentJ+1)*n+(currentK+1)*p*p]);

            if (((phip < 0 && particles[a].positive) || (phip > 0 && !particles[a].positive) ) && (abs(phip) > particles[a].r)){ // correct interface
                vector<double> phiCorrected = 
                    correctInterface(particles[a], x[currentI], x[currentI+1], y[currentJ], y[currentJ+1], z[currentK], z[currentK],
                    phi[(currentI)+(currentJ)*n+(currentK)*p*p], phi[(currentI+1)+(currentJ)*n+(currentK)*p*p], phi[(currentI+1)+(currentJ+1)*n+(currentK)*p*p],
                    phi[(currentI)+(currentJ+1)*n+(currentK)*p*p], phi[(currentI)+(currentJ)*n+(currentK+1)*p*p], phi[(currentI+1)+(currentJ)*n+(currentK+1)*p*p],
                    phi[(currentI+1)+(currentJ+1)*n+(currentK+1)*p*p], phi[(currentI)+(currentJ+1)*n+(currentK+1)*p*p]);

                    phi[(currentI)+(currentJ)*n+(currentK)*p*p] = phiCorrected[0];
                    phi[(currentI+1)+(currentJ)*n+(currentK)*p*p] = phiCorrected[1];
                    phi[(currentI+1)+(currentJ+1)*n+(currentK)*p*p] = phiCorrected[2];
                    phi[(currentI)+(currentJ+1)*n+(currentK)*p*p] = phiCorrected[3];
                    phi[(currentI)+(currentJ)*n+(currentK+1)*p*p] = phiCorrected[4];
                    phi[(currentI+1)+(currentJ)*n+(currentK+1)*p*p] = phiCorrected[5];
                    phi[(currentI+1)+(currentJ+1)*n+(currentK+1)*p*p] = phiCorrected[6];
                    phi[(currentI)+(currentJ+1)*n+(currentK+1)*p*p] = phiCorrected[7];
                    
            } 
        }

        if (it%10 == 0){
            cout << "Iteration: " << it << endl;
            chrono::steady_clock::time_point currentTime = chrono::steady_clock::now();
            cout << "Elapsed time: " << (chrono::duration_cast<chrono::seconds>(currentTime - startTime).count())  << " s." << endl;
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