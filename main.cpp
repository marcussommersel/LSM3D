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

    const int m = 50;
    const int n = 50; 
    const int p = 50;

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

    double initialVolume = volume(phi, dx, dy, dz); 
    
    // auto [ax, ay, az] = simpleVelocity(m, n, p);

    double CFL = 0.9;
    double T = 3.0;
    // double T = 0.75;

    double dtau = 0.5*dx;
    int reinitSteps = 5;
    bool doReinit = true;
    bool doParticle = true;
    int reinitFreq = 5;
    int plotFreq = 100;
    int reseedFreq = 100;
    int itmax = 2000;

    bool halfplot = true;

    vector<string> plotTimes;
    vector<string> plotTimesParticle;
    int nParticles = 32;

    saveScalarField(to_string(0.000000) + ".txt", phi, x, y, z, m, n, p);
    plotTimes.push_back(to_string(0.000000));


    vector<Particle> particles;
    double rmin = 0.1*min(dx, min(dy, dz));
    double rmax = 0.5*max(dx, max(dy, dz));
    double bmin = rmin;
    double bmax = 3.0*max(dx, max(dy, dz));

    // Initializing particles
    if (doParticle){
        Derivative norm = normal(phi, dx, dy, dz, m, n, p); // parse phi and only compute some normals instead?
        for (int k = 0; k < p; ++k){
            for (int j = 0; j < n; ++j){
                for (int i = 0; i < m; ++i){
                    if (abs(phi[i + j*n + k*p*p]) < 3*max(dx, max(dy,dz))){
                        vector<Particle> newParticles = initializeParticles(x[i], y[j], z[k], dx, dy, dz, phi, norm, i, j, k, m, n, p, nParticles);
                        particles.insert(particles.end(), newParticles.begin(), newParticles.end()); // fix not double cells
                    }
                }
            }
        }

        plotParticles(to_string(0.000000) + "particle.txt" , particles);
        plotTimesParticle.push_back(to_string(0.000000) + "particle");
    }

    double dt = 0;
    double t = 0;
    int numIt = 0;

    cout << "Setup complete." << endl;
    chrono::steady_clock::time_point currentTime = chrono::steady_clock::now();
    cout << "Elapsed time: " << (chrono::duration_cast<chrono::seconds>(currentTime - startTime).count())  << " s." << endl;

    for (int it = 1; it < itmax; ++it){

        auto [ax, ay, az] = vortexVelocity(m, n, p, x, y, z, t, T);
        dt = CFL/(vectorMax(vectorAbs(ax)/dx + vectorAbs(ay)/dy + vectorAbs(az)/dz));

        if (t + dt > T){
            dt = T - t;
        }

        t += dt;

        if (t > 1.5 && halfplot){
            saveScalarField(to_string(t) + ".txt", phi, x, y, z, m, n, p);
            plotTimes.push_back(to_string(t));
            halfplot = false;
        }

        TVDRK3_weno(phi, ax, ay, az, m, n, p, dx, dy, dz, dt);

        // Advect particles
        if (doParticle){
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

                double phip = trilinearInterpolation(particles[a].x, particles[a].y, particles[a].z, x[currentI], x[currentI+1], y[currentJ], y[currentJ+1],
                    z[currentK], z[currentK+1], phi[currentI+currentJ*n+currentK*p*p], phi[(currentI+1)+currentJ*n+currentK*p*p], phi[(currentI+1)+(currentJ+1)*n+currentK*p*p], 
                    phi[currentI+(currentJ+1)*n+currentK*p*p], phi[currentI+currentJ*n+(currentK+1)*p*p], phi[(currentI+1)+currentJ*n+(currentK+1)*p*p], 
                    phi[(currentI+1)+(currentJ+1)*n+(currentK+1)*p*p], phi[currentI+(currentJ+1)*n+(currentK+1)*p*p]);

                // Delete particles
                if ((((phip > bmax) && particles[a].positive) || ((phip > -bmax) && !particles[a].positive) || abs(phip) > 3*max(dx, max(dy,dz))) && (it%reseedFreq==0)){
                    particles.erase(particles.begin() + a);
                }
                // interface correction
                else if ((phip < 0 && particles[a].positive) || (phip > 0 && !particles[a].positive)  && (abs(phip) > particles[a].r)){ // correct interface
                    vector<double> phiCorrected = 
                        correctInterface(particles[a], x[currentI], x[currentI+1], y[currentJ], y[currentJ+1], z[currentK], z[currentK+1],
                        phi[(currentI)+(currentJ)*n+(currentK)*p*p], phi[(currentI+1)+(currentJ)*n+(currentK)*p*p], phi[(currentI+1)+(currentJ+1)*n+(currentK)*p*p],
                        phi[(currentI)+(currentJ+1)*n+(currentK)*p*p], phi[(currentI)+(currentJ)*n+(currentK+1)*p*p], phi[(currentI+1)+(currentJ)*n+(currentK+1)*p*p],
                        phi[(currentI+1)+(currentJ+1)*n+(currentK+1)*p*p], phi[(currentI)+(currentJ+1)*n+(currentK+1)*p*p], phip);

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

            // Initialize new particles
            if (it%reseedFreq == 0 && it != 0){
                vector<int> cellx;
                vector<int> celly;
                vector<int> cellz;
                vector<int> cellParticles;
                cellx.push_back((int)(particles[0].x/dx));
                celly.push_back((int)(particles[0].y/dy));
                cellz.push_back((int)(particles[0].z/dz));
                cellParticles.push_back(1);

                bool found = false;
                for (int a = 1; a < particles.size(); ++a){
                    for (int b = 0; b < cellx.size(); ++b){
                        if ((int)(particles[a].x/dx) == cellx[b] && (int)(particles[a].y/dy) == celly[b] && (int)(particles[a].z/dz) == cellz[b]){
                            cellParticles[b] += 1;
                            if (cellParticles[b] > nParticles){
                                particles.erase(particles.begin() + a);
                                a -= 1;
                            }
                            found = true;
                            break;
                        }
                    }
                    if (!found){
                        cellx.push_back((int)(particles[a].x/dx));
                        celly.push_back((int)(particles[a].y/dy));
                        cellz.push_back((int)(particles[a].z/dz));
                        cellParticles.push_back(1);
                        found = false;
                    }
                }

                Derivative norm = normal(phi, dx, dy, dz, m, n, p); // parse phi and only compute some normals instead?
                for (int a = 0; a < cellx.size(); ++a){
                    int num = nParticles - cellParticles[a];
                    if (num > 0){
                        vector<Particle> newParticles = initializeParticles(x[cellx[a]], y[celly[a]], z[cellz[a]], dx, dy, dz, phi, norm, cellx[a], celly[a], cellz[a], m, n, p, num);
                        particles.insert(particles.end(), newParticles.begin(), newParticles.end());
                    }
                }
            }
        }

        if (doReinit && (it%reinitFreq == 0)){
            vector<double> phi0 = phi;
            for (int i = 0; i < reinitSteps - 1; ++i){
                euler_upwind_reinit(phi, m, n, p, dx, dy, dz, dtau, phi0);
            }
        }

        if (doParticle){
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

                // if (((phip > bmax) && particles[a].positive) || ((phip > -bmax) && !particles[a].positive)){
                //     particles.erase(particles.begin() + a);
                // } else 
                if (((phip < 0 && particles[a].positive) || (phip > 0 && !particles[a].positive) ) && (abs(phip) > particles[a].r)){ // correct interface
                    vector<double> phiCorrected = 
                        correctInterface(particles[a], x[currentI], x[currentI+1], y[currentJ], y[currentJ+1], z[currentK], z[currentK+1],
                        phi[(currentI)+(currentJ)*n+(currentK)*p*p], phi[(currentI+1)+(currentJ)*n+(currentK)*p*p], phi[(currentI+1)+(currentJ+1)*n+(currentK)*p*p],
                        phi[(currentI)+(currentJ+1)*n+(currentK)*p*p], phi[(currentI)+(currentJ)*n+(currentK+1)*p*p], phi[(currentI+1)+(currentJ)*n+(currentK+1)*p*p],
                        phi[(currentI+1)+(currentJ+1)*n+(currentK+1)*p*p], phi[(currentI)+(currentJ+1)*n+(currentK+1)*p*p], phip);

                    phi[(currentI)+(currentJ)*n+(currentK)*p*p] = phiCorrected[0];
                    phi[(currentI+1)+(currentJ)*n+(currentK)*p*p] = phiCorrected[1];
                    phi[(currentI+1)+(currentJ+1)*n+(currentK)*p*p] = phiCorrected[2];
                    phi[(currentI)+(currentJ+1)*n+(currentK)*p*p] = phiCorrected[3];
                    phi[(currentI)+(currentJ)*n+(currentK+1)*p*p] = phiCorrected[4];
                    phi[(currentI+1)+(currentJ)*n+(currentK+1)*p*p] = phiCorrected[5];
                    phi[(currentI+1)+(currentJ+1)*n+(currentK+1)*p*p] = phiCorrected[6];
                    phi[(currentI)+(currentJ+1)*n+(currentK+1)*p*p] = phiCorrected[7];
                        
                }

                // Adjust radius
                if (sign(phip)*phip > rmax){
                    particles[a].r = rmax;
                } else if (sign(phip)*phip < rmin){
                    particles[a].r = rmin;
                } else {
                    particles[a].r = sign(phip)*phip;
                }

            }
        }

        if (it%10 == 0){
            cout << "Iteration: " << it << endl;
            chrono::steady_clock::time_point currentTime = chrono::steady_clock::now();
            cout << "Elapsed time: " << (chrono::duration_cast<chrono::seconds>(currentTime - startTime).count())  << " s." << endl;
            cout << "t = " << t << endl;
        }
        if (it%plotFreq == 0){
            saveScalarField(to_string(t) + ".txt", phi, x, y, z, m, n, p);
            plotTimes.push_back(to_string(t));

            if (doParticle){
                plotParticles(to_string(t) + "particle.txt" , particles);
                plotTimesParticle.push_back(to_string(t) + "particle");
            }
        }

        if (t == T){
            numIt = it;
            break;
        }

    }

    double endVolume = volume(phi, dx, dy, dz);
    double volumeChange = 100*(initialVolume-endVolume)/initialVolume;

    saveScalarField(to_string(T) + ".txt", phi, x, y, z, m, n, p);
    plotTimes.push_back(to_string(T));

    {
        ofstream file;
        file.open("plotTimes.txt");
        if (!file.is_open()){cerr << "could not open file." << endl;}
        for (int i = 0; i < plotTimes.size(); ++i){
            file << plotTimes[i] << endl;
        }
        file.close();
    }

    if (doParticle){
        ofstream file;
        file.open("plotTimesParticle.txt");
        if (!file.is_open()){cerr << "could not open file." << endl;}
        for (int i = 0; i < plotTimesParticle.size(); ++i){
            file << plotTimesParticle[i] << endl;
        }
        file.close();
    }

    {
        ofstream file;
        file.open("log.txt");
        if (!file.is_open()){cerr << "could not open file." << endl;}
        file << "Iterations: " << numIt << endl;
        chrono::steady_clock::time_point currentTime = chrono::steady_clock::now();
        file << "Elapsed time: " << (chrono::duration_cast<chrono::seconds>(currentTime - startTime).count())  << " s." << endl;
        file << "Initial volume: " << initialVolume << endl;
        file << "End volume: " << endVolume << endl;
        file << "Volume change: " << volumeChange << " %" << endl; 
        file.close();
    }

    return 0;
}