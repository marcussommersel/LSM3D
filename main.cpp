#include <vector>
#include <iostream>
#include <ctime>
#include <chrono>
#include <fstream>
#include "initialization.h"
#include "schemes.h"
#include "vectorUtilities.h"
#include "particleLSM.h"
#include "testCases.h"
using namespace std;

int main(){ 

    // start of setup.
    chrono::steady_clock::time_point startTime = chrono::steady_clock::now();

    cout << "Program start." << endl;

    // grid nodes in each spatial direction
    const int m = 200; // x-direction
    const int n = 200; // y-direction
    const int p = 200; // z-direction

    // size of domain
    const double xStart = 0.0;
    const double xEnd = 1.0;
    const double yStart = 0.0;
    const double yEnd = 1.0;
    const double zStart = 0.0;
    const double zEnd = 1.0;

    // vector of nodes in each direction
    vector<double> x = linspace(xStart, xEnd, m);
    vector<double> y = linspace(yStart, yEnd, n);
    vector<double> z = linspace(zStart, zEnd, p);

    // grid spacing in each direction
    double dx = x[1] - x[0];
    double dy = y[1] - y[0];
    double dz = z[1] - z[0];

    // varius parameters that may be changed
    double dtau = 0.5*dx; // size of pseudo-time step
    bool doReinit = true; // true if reinitialization should be done
    bool doParticle = false; // true if particles should be used for level set method
    bool saveParticles = false; // true if particles should be saved for plotting
    int nParticles = 64; // number of particles of each type (positive and negative) in each cell
    int reinitFreq = 1; // how often reinitialization should be conducted
    int reinitSteps = 5; // how many pseudo-time steps should be done
    int plotFreq = 20000; // how often signed distance field is saved, not including first and last time step
    int reseedFreq = 100; // how often particles are reseeded
    int itmax = 20000; // maximum number of iterations
    double CFL = 0.9; // CFL-number
    bool halfplot = true; // true if signed distance field should be saved at 0.5 t/T

    // The three implemented test cases
    string testcase = "vortex";
    // string testcase = "sheared";
    // string testcase = "simple";
    string savePath = "figures/";

    Point c;
    double r;
    double T;
    if (testcase == "vortex"){
        c = Point(0.35,0.35,0.35);
        r = 0.15;
        T = 3.0;
    } else if (testcase == "sheared"){
        c = Point(0.5,0.75,0.5);
        r = 0.15;
        T = 10.0;
    } else if (testcase == "simple"){
        c = Point(0.35,0.35,0.35);
        r = 0.15;
        T = 0.4;
    }

    vector<double> phi;

    // initial signed distance field
    signedDistanceField(phi, x, y, z, r, c, m, n, p);

    // measures of error
    double initialVolume = volume(phi, dx, dy, dz); 
    vector<double> phi0 = phi;
    double MError0 = massError(phi, dx, dy, dz, m, n, p);
    double MError = 0;

    vector<string> plotTimes;
    vector<string> plotTimesParticle;
    
    // save signed distance field
    saveScalarField(savePath + to_string(0.000000) + ".txt", phi, x, y, z, m, n, p);
    plotTimes.push_back(to_string(0.000000));

    // parameters used in particle level set method
    vector<Particle> particles;
    double rmin = 0.1*min(dx, min(dy, dz));
    double rmax = 0.5*max(dx, max(dy, dz));
    double bmin = rmin;
    double bmax = 3.0*max(dx, max(dy, dz));

    // initializing particles
    cout << "Initializing particles." << endl;
    if (doParticle){
        Derivative norm = normal(phi, dx, dy, dz, m, n, p);
        for (int k = 0; k < p; ++k){
            for (int j = 0; j < n; ++j){
                for (int i = 0; i < m; ++i){
                    if (abs(phi[i + j*n + k*p*p]) < 3*max(dx, max(dy,dz))){
                        vector<Particle> newParticles = initializeParticles(x[i], y[j], z[k], dx, dy, dz, x, y, z, phi, norm, m, n, p, nParticles);
                        particles.insert(particles.end(), newParticles.begin(), newParticles.end());
                    }
                }
            }
        }

        cout << "Initialization finished." << endl;

        if (saveParticles){
            plotParticles(savePath + to_string(0.000000) + "particle.txt" , particles);
            plotTimesParticle.push_back(to_string(0.000000) + "particle");
        }
    }

    vector<double> ax;
    vector<double> ay;
    vector<double> az;

    double t = 0;

    // initial velocity field
    if (testcase == "vortex"){
        Velocity a = vortexVelocity(m, n, p, x, y, z, t, T);
        ax = a.x;
        ay = a.y;
        az = a.z;
    } else if (testcase == "sheared"){
        Velocity a = shearedSphereVelocity(m, n, p, x, y, z, t, T);
        ax = a.x;
        ay = a.y;
        az = a.z;
    } else if (testcase == "simple"){
        Velocity a = simpleVelocity(m, n, p);
        ax = a.x;
        ay = a.y;
        az = a.z;
    }

    double dtmax = CFL/(vectorMax(vectorAbs(ax)/dx + vectorAbs(ay)/dy + vectorAbs(az)/dz)); // max time step
    double dt;
    int numIt = 0;

    cout << "Setup complete." << endl;
    chrono::steady_clock::time_point currentTime = chrono::steady_clock::now();
    cout << "Elapsed time: " << (chrono::duration_cast<chrono::seconds>(currentTime - startTime).count())  << " s." << endl << endl;

    // main iteration loop
    for (int it = 1; it < itmax; ++it){

        // statement to make sure loop has one iteration at 0.5 t/T
        if (t < 0.5*T){
            dt = min(dtmax, 0.5*T - t);
        } else if (t >= 0.5*T){
            dt = min(dtmax, T - t);
        }

        t += dt;

        // advection of level set method
        TVDRK3_weno(phi, ax, ay, az, m, n, p, dx, dy, dz, dt);

        // advection of particles
        if (doParticle){
            for (int a = 0; a < particles.size(); ++a){
                
                int i = (int)(particles[a].x/dx);
                int j = (int)(particles[a].y/dy);
                int k = (int)(particles[a].z/dz);

                double Up = trilinearInterpolation(particles[a].x, particles[a].y, particles[a].z, 
                    x[i], x[i+1], y[j], y[j+1], z[k], z[k+1], 
                    ax[i+j*n+k*p*p], 
                    ax[(i+1)+j*n+k*p*p], 
                    ax[(i+1)+(j+1)*n+k*p*p], 
                    ax[i+(j+1)*n+k*p*p], 
                    ax[i+j*n+(k+1)*p*p], 
                    ax[(i+1)+j*n+(k+1)*p*p], 
                    ax[(i+1)+(j+1)*n+(k+1)*p*p], 
                    ax[i+(j+1)*n+(k+1)*p*p]);

                double Vp = trilinearInterpolation(particles[a].x, particles[a].y, particles[a].z,
                    x[i], x[i+1], y[j], y[j+1], z[k], z[k+1], 
                    ay[i+j*n+k*p*p], 
                    ay[(i+1)+j*n+k*p*p], 
                    ay[(i+1)+(j+1)*n+k*p*p], 
                    ay[i+(j+1)*n+k*p*p], 
                    ay[i+j*n+(k+1)*p*p], 
                    ay[(i+1)+j*n+(k+1)*p*p], 
                    ay[(i+1)+(j+1)*n+(k+1)*p*p], 
                    ay[i+(j+1)*n+(k+1)*p*p]);

                double Wp = trilinearInterpolation(particles[a].x, particles[a].y, particles[a].z, 
                    x[i], x[i+1], y[j], y[j+1], z[k], z[k+1], 
                    az[i+j*n+k*p*p], 
                    az[(i+1)+j*n+k*p*p], 
                    az[(i+1)+(j+1)*n+k*p*p], 
                    az[i+(j+1)*n+k*p*p], 
                    az[i+j*n+(k+1)*p*p], 
                    az[(i+1)+j*n+(k+1)*p*p], 
                    az[(i+1)+(j+1)*n+(k+1)*p*p], 
                    az[i+(j+1)*n+(k+1)*p*p]);

                particles[a].x = particles[a].x + dt*Up;
                particles[a].y = particles[a].y + dt*Vp;
                particles[a].z = particles[a].z + dt*Wp;

                i = (int)(particles[a].x/dx);
                j = (int)(particles[a].y/dy);
                k = (int)(particles[a].z/dz);

                double phip = trilinearInterpolation(particles[a].x, particles[a].y, particles[a].z, 
                    x[i], x[i+1], y[j], y[j+1], z[k], z[k+1], 
                    phi[i+j*n+k*p*p], 
                    phi[(i+1)+j*n+k*p*p], 
                    phi[(i+1)+(j+1)*n+k*p*p], 
                    phi[i+(j+1)*n+k*p*p], 
                    phi[i+j*n+(k+1)*p*p], 
                    phi[(i+1)+j*n+(k+1)*p*p], 
                    phi[(i+1)+(j+1)*n+(k+1)*p*p], 
                    phi[i+(j+1)*n+(k+1)*p*p]);

                // delete particles
                if (abs(phip) - particles[a].r > bmax){
                    particles.erase(particles.begin() + a);
                }
                // interface correction
                else if ((phip < 0 && particles[a].positive) || (phip > 0 && !particles[a].positive)  && (abs(phip) > particles[a].r)){
    
                    vector<double> phiCorrected = 
                        correctInterface(particles[a], x[i], x[i+1], y[j], y[j+1], z[k], z[k+1],
                        phi[(i)+(j)*n+(k)*p*p], 
                        phi[(i+1)+(j)*n+(k)*p*p], 
                        phi[(i+1)+(j+1)*n+(k)*p*p],
                        phi[(i)+(j+1)*n+(k)*p*p], 
                        phi[(i)+(j)*n+(k+1)*p*p], 
                        phi[(i+1)+(j)*n+(k+1)*p*p],
                        phi[(i+1)+(j+1)*n+(k+1)*p*p], 
                        phi[(i)+(j+1)*n+(k+1)*p*p], 
                        phip);

                    phi[(i)+(j)*n+(k)*p*p] = phiCorrected[0];
                    phi[(i+1)+(j)*n+(k)*p*p] = phiCorrected[1];
                    phi[(i+1)+(j+1)*n+(k)*p*p] = phiCorrected[2];
                    phi[(i)+(j+1)*n+(k)*p*p] = phiCorrected[3];
                    phi[(i)+(j)*n+(k+1)*p*p] = phiCorrected[4];
                    phi[(i+1)+(j)*n+(k+1)*p*p] = phiCorrected[5];
                    phi[(i+1)+(j+1)*n+(k+1)*p*p] = phiCorrected[6];
                    phi[(i)+(j+1)*n+(k+1)*p*p] = phiCorrected[7];
                        
                }

                phip = trilinearInterpolation(particles[a].x, particles[a].y, particles[a].z, 
                    x[i], x[i+1], y[j], y[j+1], z[k], z[k+1], 
                    phi[i+j*n+k*p*p], 
                    phi[(i+1)+j*n+k*p*p], 
                    phi[(i+1)+(j+1)*n+k*p*p], 
                    phi[i+(j+1)*n+k*p*p], 
                    phi[i+j*n+(k+1)*p*p], 
                    phi[(i+1)+j*n+(k+1)*p*p], 
                    phi[(i+1)+(j+1)*n+(k+1)*p*p], 
                    phi[i+(j+1)*n+(k+1)*p*p]);

                // adjust radius
                if (sign(phip)*phip > rmax){
                    particles[a].r = rmax;
                } else if (sign(phip)*phip < rmin){
                    particles[a].r = rmin;
                } else {
                    particles[a].r = sign(phip)*phip;
                }

            }

            // initialize new particles
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

                Derivative norm = normal(phi, dx, dy, dz, m, n, p);
                for (int a = 0; a < cellx.size(); ++a){
                    int num = nParticles - cellParticles[a];
                    if (num > 0){
                        vector<Particle> newParticles = initializeParticles(x[cellx[a]], y[celly[a]], z[cellz[a]], dx, dy, dz, x, y, z, phi, norm, m, n, p, num);
                        particles.insert(particles.end(), newParticles.begin(), newParticles.end());
                    }
                }
            }
        }

        // reinitialization loop
        if (doReinit && (it%reinitFreq == 0)){
            vector<double> phi0 = phi;
            for (int i = 0; i < reinitSteps - 1; ++i){
                // euler_godunov_reinit(phi, m, n, p, dx, dy, dz, dtau, phi0);
                TVDRK3_godunov_reinit(phi, m, n, p, dx, dy, dz, dtau, phi0);
            }
        }

        // mass error
        MError += abs(massError(phi, dx, dy, dz, m, n, p) - MError0)*dt;

        // printing to console
        if (it%10 == 0){
            cout << "Iteration: " << it << endl;
            chrono::steady_clock::time_point currentTime = chrono::steady_clock::now();
            cout << "Elapsed time: " << (chrono::duration_cast<chrono::seconds>(currentTime - startTime).count())  << " s." << endl;
            cout << "t = " << t << endl;
        }
        if (it%plotFreq == 0){
            cout << "Saving scalar field." << endl;
            saveScalarField(savePath + to_string(t) + ".txt", phi, x, y, z, m, n, p);
            plotTimes.push_back(to_string(t));

            if (doParticle && saveParticles){
                plotParticles(savePath + to_string(t) + "particle.txt" , particles);
                plotTimesParticle.push_back(to_string(t) + "particle");
            }
            cout << "Done saving." << endl;
        }

        // plotting at 0.5 t/T
        if (t/T == 0.5 && halfplot){
            cout << "Saving scalar field." << endl;
            saveScalarField(savePath + to_string(t) + ".txt", phi, x, y, z, m, n, p);
            plotTimes.push_back(to_string(t));
            cout << "Done saving." << endl;
        }

        // finish iterations if t = T
        if (t == T){
            numIt = it;
            break;
        }

        // find velocity for next time step
        if (testcase == "vortex"){
            Velocity a = vortexVelocity(m, n, p, x, y, z, t, T);
            ax = a.x;
            ay = a.y;
            az = a.z;
        } else if (testcase == "sheared"){
            Velocity a = shearedSphereVelocity(m, n, p, x, y, z, t, T);
            ax = a.x;
            ay = a.y;
            az = a.z;
        } else if (testcase == "simple"){
            Velocity a = simpleVelocity(m, n, p);
            ax = a.x;
            ay = a.y;
            az = a.z;
        }

    }

    // writing last time step to console
    {
        cout << "Iteration: " << numIt << endl;
        chrono::steady_clock::time_point currentTime = chrono::steady_clock::now();
        cout << "Elapsed time: " << (chrono::duration_cast<chrono::seconds>(currentTime - startTime).count())  << " s." << endl;
        cout << "t = " << t << endl;
    }

    // error measures
    double endVolume = volume(phi, dx, dy, dz);
    double volumeChange = 100*(endVolume-initialVolume)/initialVolume;
    double L1Error = interfaceError(phi0, phi, dx, dy, dz, m, n, p);
    MError = MError/t;

    // saving of signed distance field for final time step
    saveScalarField(savePath + to_string(T) + ".txt", phi, x, y, z, m, n, p);
    plotTimes.push_back(to_string(T));

    if (doParticle && saveParticles){
        plotParticles(savePath + to_string(t) + "particle.txt" , particles);
        plotTimesParticle.push_back(to_string(t) + "particle");
    }

    {
        ofstream file;
        file.open(savePath + "plotTimes.txt");
        if (!file.is_open()){cerr << "could not open file." << endl;}
        for (int i = 0; i < plotTimes.size(); ++i){
            file << plotTimes[i] << endl;
        }
        file.close();
    }

    if (doParticle && saveParticles){
        ofstream file;
        file.open(savePath + "plotTimesParticle.txt");
        if (!file.is_open()){cerr << "could not open file." << endl;}
        for (int i = 0; i < plotTimesParticle.size(); ++i){
            file << plotTimesParticle[i] << endl;
        }
        file.close();
    }

    // printing log to file
    {
        ofstream file;
        file.open(savePath + "log.txt");
        if (!file.is_open()){cerr << "could not open file." << endl;}
        file << "Iterations: " << numIt << endl;
        chrono::steady_clock::time_point currentTime = chrono::steady_clock::now();
        file << "Elapsed time: " << (chrono::duration_cast<chrono::seconds>(currentTime - startTime).count())  << " s." << endl;
        file << "Initial volume: " << initialVolume << endl;
        file << "End volume: " << endVolume << endl;
        file << "Volume change: " << volumeChange << " %" << endl; 
        file << "Interface error: " << L1Error << endl;
        file << "Average area error: " << MError << endl;
        file.close();
    }

    return 0;
}
