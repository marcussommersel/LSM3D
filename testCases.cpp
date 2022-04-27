#include "testCases.h"

Velocity vortexVelocity(int M, int N, int P, vector<double> X, vector<double> Y, vector<double> Z, double t, double T){

    vector<double> U;
    vector<double> V;
    vector<double> W;

    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                U.push_back(2*sin(PI*X[i])*sin(PI*X[i])*sin(2*PI*Y[j])*sin(2*PI*Z[k])*cos(PI*t/T));
                V.push_back(-sin(2*PI*X[i])*sin(PI*Y[j])*sin(PI*Y[j])*sin(2*PI*Z[k])*cos(PI*t/T));
                W.push_back(-sin(2*PI*X[i])*sin(2*PI*Y[j])*sin(PI*Z[k])*sin(PI*Z[k])*cos(PI*t/T));
            }
        }
    }
    return Velocity {U, V, W};
}

Velocity shearedSphereVelocity(int M, int N, int P, vector<double> X, vector<double> Y, vector<double> Z, double t, double T){

    vector<double> U;
    vector<double> V;
    vector<double> W;

    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                // U.push_back(sin(PI*X[i])*cos(PI*Y[j])*cos(PI*t/T)); // Morgan et al.
                // V.push_back(-cos(PI*X[i])*sin(PI*Y[j])*cos(PI*t/T));
                // W.push_back(0.0);

                U.push_back(-2*sin(PI*X[i])*sin(PI*X[i])*sin(PI*Y[j])*cos(PI*Y[j])*cos(PI*t/T)); // Rider and Kothe
                V.push_back(2*sin(PI*Y[j])*sin(PI*Y[j])*sin(PI*X[i])*cos(PI*X[i])*cos(PI*t/T));
                W.push_back(0.0);
            }
        }
    }
    return Velocity {U, V, W};
}

Velocity simpleVelocity(int M, int N, int P){
    vector<double> U = linspace(1, 1, M*N*P);
    vector<double> V = linspace(1, 1, M*N*P);
    vector<double> W = linspace(1, 1, M*N*P);
    return Velocity {U, V, W};
}

double volume(vector<double> &phi, double dx, double dy, double dz){
    double epsilon = 1.5*dx;
    double V = 0;
    for (int i = 0; i < phi.size(); ++i){
        double H;
        if (phi[i] < -epsilon){
            H = 0.0;
        } else if (-epsilon <= phi[i] && phi[i] <= epsilon){
            H = 0.5 + phi[i]/(2*epsilon) + 1/(2*PI)*sin(PI*phi[i]/epsilon);
        } else if (epsilon < phi[i]){
            H = 1.0;
        }
        V += (1-H)*dx*dy*dz;
    }
    return V;
}