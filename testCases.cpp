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
                U.push_back(sin(PI*X[i])*cos(PI*Y[j])*cos(PI*t/T));
                V.push_back(-cos(PI*X[i])*sin(PI*Y[j])*cos(PI*t/T));
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

double surfaceArea(vector<double> &phi, double dx, double dy, double dz, double M, double N, double P){
    double A = 0;
    double epsilon = 1.5*dx;
    double phix;
    double phiy;
    double phiz;
    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                if (i==0 || i==(M-1) || j==0 || j==(N-1) || k==0 || k==(P-1)){
                    continue;
                }
                double phix = (phi[(i+1)+j*N+k*P*P] - phi[(i-1)+j*N+k*P*P])/(2*dx);
                if (phix == 0){
                    phix = (phi[(i+1)+j*N+k*P*P] - phi[i+j*N+k*P*P])/(dx);
                }
                double phiy = (phi[i+(j+1)*N+k*P*P] - phi[i+(j-1)*N+k*P*P])/(2*dy);
                if (phiy == 0){
                    phiy = (phi[i+(j+1)*N+k*P*P] - phi[i+j*N+k*P*P])/(dy);
                }
                double phiz = (phi[i+j*N+(k+1)*P*P] - phi[i+j*N+(k-1)*P*P])/(2*dz);
                if (phiz == 0){
                    phiz = (phi[i+j*N+(k+1)*P*P] - phi[i+j*N+k*P*P])/(dz);
                }
                double sigma;
                if (phi[i + j*N + k*P*P] < -epsilon){
                    sigma = 0.0;
                } else if (-epsilon <= phi[i + j*N + k*P*P] && phi[i + j*N + k*P*P] <= epsilon){
                    sigma = (1/(2*epsilon) + 1/(2*epsilon)*cos(phi[i + j*N + k*P*P]*PI/epsilon));
                } else if (epsilon < phi[i + j*N + k*P*P]){
                    sigma = 0.0;
                }
                A += (sigma)*sqrt(phix*phix + phiy*phiy + phiz*phiz)*dx*dy*dz;
            }
        }
    }
    return A;
}

double interfaceError(vector<double> &phi0, vector<double> &phi, double dx, double dy, double dz, double M, double N, double P){
    double epsilon = 1.5*dx;
    double A = surfaceArea(phi0, dx, dy, dz, M, N, P);
    double L1 = 0;
    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                L1 += abs(1.0*(phi0[i + j*N + k*P*P] < 0) - 1.0*(phi[i + j*N + k*P*P] < 0))*dx*dy*dz;
            }
        }
    }
    return L1/A;
}

double massError(vector<double> &phi, double dx, double dy, double dz, double M, double N, double P){
    double error = 0;
    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                error += abs(1.0*(phi[i + j*N + k*P*P] < 0))*dx*dy*dz;
            }
        }
    }
    return error;
}
