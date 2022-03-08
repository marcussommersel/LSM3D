#include "testCases.h"

Velocity vortexVelocity(int M, int N, int P, vector<double> X, vector<double> Y, vector<double> Z){

    vector<double> U;
    vector<double> V;
    vector<double> W;

    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                U.push_back(2*sin(PI*X[i])*sin(PI*X[i])*sin(2*PI*Y[j])*sin(2*PI*Z[k]));
                V.push_back(-sin(2*PI*X[i])*sin(PI*Y[j])*sin(PI*Y[j])*sin(2*PI*Z[k]));
                W.push_back(-sin(2*PI*X[i])*sin(2*PI*Y[j])*sin(PI*Z[k])*sin(PI*Z[k]));
            }
        }
    }
    return Velocity {U, V, W};

// u(x, y, z) = 2 sin2(π x) sin(2π y) sin(2π z),
// v(x, y, z) = −sin(2π x) sin2(π y) sin(2π z),
// w(x, y, z) = −sin(2π x) sin(2π y) sin2(π z)

}

Velocity simpleVelocity(int M, int N, int P){
    vector<double> U = linspace(1, 1, M*N*P);
    vector<double> V = linspace(1, 1, M*N*P);
    vector<double> W = linspace(1, 1, M*N*P);
    return Velocity {U, V, W};
}