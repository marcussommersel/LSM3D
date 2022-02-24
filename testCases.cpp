#include "testCases.h"

// auto vortexVelocity(int M, int N, int P){

// }

Velocity simpleVelocity(int M, int N, int P){
    vector<double> U = linspace(1, 1, M*N*P);
    vector<double> V = linspace(1, 1, M*N*P);
    vector<double> W = linspace(1, 1, M*N*P);
    return Velocity {U, V, W};
}