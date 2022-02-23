#include "testCases.h"

// auto vortexVelocity(int M, int N, int P){

// }

auto simpleVelocity(int M, int N, int P){
    vector<double> U = linspace(1, 1, M);
    vector<double> V = linspace(1, 1, N);
    vector<double> W = linspace(1, 1, P);
    struct result {vector<double> x; vector<double> y; vector<double> z;};
    return result {U, V, W};
}