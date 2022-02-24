#include "schemes.h"

Derivative upwind(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int const N, int const P, double dx, double dy, double dz){

    // const size_t size = M*N*P;
    vector<double> phix;
    vector<double> phiy;
    vector<double> phiz;

    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                
                if (i==0 || i==M-1 || j==0 || j==N-1 || k==0 || k==P-1){
                    phix.push_back(0);
                    phiy.push_back(0);
                    phiz.push_back(0);
                    continue;
                }

                if (AX[i + j*N + k*P*P] >= 0){
                    phix.push_back((arr[i + j*N + k*P*P] - arr[(i - 1) + j*N + k*P*P])/dx);
                } else if (AX[i + j*N + k*P*P] < 0){
                    phix.push_back((arr[(i + 1) + j*N + k*P*P] - arr[i + j*N + k*P*P])/dx);
                }

                if (AY[i + j*N + k*P*P] >= 0){
                    phiy.push_back((arr[i + j*N + k*P*P] - arr[i + (j - 1)*N + k*P*P])/dy);
                } else if (AY[i + j*N + k*P*P] < 0){
                    phiy.push_back((arr[i + (j + 1)*N + k*P*P] - arr[i + j*N + k*P*P])/dy);
                }

                if (AZ[i + j*N + k*P*P] >= 0){
                    phiz.push_back((arr[i + j*N + k*P*P] - arr[i + j*N + (k - 1)*P*P])/dz);
                } else if (AX[i + j*N + k*P*P] < 0){
                    phiz.push_back((arr[i + j*N + (k + 1)*P*P] - arr[i + j*N + k*P*P])/dz);
                }
            }
        }
    }
    cout << phix.size() << " " << phiy.size() << " " << phiz.size() << endl;
    return Derivative{phix, phiy, phiz};
}

void euler_upwind(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt){
    
    auto [phix, phiy, phiz] = upwind(arr, AX, AY, AZ, M, N, P, dx, dy, dz);
    arr = arr - dt*(AX*phix + AY*phiy + AZ*phiz);
}

