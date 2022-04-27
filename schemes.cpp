#include "schemes.h"

Derivative upwind(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz){

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
                } else if (AZ[i + j*N + k*P*P] < 0){
                    phiz.push_back((arr[i + j*N + (k + 1)*P*P] - arr[i + j*N + k*P*P])/dz);
                }
            }
        }
    }
    cout << AX.size() << " " << AY.size() << " " << AZ.size() << endl;
    return Derivative{phix, phiy, phiz};
}

Derivative upwind2(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz){

    vector<double> phix;
    vector<double> phiy;
    vector<double> phiz;

    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                
                if (i==0 || i==M-1 || i==1 || i==M-2 || j==0 || j==N-1 || j==1 || j==N-2 || k==0 || k==P-1 || k==1 || k==P-2){
                    phix.push_back(0);
                    phiy.push_back(0);
                    phiz.push_back(0);
                    continue;
                }

                if (AX[i+j*N+k*P*P] >= 0){
                    phix.push_back((3*arr[i+j*N+k*P*P] - 4*arr[i-1+j*N+k*P*P] + arr[i-2+j*N+k*P*P])/(2*dx));
                } else if (AX[i+j*N+k*P*P] < 0) {
                    phix.push_back((-arr[i+2+j*N+k*P*P] + 4*arr[i+1+j*N+k*P*P] - 3*arr[i+j*N+k*P*P])/(2*dx));
                }
                if (AY[i+j*N+k*P*P] >= 0){
                    phiy.push_back((3*arr[i+j*N+k*P*P] - 4*arr[i+(j-1)*N+k*P*P] + arr[i+(j-2)*N+k*P*P])/(2*dy));
                } else if (AY[i+j*N+k*P*P] < 0) {
                    phiy.push_back((-arr[i+(j+2)*N+k*P*P] + 4*arr[i+(j+1)*N+k*P*P] - 3*arr[i+j*N+k*P*P])/(2*dy));
                }
                if (AZ[i+j*N+k*P*P] >= 0){
                    phiz.push_back((3*arr[i+j*N+k*P*P] - 4*arr[i+j*N+(k-1)*P*P] + arr[i+j*N+(k-2)*P*P])/(2*dz));
                } else if (AZ[i+j*N+k*P*P] < 0) {
                    phiz.push_back((-arr[i+j*N+(k+2)*P*P] + 4*arr[i+j*N+(k+1)*P*P] - 3*arr[i+j*N+k*P*P])/(2*dz));
                }
            }
        }
    } 
    return Derivative{phix, phiy, phiz};
}

Derivative weno(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int const N, int const P, double dx, double dy, double dz){

    vector<double> phix;
    vector<double> phiy;
    vector<double> phiz;

    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                if (i == 0 || i == 1 || i==2 || i==(M-3) || i == (M-2) || i == (M-1)){
                    phix.push_back(0);
                    phiy.push_back(0);
                    phiz.push_back(0);
                    continue;
                }
                if (j == 0 || j == 1 || j==2 || j==(N-3) || j == (N-2) || j == (N-1)){
                    phix.push_back(0);
                    phiy.push_back(0);
                    phiz.push_back(0);
                    continue;
                }
                if (k == 0 || k == 1 || k==2 || k==(P-3) || k == (P-2) || k == (P-1)){
                    phix.push_back(0);
                    phiy.push_back(0);
                    phiz.push_back(0);
                    continue;
                }

                double v1;
                double v2;
                double v3;
                double v4;
                double v5;                

                {
                if (AX[i + j*N + k*P*P] >= 0){
                    v1 = (arr[(i-2) + j*N + k*P*P] - arr[(i-3) + j*N + k*P*P])/dx;
                    v2 = (arr[(i-1) + j*N + k*P*P] - arr[(i-2) + j*N + k*P*P])/dx;
                    v3 = (arr[(i) + j*N + k*P*P] - arr[(i-1) + j*N + k*P*P])/dx;
                    v4 = (arr[(i+1) + j*N + k*P*P] - arr[(i) + j*N + k*P*P])/dx;
                    v5 = (arr[(i+2) + j*N + k*P*P] - arr[(i+1) + j*N + k*P*P])/dx;
                } else if (AX[i + j*N + k*P*P] < 0){
                    v1 = (arr[(i-1) + j*N + k*P*P] - arr[(i-2) + j*N + k*P*P])/dx;
                    v2 = (arr[(i) + j*N + k*P*P] - arr[(i-1) + j*N + k*P*P])/dx;
                    v3 = (arr[(i+1) + j*N + k*P*P] - arr[(i) + j*N + k*P*P])/dx;
                    v4 = (arr[(i+2) + j*N + k*P*P] - arr[(i+1) + j*N + k*P*P])/dx;
                    v5 = (arr[(i+3) + j*N + k*P*P] - arr[(i+2) + j*N + k*P*P])/dx;
                }
                double S1 = 13/12*(v1 - 2*v2 + v3)*(v1 - 2*v2 + v3) + 1/4*(v1 - 4*v2 + v3)*(v1 - 4*v2 + v3);
                double S2 = 13/12*(v2 - 2*v3 + v4)*(v2 - 2*v3 + v4) + 1/4*(v2 - v4)*(v2 - v4);
                double S3 = 13/12*(v3 - 2*v4 + v5)*(v3 - 2*v4 + v5) + 1/4*(3*v3 - 4*v4 + v5)*(3*v3 - 4*v4 + v5);

                double epsilon = pow(10, -6)*max(max(max(max(v1*v1, v2*v2), v3*v3), v4*v4), v5*v5) + pow(10, -99);

                double alpha1 = 0.1/((S1 + epsilon)*(S1 + epsilon));
                double alpha2 = 0.6/((S2 + epsilon)*(S2 + epsilon));
                double alpha3 = 0.3/((S3 + epsilon)*(S3 + epsilon));

                double omega1 = alpha1/(alpha1 + alpha2 + alpha3);
                double omega2 = alpha2/(alpha1 + alpha2 + alpha3);
                double omega3 = alpha3/(alpha1 + alpha2 + alpha3);

                double phix1 = v1/3 - 7*v2/6 + 11*v3/6;
                double phix2 = -v2/6 + 5*v3/6 + v4/3;
                double phix3 = v3/3 + 5*v4/6 - v5/6;

                phix.push_back(omega1*phix1 + omega2*phix2 + omega3*phix3);
                }

                {
                if (AY[i + j*N + k*P*P] >= 0){
                    v1 = (arr[i + (j-2)*N + k*P*P] - arr[i + (j-3)*N + k*P*P])/dy;
                    v2 = (arr[i + (j-1)*N + k*P*P] - arr[i + (j-2)*N + k*P*P])/dy;
                    v3 = (arr[i + j*N + k*P*P] - arr[i + (j-1)*N + k*P*P])/dy;
                    v4 = (arr[i + (j+1)*N + k*P*P] - arr[i + j*N + k*P*P])/dy;
                    v5 = (arr[i + (j+2)*N + k*P*P] - arr[i + (j+1)*N + k*P*P])/dy;
                } else if (AY[i + j*N + k*P*P] < 0){
                    v1 = (arr[i + (j-1)*N + k*P*P] - arr[i + (j-2)*N + k*P*P])/dy;
                    v2 = (arr[i + j*N + k*P*P] - arr[i + (j-1)*N + k*P*P])/dy;
                    v3 = (arr[i + (j+1)*N + k*P*P] - arr[i + j*N + k*P*P])/dy;
                    v4 = (arr[i + (j+2)*N + k*P*P] - arr[i + (j+1)*N + k*P*P])/dy;
                    v5 = (arr[i + (j+3)*N + k*P*P] - arr[i + (j+2)*N + k*P*P])/dy;
                }
                double S1 = 13/12*(v1 - 2*v2 + v3)*(v1 - 2*v2 + v3) + 1/4*(v1 - 4*v2 + v3)*(v1 - 4*v2 + v3);
                double S2 = 13/12*(v2 - 2*v3 + v4)*(v2 - 2*v3 + v4) + 1/4*(v2 - v4)*(v2 - v4);
                double S3 = 13/12*(v3 - 2*v4 + v5)*(v3 - 2*v4 + v5) + 1/4*(3*v3 - 4*v4 + v5)*(3*v3 - 4*v4 + v5);

                double epsilon = pow(10, -6)*max(max(max(max(v1*v1, v2*v2), v3*v3), v4*v4), v5*v5) + pow(10, -99);

                double alpha1 = 0.1/((S1 + epsilon)*(S1 + epsilon));
                double alpha2 = 0.6/((S2 + epsilon)*(S2 + epsilon));
                double alpha3 = 0.3/((S3 + epsilon)*(S3 + epsilon));

                double omega1 = alpha1/(alpha1 + alpha2 + alpha3);
                double omega2 = alpha2/(alpha1 + alpha2 + alpha3);
                double omega3 = alpha3/(alpha1 + alpha2 + alpha3);

                double phiy1 = v1/3 - 7*v2/6 + 11*v3/6;
                double phiy2 = -v2/6 + 5*v3/6 + v4/3;
                double phiy3 = v3/3 + 5*v4/6 - v5/6;

                phiy.push_back(omega1*phiy1 + omega2*phiy2 + omega3*phiy3);

                }

                {
                if (AZ[i + j*N + k*P*P] >= 0){
                    v1 = (arr[i + j*N + (k-2)*P*P] - arr[i + j*N + (k-3)*P*P])/dz;
                    v2 = (arr[i + j*N + (k-1)*P*P] - arr[i + j*N + (k-2)*P*P])/dz;
                    v3 = (arr[i + j*N + k*P*P] - arr[i + j*N + (k-1)*P*P])/dz;
                    v4 = (arr[i + j*N + (k+1)*P*P] - arr[i + j*N + k*P*P])/dz;
                    v5 = (arr[i + j*N + (k+2)*P*P] - arr[i + j*N + (k+1)*P*P])/dz;
                } else if (AZ[i + j*N + k*P*P] < 0){
                    v1 = (arr[i + j*N + (k-1)*P*P] - arr[i + j*N + (k-2)*P*P])/dz;
                    v2 = (arr[i + j*N + k*P*P] - arr[i + j*N + (k-1)*P*P])/dz;
                    v3 = (arr[i + j*N + (k+1)*P*P] - arr[i + j*N + k*P*P])/dz;
                    v4 = (arr[i + j*N + (k+2)*P*P] - arr[i + j*N + (k+1)*P*P])/dz;
                    v5 = (arr[i + j*N + (k+3)*P*P] - arr[i + j*N + (k+2)*P*P])/dz;
                }
                double S1 = 13/12*(v1 - 2*v2 + v3)*(v1 - 2*v2 + v3) + 1/4*(v1 - 4*v2 + v3)*(v1 - 4*v2 + v3);
                double S2 = 13/12*(v2 - 2*v3 + v4)*(v2 - 2*v3 + v4) + 1/4*(v2 - v4)*(v2 - v4);
                double S3 = 13/12*(v3 - 2*v4 + v5)*(v3 - 2*v4 + v5) + 1/4*(3*v3 - 4*v4 + v5)*(3*v3 - 4*v4 + v5);

                double epsilon = pow(10, -6)*max(max(max(max(v1*v1, v2*v2), v3*v3), v4*v4), v5*v5) + pow(10, -99);

                double alpha1 = 0.1/((S1 + epsilon)*(S1 + epsilon));
                double alpha2 = 0.6/((S2 + epsilon)*(S2 + epsilon));
                double alpha3 = 0.3/((S3 + epsilon)*(S3 + epsilon));

                double omega1 = alpha1/(alpha1 + alpha2 + alpha3);
                double omega2 = alpha2/(alpha1 + alpha2 + alpha3);
                double omega3 = alpha3/(alpha1 + alpha2 + alpha3);

                double phiz1 = v1/3 - 7*v2/6 + 11*v3/6;
                double phiz2 = -v2/6 + 5*v3/6 + v4/3;
                double phiz3 = v3/3 + 5*v4/6 - v5/6;

                phiz.push_back(omega1*phiz1 + omega2*phiz2 + omega3*phiz3);
                }
            }
        }
    }
    return Derivative{phix, phiy, phiz};
}

Derivative godunov(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int const N, int const P, double dx, double dy, double dz){

    vector<double> phix;
    vector<double> phiy;
    vector<double> phiz;

    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){

                if (i == 0 || i == (M-1) || j == 0 || j == (N-1) || k == 0 || k == (P-1) ){
                    phix.push_back(1.0);
                    phiy.push_back(1.0);
                    phiz.push_back(1.0);
                    continue;
                }

                double phix_m = (arr[i + j*N + k*P*P] - arr[(i-1) + j*N + k*P*P])/dx;
                double phix_p = (arr[(i+1) + j*N + k*P*P] - arr[i + j*N + k*P*P])/dx;

                double phiy_m = (arr[i + j*N + k*P*P] - arr[i + (j-1)*N + k*P*P])/dy;
                double phiy_p = (arr[i + (j+1)*N + k*P*P] - arr[i + j*N + k*P*P])/dy;

                double phiz_m = (arr[i + j*N + k*P*P] - arr[i + j*N + (k-1)*P*P])/dz;
                double phiz_p = (arr[i + j*N + (k+1)*P*P] - arr[i + j*N + k*P*P])/dz;

                if (AX[i + j*N + k*P*P] >= 0){
                    phix.push_back(sqrt(max(max(phix_m, 0.0)*max(phix_m, 0.0), min(phix_p, 0.0)*min(phix_p, 0.0))));
                } else if (AX[i + j*N + k*P*P] < 0){
                    phix.push_back(sqrt(max(min(phix_m, 0.0)*min(phix_m, 0.0), max(phix_p, 0.0)*max(phix_p, 0.0))));
                }

                if (AY[i + j*N + k*P*P] >= 0){
                    phiy.push_back(sqrt(max(max(phiy_m, 0.0)*max(phiy_m, 0.0), min(phiy_p, 0.0)*min(phiy_p, 0.0))));
                } else if (AY[i + j*N + k*P*P] < 0){
                    phiy.push_back(sqrt(max(min(phiy_m, 0.0)*min(phiy_m, 0.0), max(phiy_p, 0.0)*max(phiy_p, 0.0))));
                }

                if (AZ[i + j*N + k*P*P] >= 0){
                    phiz.push_back(sqrt(max(max(phiz_m, 0.0)*max(phiz_m, 0.0), min(phiz_p, 0.0)*min(phiz_p, 0.0))));
                } else if (AZ[i + j*N + k*P*P] < 0){
                    phiz.push_back(sqrt(max(min(phiz_m, 0.0)*min(phiz_m, 0.0), max(phiz_p, 0.0)*max(phiz_p, 0.0))));
                }
            }
        }
    }
    return Derivative{phix, phiy, phiz};
}

void euler_upwind(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt){
    auto [phix, phiy, phiz] = upwind(arr, AX, AY, AZ, M, N, P, dx, dy, dz);
    arr = arr - dt*(AX*phix + AY*phiy + AZ*phiz);
}

void TVDRK3_upwind(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt){

    vector<double> n1;
    vector<double> n2;
    vector<double> n3_2;

    {
        auto [phix, phiy, phiz] = upwind(arr, AX, AY, AZ, M, N, P, dx, dy, dz);
        n1 =  arr - dt*(AX*phix + AY*phiy + AZ*phiz);
    }

    {
        auto [phix, phiy, phiz] = upwind(n1, AX, AY, AZ, M, N, P, dx, dy, dz);
        n2 =  n1 - dt*(AX*phix + AY*phiy + AZ*phiz);
    }

    vector<double> n1_2 = 3/4*arr + 1/4*n2;

    {
        auto [phix, phiy, phiz] = upwind(n1_2, AX, AY, AZ, M, N, P, dx, dy, dz);
        n3_2 =  n1_2 - dt*(AX*phix + AY*phiy + AZ*phiz);
    }

    arr = 1/3*arr + 2/3*n3_2;

}

void TVDRK3_weno(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt){

    {
    vector<double> n1;
    vector<double> n2;
    vector<double> n3_2;

    {
        auto [phix, phiy, phiz] = weno(arr, AX, AY, AZ, M, N, P, dx, dy, dz);
        n1 =  arr - dt*(AX*phix + AY*phiy + AZ*phiz);
    }

    {
        auto [phix, phiy, phiz] = weno(n1, AX, AY, AZ, M, N, P, dx, dy, dz);
        n2 =  n1 - dt*(AX*phix + AY*phiy + AZ*phiz);
        // vector<double>().swap(n1);
    }

    vector<double> n1_2 = 3.0/4*arr + 1.0/4*n2;
    // vector<double>().swap(n2);

    {
        auto [phix, phiy, phiz] = weno(n1_2, AX, AY, AZ, M, N, P, dx, dy, dz);
        n3_2 =  n1_2 - dt*(AX*phix + AY*phiy + AZ*phiz);
        // vector<double>().swap(n1_2);
    }

    arr = 1.0/3*arr + 2.0/3*n3_2;
    
    }
    
    // vector<double>().swap(n3_2);

    for(int i = 0; i < M; ++i){
        for (int j = 0; j < N; ++j){
            // arr[2 + i*N + j*P*P] = arr[3 + i*N + j*P*P] - (arr[4 + i*N + j*P*P] - arr[3 + i*N + j*P*P]);
            // arr[1 + i*N + j*P*P] = arr[2 + i*N + j*P*P] - (arr[3 + i*N + j*P*P] - arr[2 + i*N + j*P*P]);
            // arr[0 + i*N + j*P*P] = arr[1 + i*N + j*P*P] - (arr[2 + i*N + j*P*P] - arr[1 + i*N + j*P*P]);
            // arr[(M-3) + i*N + j*P*P] = arr[(M-4) + i*N + j*P*P] - (arr[(M-5) + i*N + j*P*P] - arr[(M-4) + i*N + j*P*P]);
            // arr[(M-2) + i*N + j*P*P] = arr[(M-3) + i*N + j*P*P] - (arr[(M-4) + i*N + j*P*P] - arr[(M-3) + i*N + j*P*P]);
            // arr[(M-1) + i*N + j*P*P] = arr[(M-2) + i*N + j*P*P] - (arr[(M-3) + i*N + j*P*P] - arr[(M-2) + i*N + j*P*P]);

            // arr[j + (2)*N + i*P*P] = arr[j + (3)*N + i*P*P] - (arr[j + (4)*N + i*P*P] - arr[j + (3)*N + i*P*P]);
            // arr[j + (1)*N + i*P*P] = arr[j + (2)*N + i*P*P] - (arr[j + (3)*N + i*P*P] - arr[j + (2)*N + i*P*P]);
            // arr[j + (0)*N + i*P*P] = arr[j + (1)*N + i*P*P] - (arr[j + (2)*N + i*P*P] - arr[j + (1)*N + i*P*P]);
            // arr[j + (N-3)*N + i*P*P] = arr[j + (N-4)*N + i*P*P] - (arr[j + (N-5)*N + i*P*P] - arr[j + (N-4)*N + i*P*P]);
            // arr[j + (N-2)*N + i*P*P] = arr[j + (N-3)*N + i*P*P] - (arr[j + (N-4)*N + i*P*P] - arr[j + (N-3)*N + i*P*P]);
            // arr[j + (N-1)*N + i*P*P] = arr[j + (N-2)*N + i*P*P] - (arr[j + (N-3)*N + i*P*P] - arr[j + (N-2)*N + i*P*P]);

            // arr[j + i*N + (2)*P*P] = arr[j + i*N + (3)*P*P] - (arr[j + i*N + (4)*P*P] - arr[j + i*N + (3)*P*P]);
            // arr[j + i*N + (1)*P*P] = arr[j + i*N + (2)*P*P] - (arr[j + i*N + (3)*P*P] - arr[j + i*N + (2)*P*P]);
            // arr[j + i*N + (0)*P*P] = arr[j + i*N + (1)*P*P] - (arr[j + i*N + (2)*P*P] - arr[j + i*N + (1)*P*P]);
            // arr[j + i*N + (N-3)*P*P] = arr[j + i*N + (N-4)*P*P] - (arr[j + i*N + (N-5)*P*P] - arr[j + i*N + (N-4)*P*P]);
            // arr[j + i*N + (N-2)*P*P] = arr[j + i*N + (N-3)*P*P] - (arr[j + i*N + (N-4)*P*P] - arr[j + i*N + (N-3)*P*P]);
            // arr[j + i*N + (N-1)*P*P] = arr[j + i*N + (N-2)*P*P] - (arr[j + i*N + (N-3)*P*P] - arr[j + i*N + (N-2)*P*P]);

            arr[2 + i*N + j*P*P] = arr[3 + i*N + j*P*P] + (arr[4 + i*N + j*P*P] - arr[3 + i*N + j*P*P]);
            arr[1 + i*N + j*P*P] = arr[2 + i*N + j*P*P] + (arr[3 + i*N + j*P*P] - arr[2 + i*N + j*P*P]);
            arr[0 + i*N + j*P*P] = arr[1 + i*N + j*P*P] + (arr[2 + i*N + j*P*P] - arr[1 + i*N + j*P*P]);
            arr[(M-3) + i*N + j*P*P] = arr[(M-4) + i*N + j*P*P] + (arr[(M-5) + i*N + j*P*P] - arr[(M-4) + i*N + j*P*P]);
            arr[(M-2) + i*N + j*P*P] = arr[(M-3) + i*N + j*P*P] + (arr[(M-4) + i*N + j*P*P] - arr[(M-3) + i*N + j*P*P]);
            arr[(M-1) + i*N + j*P*P] = arr[(M-2) + i*N + j*P*P] + (arr[(M-3) + i*N + j*P*P] - arr[(M-2) + i*N + j*P*P]);

            arr[j + (2)*N + i*P*P] = arr[j + (3)*N + i*P*P] + (arr[j + (4)*N + i*P*P] - arr[j + (3)*N + i*P*P]);
            arr[j + (1)*N + i*P*P] = arr[j + (2)*N + i*P*P] + (arr[j + (3)*N + i*P*P] - arr[j + (2)*N + i*P*P]);
            arr[j + (0)*N + i*P*P] = arr[j + (1)*N + i*P*P] + (arr[j + (2)*N + i*P*P] - arr[j + (1)*N + i*P*P]);
            arr[j + (N-3)*N + i*P*P] = arr[j + (N-4)*N + i*P*P] + (arr[j + (N-5)*N + i*P*P] - arr[j + (N-4)*N + i*P*P]);
            arr[j + (N-2)*N + i*P*P] = arr[j + (N-3)*N + i*P*P] + (arr[j + (N-4)*N + i*P*P] - arr[j + (N-3)*N + i*P*P]);
            arr[j + (N-1)*N + i*P*P] = arr[j + (N-2)*N + i*P*P] + (arr[j + (N-3)*N + i*P*P] - arr[j + (N-2)*N + i*P*P]);

            arr[j + i*N + (2)*P*P] = arr[j + i*N + (3)*P*P] + (arr[j + i*N + (4)*P*P] - arr[j + i*N + (3)*P*P]);
            arr[j + i*N + (1)*P*P] = arr[j + i*N + (2)*P*P] + (arr[j + i*N + (3)*P*P] - arr[j + i*N + (2)*P*P]);
            arr[j + i*N + (0)*P*P] = arr[j + i*N + (1)*P*P] + (arr[j + i*N + (2)*P*P] - arr[j + i*N + (1)*P*P]);
            arr[j + i*N + (N-3)*P*P] = arr[j + i*N + (N-4)*P*P] + (arr[j + i*N + (N-5)*P*P] - arr[j + i*N + (N-4)*P*P]);
            arr[j + i*N + (N-2)*P*P] = arr[j + i*N + (N-3)*P*P] + (arr[j + i*N + (N-4)*P*P] - arr[j + i*N + (N-3)*P*P]);
            arr[j + i*N + (N-1)*P*P] = arr[j + i*N + (N-2)*P*P] + (arr[j + i*N + (N-3)*P*P] - arr[j + i*N + (N-2)*P*P]);

        }
    }
}

void euler_weno(vector<double> &arr, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt){
    auto [phix, phiy, phiz] = weno(arr, AX, AY, AZ, M, N, P, dx, dy, dz);
    arr = arr - dt*(AX*phix + AY*phiy + AZ*phiz);

    for(int i = 0; i < M; ++i){
        for (int j = 0; j < N; ++j){
            arr[2 + i*N + j*P*P] = arr[3 + i*N + j*P*P] - (arr[4 + i*N + j*P*P] - arr[3 + i*N + j*P*P]);
            arr[1 + i*N + j*P*P] = arr[2 + i*N + j*P*P] - (arr[3 + i*N + j*P*P] - arr[2 + i*N + j*P*P]);
            arr[0 + i*N + j*P*P] = arr[1 + i*N + j*P*P] - (arr[2 + i*N + j*P*P] - arr[1 + i*N + j*P*P]);
            arr[(M-3) + i*N + j*P*P] = arr[(M-4) + i*N - j*P*P] + (arr[(M-5) + i*N + j*P*P] - arr[(M-4) + i*N + j*P*P]);
            arr[(M-2) + i*N + j*P*P] = arr[(M-3) + i*N - j*P*P] + (arr[(M-4) + i*N + j*P*P] - arr[(M-3) + i*N + j*P*P]);
            arr[(M-1) + i*N + j*P*P] = arr[(M-2) + i*N - j*P*P] + (arr[(M-3) + i*N + j*P*P] - arr[(M-2) + i*N + j*P*P]);

            arr[j + (2)*N + i*P*P] = arr[j + (3)*N + i*P*P] - (arr[j + (4)*N + i*P*P] - arr[j + (3)*N + i*P*P]);
            arr[j + (1)*N + i*P*P] = arr[j + (2)*N + i*P*P] - (arr[j + (3)*N + i*P*P] - arr[j + (2)*N + i*P*P]);
            arr[j + (0)*N + i*P*P] = arr[j + (1)*N + i*P*P] - (arr[j + (2)*N + i*P*P] - arr[j + (1)*N + i*P*P]);
            arr[j + (N-3)*N + i*P*P] = arr[j + (N-4)*N + i*P*P] - (arr[j + (N-5)*N + i*P*P] - arr[j + (N-4)*N + i*P*P]);
            arr[j + (N-2)*N + i*P*P] = arr[j + (N-3)*N + i*P*P] - (arr[j + (N-4)*N + i*P*P] - arr[j + (N-3)*N + i*P*P]);
            arr[j + (N-1)*N + i*P*P] = arr[j + (N-2)*N + i*P*P] - (arr[j + (N-3)*N + i*P*P] - arr[j + (N-2)*N + i*P*P]);

            arr[j + i*N + (2)*P*P] = arr[j + i*N + (3)*P*P] - (arr[j + i*N + (4)*P*P] - arr[j + i*N + (3)*P*P]);
            arr[j + i*N + (1)*P*P] = arr[j + i*N + (2)*P*P] - (arr[j + i*N + (3)*P*P] - arr[j + i*N + (2)*P*P]);
            arr[j + i*N + (0)*P*P] = arr[j + i*N + (1)*P*P] - (arr[j + i*N + (2)*P*P] - arr[j + i*N + (1)*P*P]);
            arr[j + i*N + (N-3)*P*P] = arr[j + i*N + (N-4)*P*P] - (arr[j + i*N + (N-5)*P*P] - arr[j + i*N + (N-4)*P*P]);
            arr[j + i*N + (N-2)*P*P] = arr[j + i*N + (N-3)*P*P] - (arr[j + i*N + (N-4)*P*P] - arr[j + i*N + (N-3)*P*P]);
            arr[j + i*N + (N-1)*P*P] = arr[j + i*N + (N-2)*P*P] - (arr[j + i*N + (N-3)*P*P] - arr[j + i*N + (N-2)*P*P]);
        }
    }
}

void TVDRK3_godunov_reinit(vector<double> &arr, vector<double> X, vector<double> Y, vector<double> Z, int M, int N, int P, double dx, double dy, double dz, double dt, vector<double> phi0){

    vector<double> n1;
    vector<double> n2;
    vector<double> n3_2;

    vector<double> S = phi0/(vectorSqrt(phi0*phi0 + max(max(dx,dy),dz)*max(max(dx,dy),dz)));

    {
        auto [phix, phiy, phiz] = godunov(arr, S, S, S, M, N, P, dx, dy, dz);
        // n1 =  arr - dt*S*(vectorCBRT(phix*phix*phix + phiy*phiy*phiy + phiz*phiz*phiz) - 1.0);// feil her?
        n1 =  arr - dt*S*(vectorSqrt(phix*phix + phiy*phiy + phiz*phiz) - 1.0);// feil her?
    }

    {
        auto [phix, phiy, phiz] = godunov(n1, S, S, S, M, N, P, dx, dy, dz);
        // n2 =  n1 - dt*S*(vectorCBRT(phix*phix*phix + phiy*phiy*phiy + phiz*phiz*phiz) - 1.0);// feil her?
        n2 =  n1 - dt*S*(vectorSqrt(phix*phix + phiy*phiy + phiz*phiz) - 1.0);// feil her? skal være ^2?
    }

    vector<double> n1_2 = (3.0/4*arr + 1.0/4*n2);

    {
        auto [phix, phiy, phiz] = godunov(n1_2, S, S, S, M, N, P, dx, dy, dz);
        // n3_2 =  n1_2 - dt*S*(vectorCBRT(phix*phix*phix + phiy*phiy*phiy + phiz*phiz*phiz) - 1.0); // feil her?
        n3_2 =  n1_2 - dt*S*(vectorSqrt(phix*phix + phiy*phiy + phiz*phiz) - 1.0); // feil her?
    }

    arr = 1.0/3*arr + 2.0/3*n3_2;

    for(int i = 0; i < M; ++i){
        for (int j = 0; j < N; ++j){ // feil under
            arr[2 + i*N + j*P*P] = arr[3 + i*N + j*P*P] + (arr[4 + i*N + j*P*P] - arr[3 + i*N + j*P*P]);
            arr[1 + i*N + j*P*P] = arr[2 + i*N + j*P*P] + (arr[3 + i*N + j*P*P] - arr[2 + i*N + j*P*P]);
            arr[0 + i*N + j*P*P] = arr[1 + i*N + j*P*P] + (arr[2 + i*N + j*P*P] - arr[1 + i*N + j*P*P]);
            arr[(M-3) + i*N + j*P*P] = arr[(M-4) + i*N + j*P*P] + (arr[(M-5) + i*N + j*P*P] - arr[(M-4) + i*N + j*P*P]);
            arr[(M-2) + i*N + j*P*P] = arr[(M-3) + i*N + j*P*P] + (arr[(M-4) + i*N + j*P*P] - arr[(M-3) + i*N + j*P*P]);
            arr[(M-1) + i*N + j*P*P] = arr[(M-2) + i*N + j*P*P] + (arr[(M-3) + i*N + j*P*P] - arr[(M-2) + i*N + j*P*P]);

            arr[j + (2)*N + i*P*P] = arr[j + (3)*N + i*P*P] + (arr[j + (4)*N + i*P*P] - arr[j + (3)*N + i*P*P]);
            arr[j + (1)*N + i*P*P] = arr[j + (2)*N + i*P*P] + (arr[j + (3)*N + i*P*P] - arr[j + (2)*N + i*P*P]);
            arr[j + (0)*N + i*P*P] = arr[j + (1)*N + i*P*P] + (arr[j + (2)*N + i*P*P] - arr[j + (1)*N + i*P*P]);
            arr[j + (N-3)*N + i*P*P] = arr[j + (N-4)*N + i*P*P] + (arr[j + (N-5)*N + i*P*P] - arr[j + (N-4)*N + i*P*P]);
            arr[j + (N-2)*N + i*P*P] = arr[j + (N-3)*N + i*P*P] + (arr[j + (N-4)*N + i*P*P] - arr[j + (N-3)*N + i*P*P]);
            arr[j + (N-1)*N + i*P*P] = arr[j + (N-2)*N + i*P*P] + (arr[j + (N-3)*N + i*P*P] - arr[j + (N-2)*N + i*P*P]);

            arr[j + i*N + (2)*P*P] = arr[j + i*N + (3)*P*P] + (arr[j + i*N + (4)*P*P] - arr[j + i*N + (3)*P*P]);
            arr[j + i*N + (1)*P*P] = arr[j + i*N + (2)*P*P] + (arr[j + i*N + (3)*P*P] - arr[j + i*N + (2)*P*P]);
            arr[j + i*N + (0)*P*P] = arr[j + i*N + (1)*P*P] + (arr[j + i*N + (2)*P*P] - arr[j + i*N + (1)*P*P]);
            arr[j + i*N + (N-3)*P*P] = arr[j + i*N + (N-4)*P*P] + (arr[j + i*N + (N-5)*P*P] - arr[j + i*N + (N-4)*P*P]);
            arr[j + i*N + (N-2)*P*P] = arr[j + i*N + (N-3)*P*P] + (arr[j + i*N + (N-4)*P*P] - arr[j + i*N + (N-3)*P*P]);
            arr[j + i*N + (N-1)*P*P] = arr[j + i*N + (N-2)*P*P] + (arr[j + i*N + (N-3)*P*P] - arr[j + i*N + (N-2)*P*P]);


            // arr[2 + i*N + j*P*P] = arr[3 + i*N + j*P*P] - (arr[4 + i*N + j*P*P] - arr[3 + i*N + j*P*P]);
            // arr[1 + i*N + j*P*P] = arr[2 + i*N + j*P*P] - (arr[3 + i*N + j*P*P] - arr[2 + i*N + j*P*P]);
            // arr[0 + i*N + j*P*P] = arr[1 + i*N + j*P*P] - (arr[2 + i*N + j*P*P] - arr[1 + i*N + j*P*P]);
            // arr[(M-3) + i*N + j*P*P] = arr[(M-4) + i*N + j*P*P] - (arr[(M-5) + i*N + j*P*P] - arr[(M-4) + i*N + j*P*P]);
            // arr[(M-2) + i*N + j*P*P] = arr[(M-3) + i*N + j*P*P] - (arr[(M-4) + i*N + j*P*P] - arr[(M-3) + i*N + j*P*P]);
            // arr[(M-1) + i*N + j*P*P] = arr[(M-2) + i*N + j*P*P] - (arr[(M-3) + i*N + j*P*P] - arr[(M-2) + i*N + j*P*P]);

            // arr[j + (2)*N + i*P*P] = arr[j + (3)*N + i*P*P] - (arr[j + (4)*N + i*P*P] - arr[j + (3)*N + i*P*P]);
            // arr[j + (1)*N + i*P*P] = arr[j + (2)*N + i*P*P] - (arr[j + (3)*N + i*P*P] - arr[j + (2)*N + i*P*P]);
            // arr[j + (0)*N + i*P*P] = arr[j + (1)*N + i*P*P] - (arr[j + (2)*N + i*P*P] - arr[j + (1)*N + i*P*P]);
            // arr[j + (N-3)*N + i*P*P] = arr[j + (N-4)*N + i*P*P] - (arr[j + (N-5)*N + i*P*P] - arr[j + (N-4)*N + i*P*P]);
            // arr[j + (N-2)*N + i*P*P] = arr[j + (N-3)*N + i*P*P] - (arr[j + (N-4)*N + i*P*P] - arr[j + (N-3)*N + i*P*P]);
            // arr[j + (N-1)*N + i*P*P] = arr[j + (N-2)*N + i*P*P] - (arr[j + (N-3)*N + i*P*P] - arr[j + (N-2)*N + i*P*P]);

            // arr[j + i*N + (2)*P*P] = arr[j + i*N + (3)*P*P] - (arr[j + i*N + (4)*P*P] - arr[j + i*N + (3)*P*P]);
            // arr[j + i*N + (1)*P*P] = arr[j + i*N + (2)*P*P] - (arr[j + i*N + (3)*P*P] - arr[j + i*N + (2)*P*P]);
            // arr[j + i*N + (0)*P*P] = arr[j + i*N + (1)*P*P] - (arr[j + i*N + (2)*P*P] - arr[j + i*N + (1)*P*P]);
            // arr[j + i*N + (N-3)*P*P] = arr[j + i*N + (N-4)*P*P] - (arr[j + i*N + (N-5)*P*P] - arr[j + i*N + (N-4)*P*P]);
            // arr[j + i*N + (N-2)*P*P] = arr[j + i*N + (N-3)*P*P] - (arr[j + i*N + (N-4)*P*P] - arr[j + i*N + (N-3)*P*P]);
            // arr[j + i*N + (N-1)*P*P] = arr[j + i*N + (N-2)*P*P] - (arr[j + i*N + (N-3)*P*P] - arr[j + i*N + (N-2)*P*P]);
        }
    }
}

void euler_upwind_reinit(vector<double> &arr, int M, int N, int P, double dx, double dy, double dz, double dt, const vector<double> &phi0){
    vector<double> phi;
    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                if (i==0 || i==(M-1) || i==1 || i==M-2 || i==2 || i==M-3 || j==0 || j==(N-1) || j==1 || j==N-2 || j==2 || j==N-3 || k==0 || k==(P-1) || k==1 || k==P-2 || k==2 || k==P-3){
                    phi.push_back(arr[i+j*N+k*P*P]);
                    continue;
                }

                double a = (arr[i+j*N+k*P*P]-arr[(i-1)+j*N+k*P*P])/dx;
                double b = (arr[(i+1)+j*N+k*P*P]-arr[i+j*N+k*P*P])/dx;
                double c = (arr[i+j*N+k*P*P]-arr[i+(j-1)*N+k*P*P])/dy;
                double d = (arr[i+(j+1)*N+k*P*P]-arr[i+j*N+k*P*P])/dy;
                double e = (arr[i+j*N+k*P*P]-arr[i+j*N+(k-1)*P*P])/dz;
                double f = (arr[i+j*N+(k+1)*P*P]-arr[i+j*N+k*P*P])/dz;

                double G;
                if (phi0[i+j*N+k*P*P] > 0){
                    G = sqrt(max(max(a,0.0)*max(a,0.0), min(b,0.0)*min(b,0.0))
                        + max(max(c,0.0)*max(c,0.0), min(d,0.0)*min(d,0.0))
                        + max(max(e,0.0)*max(e,0.0), min(f,0.0)*min(f,0.0))) - 1;
                } else if (phi0[i+j*N+k*P*P] < 0){
                    G = sqrt(max(min(a,0.0)*min(a,0.0), max(b,0.0)*max(b,0.0))
                    + max(min(c,0.0)*min(c,0.0), max(d,0.0)*max(d,0.0))
                    + max(min(e,0.0)*min(e,0.0), max(f,0.0)*max(f,0.0))) - 1;
                }
                phi.push_back(arr[i + j*N + k*P*P] - dt*sign(phi0[i+j*N+k*P*P])*G);
            }
        }
    }
    arr = phi;
}

void second_Order_Reinit(vector<double> &arr, vector<double> X, vector<double> Y, vector<double> Z, int M, int N, int P, double dx, double dy, double dz, double dt, const vector<double> &phi0){
    vector<double> phi;
    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                if (i==0 || i==1 || i==(M-1) || i==(M-2) || j==0 || j==1 || j==(N-1) || j==(N-2) || k==0 || k==1 || k==(P-1) || k==(P-2)){
                    phi.push_back(arr[i+j*N+k*P*P]);
                    continue;
                }
                
                // Nytt
                if (arr[i + j*N + k*P*P]*arr[(i-1) + j*N + k*P*P] < 0 || arr[i + j*N + k*P*P]*arr[(i+1) + j*N + k*P*P] < 0 || 
                    arr[i + j*N + k*P*P]*arr[i + (j-1)*N + k*P*P] < 0 || arr[i + j*N + k*P*P]*arr[i + (j+1)*N + k*P*P] < 0 ||
                    arr[i + j*N + k*P*P]*arr[i + j*N + (k-1)*P*P] < 0 || arr[i + j*N + k*P*P]*arr[i + j*N + (k+1)*P*P] < 0){
                    
                    double D = 2*dx*phi0[i+j*N+k*P*P]/cbrt(pow((phi0[(i+1)+j*N+k*P*P]-phi0[(i-1)+j*N+k*P*P]),3)
                    + pow((phi0[i+(j+1)*N+k*P*P]-phi0[i+(j-1)*N+k*P*P]),3) + pow((phi0[i+j*N+(k+1)*P*P]-phi0[i+j*N+(k-1)*P*P]),3));
                    phi.push_back(arr[i+j*N+k*P*P] - ((dt/dx)*sign(phi0[i+j*N+k*P*P])*abs(arr[i+j*N+k*P*P])-D));
                
                } else { // slutt nytt
                    double a_x, a_y, a_z, b_x, b_y, b_z;

                    {
                    // x-direction
                    vector<double> xp;
                    vector<double> yp;

                    if ((X[i]*X[i+1] < 0) || (X[i-1]*X[i] < 0)){
                        double f1_0_1 = (arr[(i-1)+j*N+k*P*P]-arr[(i)+j*N+k*P*P])/(X[i-1] - X[i]);
                        double f1_1_2 = (arr[(i)+j*N+k*P*P]-arr[(i+1)+j*N+k*P*P])/(X[i] - X[i+1]);
                        double f1_2_3 = (arr[(i+1)+j*N+k*P*P]-arr[(i+2)+j*N+k*P*P])/(X[i+1] - X[i+2]);
                        double f2_0_1_2 = (f1_1_2 - f1_0_1)/(X[i+1] - X[i-1]);
                        double f2_1_2_3 = (f1_2_3 - f1_1_2)/(X[i+2] - X[i]);
                        double f3_0_1_2_3 = (f2_1_2_3 - f2_0_1_2)/(X[i+2] - X[i-1]);
                        double x0 = arr[(i)+j*N+k*P*P] + f1_0_1*(-X[i-1]) + f2_0_1_2*(X[i-1]*X[i]) + f3_0_1_2_3*(-X[i-1]*X[i]*X[i+1]); //interface location

                        if (x0 < X[i]){
                            xp.push_back(X[i-1]);
                            xp.push_back(x0);
                            xp.push_back(X[i]);
                            xp.push_back(X[i+1]);
                            xp.push_back(X[i+2]);

                            yp.push_back(arr[(i-1)+j*N+k*P*P]);
                            yp.push_back(0);
                            yp.push_back(arr[i+j*N+k*P*P]);
                            yp.push_back(arr[(i+1)+j*N+k*P*P]);
                            yp.push_back(arr[(i+2)+j*N+k*P*P]);
                        } else {
                            xp.push_back(X[i-2]);
                            xp.push_back(X[i-1]);
                            xp.push_back(X[i]);
                            xp.push_back(x0);
                            xp.push_back(X[i+1]);

                            yp.push_back(arr[(i-2)+j*N+k*P*P]);
                            yp.push_back(arr[(i-1)+j*N+k*P*P]);
                            yp.push_back(arr[i+j*N+k*P*P]);
                            yp.push_back(0);
                            yp.push_back(arr[(i+1)+j*N+k*P*P]);   
                        }
                    } else {
                        xp.push_back(X[i-2]);
                        xp.push_back(X[i-1]);
                        xp.push_back(X[i]);
                        xp.push_back(X[i+1]);
                        xp.push_back(X[i+2]);

                        yp.push_back(arr[(i-2)+j*N+k*P*P]);
                        yp.push_back(arr[(i-1)+j*N+k*P*P]);
                        yp.push_back(arr[i+j*N+k*P*P]);
                        yp.push_back(arr[(i+1)+j*N+k*P*P]);
                        yp.push_back(arr[(i+2)+j*N+k*P*P]);
                    }
                    vector<double> theta1;
                    vector<double> theta2;
                    for (int a = 0; a < 5; ++a){
                        theta1.push_back((yp[a+1] - yp[a])/(xp[a+1] - xp[a]));
                    }
                    for (int a = 0; a < 4; ++a){
                        theta2.push_back((theta1[a+1] - theta1[a])/(xp[a+2] - xp[a]));
                    }
                    double c_m = minMod(theta2[0], theta2[1]);
                    double c_p = minMod(theta2[1], theta2[2]);
                    double a_x = theta1[1] + c_m*(xp[2] - xp[1]); // D_x^-phi_i 
                    double b_x = theta1[2] + c_p*(xp[2] - xp[3]); // D_x^+phi_i
                    }

                    {                

                    // y-direction
                    vector<double> xp;
                    vector<double> yp;
                    if ((Y[j]*Y[j+1] < 0) || (Y[j-1]*Y[j] < 0)){
                        double f1_0_1 = (arr[i+(j-1)*N+k*P*P]-arr[i+j*N+k*P*P])/(Y[j-1] - Y[j]);
                        double f1_1_2 = (arr[i+j*N+k*P*P]-arr[i+(j+1)*N+k*P*P])/(Y[j] - Y[j+1]);
                        double f1_2_3 = (arr[i+(j+1)*N+k*P*P]-arr[i+(j+2)*N+k*P*P])/(Y[j+1] - Y[j+2]);
                        double f2_0_1_2 = (f1_1_2 - f1_0_1)/(Y[j+1] - Y[j-1]);
                        double f2_1_2_3 = (f1_2_3 - f1_1_2)/(Y[j+2] - Y[j]);
                        double f3_0_1_2_3 = (f2_1_2_3 - f2_0_1_2)/(Y[j+2] - Y[j-1]);
                        double x0 = arr[i+j*N+k*P*P] + f1_0_1*(-Y[j-1]) + f2_0_1_2*(Y[j-1]*Y[j]) + f3_0_1_2_3*(-Y[j-1]*Y[j]*Y[j+1]); //interface location

                        if (x0 < Y[j]){
                            xp.push_back(Y[j-1]);
                            xp.push_back(x0);
                            xp.push_back(Y[j]);
                            xp.push_back(Y[j+1]);
                            xp.push_back(Y[j+2]);

                            yp.push_back(arr[i+(j-1)*N+k*P*P]);
                            yp.push_back(0);
                            yp.push_back(arr[i+j*N+k*P*P]);
                            yp.push_back(arr[i+(j+1)*N+k*P*P]);
                            yp.push_back(arr[i+(j+2)*N+k*P*P]);
                        } else {
                            xp.push_back(Y[j-2]);
                            xp.push_back(Y[j-1]);
                            xp.push_back(Y[j]);
                            xp.push_back(x0);
                            xp.push_back(Y[j+1]);

                            yp.push_back(arr[i+(j-2)*N+k*P*P]);
                            yp.push_back(arr[i+(j-1)*N+k*P*P]);
                            yp.push_back(arr[i+j*N+k*P*P]);
                            yp.push_back(0);
                            yp.push_back(arr[i+(j+1)*N+k*P*P]);   
                        }
                    } else {
                        xp.push_back(Y[j-2]);
                        xp.push_back(Y[j-1]);
                        xp.push_back(Y[j]);
                        xp.push_back(Y[j+1]);
                        xp.push_back(Y[j+2]);

                        yp.push_back(arr[i+(j-2)*N+k*P*P]);
                        yp.push_back(arr[i+(j-1)*N+k*P*P]);
                        yp.push_back(arr[i+j*N+k*P*P]);
                        yp.push_back(arr[i+(j+1)*N+k*P*P]);
                        yp.push_back(arr[i+(j+2)*N+k*P*P]);
                    }
                    vector<double> theta1;
                    vector<double> theta2;
                    for (int a = 0; a < 5; ++a){
                        theta1.push_back((yp[a+1] - yp[a])/(xp[a+1] - xp[a]));
                    }
                    for (int a = 0; a < 4; ++a){
                        theta2.push_back((theta1[a+1] - theta1[a])/(xp[a+2] - xp[a]));
                    }
                    double c_m = minMod(theta2[0], theta2[1]);
                    double c_p = minMod(theta2[1], theta2[2]);
                    double a_y = theta1[1] + c_m*(xp[2] - xp[1]); // D_y^-phi_j 
                    double b_y = theta1[2] + c_p*(xp[2] - xp[3]); // D_y^+phi_j

                    }

                    {
                    // z-direction

                    vector<double> xp;
                    vector<double> yp;

                    if ((Z[k]*Z[k+1] < 0) || (Z[k-1]*Z[k] < 0)){
                        double f1_0_1 = (arr[i+j*N+(k-1)*P*P]-arr[i+j*N+k*P*P])/(Z[k-1] - Z[k]);
                        double f1_1_2 = (arr[i+j*N+k*P*P]-arr[i+j*N+(k+1)*P*P])/(Z[k] - Z[k+1]);
                        double f1_2_3 = (arr[i+j*N+(k+1)*P*P]-arr[i+j*N+(k+2)*P*P])/(Z[k+1] - Z[k+2]);
                        double f2_0_1_2 = (f1_1_2 - f1_0_1)/(Z[k+1] - Z[k-1]);
                        double f2_1_2_3 = (f1_2_3 - f1_1_2)/(Z[k+2] - Z[k]);
                        double f3_0_1_2_3 = (f2_1_2_3 - f2_0_1_2)/(Z[k+2] - Z[k-1]);
                        double x0 = arr[i+j*N+k*P*P] + f1_0_1*(-Z[k-1]) + f2_0_1_2*(Z[k-1]*Z[k]) + f3_0_1_2_3*(-Z[k-1]*Z[k]*Z[k+1]); //interface location

                        if (x0 < Z[k]){
                            xp.push_back(Z[k-1]);
                            xp.push_back(x0);
                            xp.push_back(Z[k]);
                            xp.push_back(Z[k+1]);
                            xp.push_back(Z[k+2]);

                            yp.push_back(arr[i+j*N+(k-1)*P*P]);
                            yp.push_back(0);
                            yp.push_back(arr[i+j*N+k*P*P]);
                            yp.push_back(arr[i+j*N+(k+1)*P*P]);
                            yp.push_back(arr[i+j*N+(k+2)*P*P]);
                        } else {
                            xp.push_back(Z[k-2]);
                            xp.push_back(Z[k-1]);
                            xp.push_back(Z[k]);
                            xp.push_back(x0);
                            xp.push_back(Z[k+1]);

                            yp.push_back(arr[i+j*N+(k-2)*P*P]);
                            yp.push_back(arr[i+j*N+(k-1)*P*P]);
                            yp.push_back(arr[i+j*N+k*P*P]);
                            yp.push_back(0);
                            yp.push_back(arr[i+j*N+(k+1)*P*P]);   
                        }
                    } else {
                        xp.push_back(Z[k-2]);
                        xp.push_back(Z[k-1]);
                        xp.push_back(Z[k]);
                        xp.push_back(Z[k+1]);
                        xp.push_back(Z[k+2]);

                        yp.push_back(arr[i+j*N+(k-2)*P*P]);
                        yp.push_back(arr[i+j*N+(k-1)*P*P]);
                        yp.push_back(arr[i+j*N+k*P*P]);
                        yp.push_back(arr[i+j*N+(k+1)*P*P]);
                        yp.push_back(arr[i+j*N+(k+2)*P*P]);
                    }
                    vector<double> theta1;
                    vector<double> theta2;
                    for (int a = 0; a < 5; ++a){
                        theta1.push_back((yp[a+1] - yp[a])/(xp[a+1] - xp[a]));
                    }
                    for (int a = 0; a < 4; ++a){
                        theta2.push_back((theta1[a+1] - theta1[a])/(xp[a+2] - xp[a]));
                    }
                    double c_m = minMod(theta2[0], theta2[1]);
                    double c_p = minMod(theta2[1], theta2[2]);
                    double a_z = theta1[1] + c_m*(xp[2] - xp[1]); // D_z^-phi_i 
                    double b_z = theta1[2] + c_p*(xp[2] - xp[3]); // D_z^+phi_i
                    
                    }

                    double G;
                    if (phi0[i+j*N+k*P*P] > 0){
                        G = sqrt(max(max(a_x,0.0)*max(a_x,0.0), min(b_x,0.0)*min(b_x,0.0))
                            + max(max(a_y,0.0)*max(a_y,0.0), min(b_y,0.0)*min(b_y,0.0))
                            + max(max(a_z,0.0)*max(a_z,0.0), min(b_z,0.0)*min(b_z,0.0))) - 1;
                    } else if (phi0[i+j*N+k*P*P] < 0){
                        G = sqrt(max(min(a_x,0.0)*min(a_x,0.0), max(b_x,0.0)*max(b_x,0.0))
                            + max(min(a_y,0.0)*min(a_y,0.0), max(b_y,0.0)*max(b_y,0.0))
                            + max(min(a_z,0.0)*min(a_z,0.0), max(b_z,0.0)*max(b_z,0.0))) - 1;
                    }
                    phi.push_back(arr[i + j*N + k*P*P] - dt*sign(phi0[i+j*N+k*P*P])*G);
                }
            }
        }
    }
    arr = phi;
}

int sign(double num){
    int res;
    if (num < 0){
        res = -1;
    } else if (num > 0){
        res = 1;
    } else if (num == 0){
        res = 0;
    } 
    return res;
}

double minAbs(double a, double b){
    if (abs(a) < abs(b)){
        return a;
    } else {
        return b;
    }
}

double minMod(double a, double b){
    double res;
    if (abs(a) <= abs(b) && a*b > 0){
        res = a;
    } else if (abs(a) > abs(b) && a*b > 0){
        res = b;
    } else if (a*b <= 0) {
        res = 0;
    }
    return res;
}
