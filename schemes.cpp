#include "schemes.h"

Derivative upwind(vector<double> &phi, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz){

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
                    phix.push_back((phi[i + j*N + k*P*P] - phi[(i - 1) + j*N + k*P*P])/dx);
                } else if (AX[i + j*N + k*P*P] < 0){
                    phix.push_back((phi[(i + 1) + j*N + k*P*P] - phi[i + j*N + k*P*P])/dx);
                }

                if (AY[i + j*N + k*P*P] >= 0){
                    phiy.push_back((phi[i + j*N + k*P*P] - phi[i + (j - 1)*N + k*P*P])/dy);
                } else if (AY[i + j*N + k*P*P] < 0){
                    phiy.push_back((phi[i + (j + 1)*N + k*P*P] - phi[i + j*N + k*P*P])/dy);
                }

                if (AZ[i + j*N + k*P*P] >= 0){
                    phiz.push_back((phi[i + j*N + k*P*P] - phi[i + j*N + (k - 1)*P*P])/dz);
                } else if (AZ[i + j*N + k*P*P] < 0){
                    phiz.push_back((phi[i + j*N + (k + 1)*P*P] - phi[i + j*N + k*P*P])/dz);
                }
            }
        }
    }
    cout << AX.size() << " " << AY.size() << " " << AZ.size() << endl;
    return Derivative{phix, phiy, phiz};
}

Derivative weno(vector<double> &phi, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int const N, int const P, double dx, double dy, double dz){

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
                    v1 = (phi[(i-2) + j*N + k*P*P] - phi[(i-3) + j*N + k*P*P])/dx;
                    v2 = (phi[(i-1) + j*N + k*P*P] - phi[(i-2) + j*N + k*P*P])/dx;
                    v3 = (phi[(i) + j*N + k*P*P] - phi[(i-1) + j*N + k*P*P])/dx;
                    v4 = (phi[(i+1) + j*N + k*P*P] - phi[(i) + j*N + k*P*P])/dx;
                    v5 = (phi[(i+2) + j*N + k*P*P] - phi[(i+1) + j*N + k*P*P])/dx;
                } else if (AX[i + j*N + k*P*P] < 0){
                    v1 = (phi[(i-1) + j*N + k*P*P] - phi[(i-2) + j*N + k*P*P])/dx;
                    v2 = (phi[(i) + j*N + k*P*P] - phi[(i-1) + j*N + k*P*P])/dx;
                    v3 = (phi[(i+1) + j*N + k*P*P] - phi[(i) + j*N + k*P*P])/dx;
                    v4 = (phi[(i+2) + j*N + k*P*P] - phi[(i+1) + j*N + k*P*P])/dx;
                    v5 = (phi[(i+3) + j*N + k*P*P] - phi[(i+2) + j*N + k*P*P])/dx;
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
                    v1 = (phi[i + (j-2)*N + k*P*P] - phi[i + (j-3)*N + k*P*P])/dy;
                    v2 = (phi[i + (j-1)*N + k*P*P] - phi[i + (j-2)*N + k*P*P])/dy;
                    v3 = (phi[i + j*N + k*P*P] - phi[i + (j-1)*N + k*P*P])/dy;
                    v4 = (phi[i + (j+1)*N + k*P*P] - phi[i + j*N + k*P*P])/dy;
                    v5 = (phi[i + (j+2)*N + k*P*P] - phi[i + (j+1)*N + k*P*P])/dy;
                } else if (AY[i + j*N + k*P*P] < 0){
                    v1 = (phi[i + (j-1)*N + k*P*P] - phi[i + (j-2)*N + k*P*P])/dy;
                    v2 = (phi[i + j*N + k*P*P] - phi[i + (j-1)*N + k*P*P])/dy;
                    v3 = (phi[i + (j+1)*N + k*P*P] - phi[i + j*N + k*P*P])/dy;
                    v4 = (phi[i + (j+2)*N + k*P*P] - phi[i + (j+1)*N + k*P*P])/dy;
                    v5 = (phi[i + (j+3)*N + k*P*P] - phi[i + (j+2)*N + k*P*P])/dy;
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
                    v1 = (phi[i + j*N + (k-2)*P*P] - phi[i + j*N + (k-3)*P*P])/dz;
                    v2 = (phi[i + j*N + (k-1)*P*P] - phi[i + j*N + (k-2)*P*P])/dz;
                    v3 = (phi[i + j*N + k*P*P] - phi[i + j*N + (k-1)*P*P])/dz;
                    v4 = (phi[i + j*N + (k+1)*P*P] - phi[i + j*N + k*P*P])/dz;
                    v5 = (phi[i + j*N + (k+2)*P*P] - phi[i + j*N + (k+1)*P*P])/dz;
                } else if (AZ[i + j*N + k*P*P] < 0){
                    v1 = (phi[i + j*N + (k-1)*P*P] - phi[i + j*N + (k-2)*P*P])/dz;
                    v2 = (phi[i + j*N + k*P*P] - phi[i + j*N + (k-1)*P*P])/dz;
                    v3 = (phi[i + j*N + (k+1)*P*P] - phi[i + j*N + k*P*P])/dz;
                    v4 = (phi[i + j*N + (k+2)*P*P] - phi[i + j*N + (k+1)*P*P])/dz;
                    v5 = (phi[i + j*N + (k+3)*P*P] - phi[i + j*N + (k+2)*P*P])/dz;
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

Derivative godunov(vector<double> &phi, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int const N, int const P, double dx, double dy, double dz){

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

                double phix_m = (phi[i + j*N + k*P*P] - phi[(i-1) + j*N + k*P*P])/dx;
                double phix_p = (phi[(i+1) + j*N + k*P*P] - phi[i + j*N + k*P*P])/dx;

                double phiy_m = (phi[i + j*N + k*P*P] - phi[i + (j-1)*N + k*P*P])/dy;
                double phiy_p = (phi[i + (j+1)*N + k*P*P] - phi[i + j*N + k*P*P])/dy;

                double phiz_m = (phi[i + j*N + k*P*P] - phi[i + j*N + (k-1)*P*P])/dz;
                double phiz_p = (phi[i + j*N + (k+1)*P*P] - phi[i + j*N + k*P*P])/dz;

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

void euler_upwind(vector<double> &phi, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt){
    auto [phix, phiy, phiz] = upwind(phi, AX, AY, AZ, M, N, P, dx, dy, dz);
    phi = phi - dt*(AX*phix + AY*phiy + AZ*phiz);
}

void TVDRK3_upwind(vector<double> &phi, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt){

    vector<double> n1;
    vector<double> n2;
    vector<double> n3_2;

    {
        auto [phix, phiy, phiz] = upwind(phi, AX, AY, AZ, M, N, P, dx, dy, dz);
        n1 =  phi - dt*(AX*phix + AY*phiy + AZ*phiz);
    }

    {
        auto [phix, phiy, phiz] = upwind(n1, AX, AY, AZ, M, N, P, dx, dy, dz);
        n2 =  n1 - dt*(AX*phix + AY*phiy + AZ*phiz);
    }

    vector<double> n1_2 = 3/4*phi + 1/4*n2;

    {
        auto [phix, phiy, phiz] = upwind(n1_2, AX, AY, AZ, M, N, P, dx, dy, dz);
        n3_2 =  n1_2 - dt*(AX*phix + AY*phiy + AZ*phiz);
    }

    phi = 1/3*phi + 2/3*n3_2;

}

void TVDRK3_weno(vector<double> &phi, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt){

    {
    vector<double> n1;
    vector<double> n2;
    vector<double> n3_2;

    {
        auto [phix, phiy, phiz] = weno(phi, AX, AY, AZ, M, N, P, dx, dy, dz);
        n1 =  phi - dt*(AX*phix + AY*phiy + AZ*phiz);
    }

    {
        auto [phix, phiy, phiz] = weno(n1, AX, AY, AZ, M, N, P, dx, dy, dz);
        n2 =  n1 - dt*(AX*phix + AY*phiy + AZ*phiz);
    }

    vector<double> n1_2 = 3.0/4*phi + 1.0/4*n2;

    {
        auto [phix, phiy, phiz] = weno(n1_2, AX, AY, AZ, M, N, P, dx, dy, dz);
        n3_2 =  n1_2 - dt*(AX*phix + AY*phiy + AZ*phiz);
    }

    phi = 1.0/3*phi + 2.0/3*n3_2;

    }

    for(int i = 0; i < M; ++i){
        for (int j = 0; j < N; ++j){
            phi[2 + i*N + j*P*P] = phi[3 + i*N + j*P*P] - (phi[4 + i*N + j*P*P] - phi[3 + i*N + j*P*P]);
            phi[1 + i*N + j*P*P] = phi[2 + i*N + j*P*P] - (phi[3 + i*N + j*P*P] - phi[2 + i*N + j*P*P]);
            phi[0 + i*N + j*P*P] = phi[1 + i*N + j*P*P] - (phi[2 + i*N + j*P*P] - phi[1 + i*N + j*P*P]);
            phi[(M-3) + i*N + j*P*P] = phi[(M-4) + i*N + j*P*P] - (phi[(M-5) + i*N + j*P*P] - phi[(M-4) + i*N + j*P*P]);
            phi[(M-2) + i*N + j*P*P] = phi[(M-3) + i*N + j*P*P] - (phi[(M-4) + i*N + j*P*P] - phi[(M-3) + i*N + j*P*P]);
            phi[(M-1) + i*N + j*P*P] = phi[(M-2) + i*N + j*P*P] - (phi[(M-3) + i*N + j*P*P] - phi[(M-2) + i*N + j*P*P]);

            phi[j + (2)*N + i*P*P] = phi[j + (3)*N + i*P*P] - (phi[j + (4)*N + i*P*P] - phi[j + (3)*N + i*P*P]);
            phi[j + (1)*N + i*P*P] = phi[j + (2)*N + i*P*P] - (phi[j + (3)*N + i*P*P] - phi[j + (2)*N + i*P*P]);
            phi[j + (0)*N + i*P*P] = phi[j + (1)*N + i*P*P] - (phi[j + (2)*N + i*P*P] - phi[j + (1)*N + i*P*P]);
            phi[j + (N-3)*N + i*P*P] = phi[j + (N-4)*N + i*P*P] - (phi[j + (N-5)*N + i*P*P] - phi[j + (N-4)*N + i*P*P]);
            phi[j + (N-2)*N + i*P*P] = phi[j + (N-3)*N + i*P*P] - (phi[j + (N-4)*N + i*P*P] - phi[j + (N-3)*N + i*P*P]);
            phi[j + (N-1)*N + i*P*P] = phi[j + (N-2)*N + i*P*P] - (phi[j + (N-3)*N + i*P*P] - phi[j + (N-2)*N + i*P*P]);

            phi[j + i*N + (2)*P*P] = phi[j + i*N + (3)*P*P] - (phi[j + i*N + (4)*P*P] - phi[j + i*N + (3)*P*P]);
            phi[j + i*N + (1)*P*P] = phi[j + i*N + (2)*P*P] - (phi[j + i*N + (3)*P*P] - phi[j + i*N + (2)*P*P]);
            phi[j + i*N + (0)*P*P] = phi[j + i*N + (1)*P*P] - (phi[j + i*N + (2)*P*P] - phi[j + i*N + (1)*P*P]);
            phi[j + i*N + (N-3)*P*P] = phi[j + i*N + (N-4)*P*P] - (phi[j + i*N + (N-5)*P*P] - phi[j + i*N + (N-4)*P*P]);
            phi[j + i*N + (N-2)*P*P] = phi[j + i*N + (N-3)*P*P] - (phi[j + i*N + (N-4)*P*P] - phi[j + i*N + (N-3)*P*P]);
            phi[j + i*N + (N-1)*P*P] = phi[j + i*N + (N-2)*P*P] - (phi[j + i*N + (N-3)*P*P] - phi[j + i*N + (N-2)*P*P]);

        }
    }
}

void euler_weno(vector<double> &phi, vector<double> AX, vector<double> AY, vector<double> AZ, int M, int N, int P, double dx, double dy, double dz, double dt){
    auto [phix, phiy, phiz] = weno(phi, AX, AY, AZ, M, N, P, dx, dy, dz);
    phi = phi - dt*(AX*phix + AY*phiy + AZ*phiz);

    for(int i = 0; i < M; ++i){
        for (int j = 0; j < N; ++j){
            phi[2 + i*N + j*P*P] = phi[3 + i*N + j*P*P] - (phi[4 + i*N + j*P*P] - phi[3 + i*N + j*P*P]);
            phi[1 + i*N + j*P*P] = phi[2 + i*N + j*P*P] - (phi[3 + i*N + j*P*P] - phi[2 + i*N + j*P*P]);
            phi[0 + i*N + j*P*P] = phi[1 + i*N + j*P*P] - (phi[2 + i*N + j*P*P] - phi[1 + i*N + j*P*P]);
            phi[(M-3) + i*N + j*P*P] = phi[(M-4) + i*N - j*P*P] + (phi[(M-5) + i*N + j*P*P] - phi[(M-4) + i*N + j*P*P]);
            phi[(M-2) + i*N + j*P*P] = phi[(M-3) + i*N - j*P*P] + (phi[(M-4) + i*N + j*P*P] - phi[(M-3) + i*N + j*P*P]);
            phi[(M-1) + i*N + j*P*P] = phi[(M-2) + i*N - j*P*P] + (phi[(M-3) + i*N + j*P*P] - phi[(M-2) + i*N + j*P*P]);

            phi[j + (2)*N + i*P*P] = phi[j + (3)*N + i*P*P] - (phi[j + (4)*N + i*P*P] - phi[j + (3)*N + i*P*P]);
            phi[j + (1)*N + i*P*P] = phi[j + (2)*N + i*P*P] - (phi[j + (3)*N + i*P*P] - phi[j + (2)*N + i*P*P]);
            phi[j + (0)*N + i*P*P] = phi[j + (1)*N + i*P*P] - (phi[j + (2)*N + i*P*P] - phi[j + (1)*N + i*P*P]);
            phi[j + (N-3)*N + i*P*P] = phi[j + (N-4)*N + i*P*P] - (phi[j + (N-5)*N + i*P*P] - phi[j + (N-4)*N + i*P*P]);
            phi[j + (N-2)*N + i*P*P] = phi[j + (N-3)*N + i*P*P] - (phi[j + (N-4)*N + i*P*P] - phi[j + (N-3)*N + i*P*P]);
            phi[j + (N-1)*N + i*P*P] = phi[j + (N-2)*N + i*P*P] - (phi[j + (N-3)*N + i*P*P] - phi[j + (N-2)*N + i*P*P]);

            phi[j + i*N + (2)*P*P] = phi[j + i*N + (3)*P*P] - (phi[j + i*N + (4)*P*P] - phi[j + i*N + (3)*P*P]);
            phi[j + i*N + (1)*P*P] = phi[j + i*N + (2)*P*P] - (phi[j + i*N + (3)*P*P] - phi[j + i*N + (2)*P*P]);
            phi[j + i*N + (0)*P*P] = phi[j + i*N + (1)*P*P] - (phi[j + i*N + (2)*P*P] - phi[j + i*N + (1)*P*P]);
            phi[j + i*N + (N-3)*P*P] = phi[j + i*N + (N-4)*P*P] - (phi[j + i*N + (N-5)*P*P] - phi[j + i*N + (N-4)*P*P]);
            phi[j + i*N + (N-2)*P*P] = phi[j + i*N + (N-3)*P*P] - (phi[j + i*N + (N-4)*P*P] - phi[j + i*N + (N-3)*P*P]);
            phi[j + i*N + (N-1)*P*P] = phi[j + i*N + (N-2)*P*P] - (phi[j + i*N + (N-3)*P*P] - phi[j + i*N + (N-2)*P*P]);
        }
    }
}

void TVDRK3_godunov_reinit(vector<double> &phi, int M, int N, int P, double dx, double dy, double dz, double dt, vector<double> phi0){

    vector<double> n1;
    vector<double> n2;
    vector<double> n3_2;

    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                if (i==0 || i==(M-1) || i==1 || i==M-2 || i==2 || i==M-3 || j==0 || j==(N-1) || j==1 || j==N-2 || j==2 || j==N-3 || k==0 || k==(P-1) || k==1 || k==P-2 || k==2 || k==P-3){
                    n1.push_back(phi[i+j*N+k*P*P]);
                    continue;
                }

                double a = (phi[i+j*N+k*P*P]-phi[(i-1)+j*N+k*P*P])/dx;
                double b = (phi[(i+1)+j*N+k*P*P]-phi[i+j*N+k*P*P])/dx;
                double c = (phi[i+j*N+k*P*P]-phi[i+(j-1)*N+k*P*P])/dy;
                double d = (phi[i+(j+1)*N+k*P*P]-phi[i+j*N+k*P*P])/dy;
                double e = (phi[i+j*N+k*P*P]-phi[i+j*N+(k-1)*P*P])/dz;
                double f = (phi[i+j*N+(k+1)*P*P]-phi[i+j*N+k*P*P])/dz;

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
                n1.push_back(phi[i + j*N + k*P*P] - dt*sign(phi0[i+j*N+k*P*P])*G);
            }
        }
    }

    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                if (i==0 || i==(M-1) || i==1 || i==M-2 || i==2 || i==M-3 || j==0 || j==(N-1) || j==1 || j==N-2 || j==2 || j==N-3 || k==0 || k==(P-1) || k==1 || k==P-2 || k==2 || k==P-3){
                    n2.push_back(n1[i+j*N+k*P*P]);
                    continue;
                }

                double a = (n1[i+j*N+k*P*P]-n1[(i-1)+j*N+k*P*P])/dx;
                double b = (n1[(i+1)+j*N+k*P*P]-n1[i+j*N+k*P*P])/dx;
                double c = (n1[i+j*N+k*P*P]-n1[i+(j-1)*N+k*P*P])/dy;
                double d = (n1[i+(j+1)*N+k*P*P]-n1[i+j*N+k*P*P])/dy;
                double e = (n1[i+j*N+k*P*P]-n1[i+j*N+(k-1)*P*P])/dz;
                double f = (n1[i+j*N+(k+1)*P*P]-n1[i+j*N+k*P*P])/dz;

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
                n2.push_back(n1[i + j*N + k*P*P] - dt*sign(phi0[i+j*N+k*P*P])*G);
            }
        }
    }

    vector<double> n1_2 = (3.0/4*phi + 1.0/4*n2);

    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                if (i==0 || i==(M-1) || i==1 || i==M-2 || i==2 || i==M-3 || j==0 || j==(N-1) || j==1 || j==N-2 || j==2 || j==N-3 || k==0 || k==(P-1) || k==1 || k==P-2 || k==2 || k==P-3){
                    n3_2.push_back(n1_2[i+j*N+k*P*P]);
                    continue;
                }

                double a = (n1_2[i+j*N+k*P*P]-n1_2[(i-1)+j*N+k*P*P])/dx;
                double b = (n1_2[(i+1)+j*N+k*P*P]-n1_2[i+j*N+k*P*P])/dx;
                double c = (n1_2[i+j*N+k*P*P]-n1_2[i+(j-1)*N+k*P*P])/dy;
                double d = (n1_2[i+(j+1)*N+k*P*P]-n1_2[i+j*N+k*P*P])/dy;
                double e = (n1_2[i+j*N+k*P*P]-n1_2[i+j*N+(k-1)*P*P])/dz;
                double f = (n1_2[i+j*N+(k+1)*P*P]-n1_2[i+j*N+k*P*P])/dz;

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
                n3_2.push_back(n1_2[i + j*N + k*P*P] - dt*sign(phi0[i+j*N+k*P*P])*G);
            }
        }
    }

    phi = 1.0/3*phi + 2.0/3*n3_2;

}

void euler_upwind_reinit(vector<double> &phi, int M, int N, int P, double dx, double dy, double dz, double dt, const vector<double> &phi0){
    vector<double> phiNew;
    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                if (i==0 || i==(M-1) || i==1 || i==M-2 || i==2 || i==M-3 || j==0 || j==(N-1) || j==1 || j==N-2 || j==2 || j==N-3 || k==0 || k==(P-1) || k==1 || k==P-2 || k==2 || k==P-3){
                    phiNew.push_back(phi[i+j*N+k*P*P]);
                    continue;
                }

                double a = (phi[i+j*N+k*P*P]-phi[(i-1)+j*N+k*P*P])/dx;
                double b = (phi[(i+1)+j*N+k*P*P]-phi[i+j*N+k*P*P])/dx;
                double c = (phi[i+j*N+k*P*P]-phi[i+(j-1)*N+k*P*P])/dy;
                double d = (phi[i+(j+1)*N+k*P*P]-phi[i+j*N+k*P*P])/dy;
                double e = (phi[i+j*N+k*P*P]-phi[i+j*N+(k-1)*P*P])/dz;
                double f = (phi[i+j*N+(k+1)*P*P]-phi[i+j*N+k*P*P])/dz;

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
                phiNew.push_back(phi[i + j*N + k*P*P] - dt*sign(phi0[i+j*N+k*P*P])*G);
            }
        }
    }
    phi = phiNew;
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
