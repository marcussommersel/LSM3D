#include "particleLSM.h"

// initializes particles in a cell where (x0, y0, z0) are the coordinates of the cell closest ot origo.
// dx, dy, dz are the grid spacing. X, Y, Z are all grid nodes in the computational domain.
// phi is the signed distance field. normal is the normal-vector. M, N, P is the number of grid nodes in each direction.
// numParticles are the number of particles of each type to be initialized in the cell.
vector<Particle> initializeParticles(double x0, double y0, double z0, double dx, double dy, double dz,
    vector<double> &X, vector<double> &Y, vector<double> &Z, vector<double> &phi, Derivative &normal,
    int M, int N, int P, int numParticles){

    vector<Particle> particles;
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, 100); // define the range
    double rmin = 0.1*min(dx, min(dy, dz)); // minimum particle radius
    double rmax = 0.5*max(dx, max(dy, dz)); // maximum particle radius
    double bmin = rmin;
    double bmax = 3.0*max(dx, max(dy, dz));
    double lambda = 1.0;
    double itmax = 15; // max iterations in the attraction step
    for (int p = 0; p < numParticles*2; ++p){
        int positive = p%2; // even negative, odd positive

        // random coordinate in a cell
        double x = x0 + dx*distr(gen)/100.0;
        double y = y0 + dy*distr(gen)/100.0;
        double z = z0 + dz*distr(gen)/100.0;
        double phip;

        // index position of the particle
        int i = (int)(x/dx);
        int j = (int)(y/dy);
        int k = (int)(z/dz);

        if (i < 0 || j < 0 || k < 0 || i >= M || j >= N || k >= P){continue;}

        phip = trilinearInterpolation(x, y, z, X[i], X[i+1], Y[j], Y[j+1], Z[k], Z[k+1],
            phi[i+j*N+k*P*P],
            phi[i+1+j*N+k*P*P],
            phi[i+1+(j+1)*N+k*P*P],
            phi[i+(j+1)*N+k*P*P],
            phi[i+j*N+(k+1)*P*P],
            phi[i+1+j*N+(k+1)*P*P],
            phi[i+1+(j+1)*N+(k+1)*P*P],
            phi[i+(j+1)*N+(k+1)*P*P]);

        double phiGoal = positive*(bmin + (bmax - bmin)*distr(gen)/100.0) - (1-positive)*(bmin + (bmax - bmin)*distr(gen)/100.0);

        // attraction step
        for (int it = 0; it < itmax; ++it){

            double Nxp = trilinearInterpolation(x, y, z, X[i], X[i+1], Y[j], Y[j+1], Z[k], Z[k+1],
                normal.x[i+j*N+k*P*P],
                normal.x[i+1+j*N+k*P*P],
                normal.x[i+1+(j+1)*N+k*P*P],
                normal.x[i+(j+1)*N+k*P*P],
                normal.x[i+j*N+(k+1)*P*P],
                normal.x[i+1+j*N+(k+1)*P*P],
                normal.x[i+1+(j+1)*N+(k+1)*P*P],
                normal.x[i+(j+1)*N+(k+1)*P*P]);
            x = x + lambda*(phiGoal - phip)*Nxp;

            double Nyp = trilinearInterpolation(x, y, z, X[i], X[i+1], Y[j], Y[j+1], Z[k], Z[k+1],
                normal.y[i+j*N+k*P*P],
                normal.y[i+1+j*N+k*P*P],
                normal.y[i+1+(j+1)*N+k*P*P],
                normal.y[i+(j+1)*N+k*P*P],
                normal.y[i+j*N+(k+1)*P*P],
                normal.y[i+1+j*N+(k+1)*P*P],
                normal.y[i+1+(j+1)*N+(k+1)*P*P],
                normal.y[i+(j+1)*N+(k+1)*P*P]);
            y = y + lambda*(phiGoal - phip)*Nyp;

            double Nzp = trilinearInterpolation(x, y, z, X[i], X[i+1], Y[j], Y[j+1], Z[k], Z[k+1],
                normal.z[i+j*N+k*P*P],
                normal.z[i+1+j*N+k*P*P],
                normal.z[i+1+(j+1)*N+k*P*P],
                normal.z[i+(j+1)*N+k*P*P],
                normal.z[i+j*N+(k+1)*P*P],
                normal.z[i+1+j*N+(k+1)*P*P],
                normal.z[i+1+(j+1)*N+(k+1)*P*P],
                normal.z[i+(j+1)*N+(k+1)*P*P]);
            z = z + lambda*(phiGoal - phip)*Nzp;

            i = (int)(x/dx);
            j = (int)(y/dy);
            k = (int)(z/dz);

            if (i < 0 || j < 0 || k < 0 || i >= M || j >= N || k >= P){break;}

            phip = trilinearInterpolation(x, y, z, X[i], X[i+1], Y[j], Y[j+1], Z[k], Z[k+1],
                phi[i+j*N+k*P*P],
                phi[i+1+j*N+k*P*P],
                phi[i+1+(j+1)*N+k*P*P],
                phi[i+(j+1)*N+k*P*P],
                phi[i+j*N+(k+1)*P*P],
                phi[i+1+j*N+(k+1)*P*P],
                phi[i+1+(j+1)*N+(k+1)*P*P],
                phi[i+(j+1)*N+(k+1)*P*P]);

            if ((positive && (phip >= bmin && phip <= bmax)) || (!positive && (phip <= -bmin && phip >= -bmax))){
                double r;
                if (sign(phip)*phip > rmax){
                    r = rmax;
                } else if (sign(phip)*phip < rmin){
                    r = rmin;
                } else {
                    r = sign(phip)*phip;
                }
                particles.push_back(Particle(x, y, z, r, positive));
                break;
            } else {
                lambda = lambda/2.0;
            }
        }
    }
    return particles;
}

// returns the corrected values of the signed distance field for each corner of the cells.
// the xi, yj, zk values are the coordinates of the corners of the cell, where x0 < x1, y0 < y1, z0 < z1.
// the phiijk are the values of the signed distance fiield at (i,j,k)
vector<double> correctInterface(Particle p, double x0, double x1, double y0, double y1, double z0, double z1,
    double phi000, double phi100, double phi110, double phi010, double phi001, double phi101, double phi111, double phi011, double phip){

    // distance from particle surface to cell corners
    double phip000 = sign(phip)*(p.r - sqrt(pow(x0 - p.x, 2) + pow(y0 - p.y, 2) + pow(z0 - p.z, 2)));
    double phip100 = sign(phip)*(p.r - sqrt(pow(x1 - p.x, 2) + pow(y0 - p.y, 2) + pow(z0 - p.z, 2)));
    double phip110 = sign(phip)*(p.r - sqrt(pow(x1 - p.x, 2) + pow(y1 - p.y, 2) + pow(z0 - p.z, 2)));
    double phip010 = sign(phip)*(p.r - sqrt(pow(x0 - p.x, 2) + pow(y1 - p.y, 2) + pow(z0 - p.z, 2)));
    double phip001 = sign(phip)*(p.r - sqrt(pow(x0 - p.x, 2) + pow(y0 - p.y, 2) + pow(z1 - p.z, 2)));
    double phip101 = sign(phip)*(p.r - sqrt(pow(x1 - p.x, 2) + pow(y0 - p.y, 2) + pow(z1 - p.z, 2)));
    double phip111 = sign(phip)*(p.r - sqrt(pow(x1 - p.x, 2) + pow(y1 - p.y, 2) + pow(z1 - p.z, 2)));
    double phip011 = sign(phip)*(p.r - sqrt(pow(x0 - p.x, 2) + pow(y1 - p.y, 2) + pow(z1 - p.z, 2)));

    vector<double> phi_p;
    vector<double> phi_m;
    vector<double> phi;

    phi_p.push_back(max(phi000, phip000));
    phi_p.push_back(max(phi100, phip100));
    phi_p.push_back(max(phi110, phip110));
    phi_p.push_back(max(phi010, phip010));
    phi_p.push_back(max(phi001, phip001));
    phi_p.push_back(max(phi101, phip101));
    phi_p.push_back(max(phi111, phip111));
    phi_p.push_back(max(phi011, phip011));

    phi_m.push_back(min(phi000, phip000));
    phi_m.push_back(min(phi100, phip100));
    phi_m.push_back(min(phi110, phip110));
    phi_m.push_back(min(phi010, phip010));
    phi_m.push_back(min(phi001, phip001));
    phi_m.push_back(min(phi101, phip101));
    phi_m.push_back(min(phi111, phip111));
    phi_m.push_back(min(phi011, phip011));

    for (int i = 0; i < phi_p.size(); ++i){
        if (abs(phi_p[i]) <= abs(phi_m[i])){
            phi.push_back(phi_p[i]);
        } else if (abs(phi_p[i]) > abs(phi_m[i])){
            phi.push_back(phi_m[i]);
        }
    }
    return phi;
}

// returns the interpolated value at (x, y, z)
// the xi, yj, zk values are the coordinates of the corners of the cell, where x0 < x1, y0 < y1, z0 < z1.
// the phiijk are the values of the signed distance fiield at (i,j,k)
double trilinearInterpolation(double x, double y, double z, double x0, double x1, double y0, double y1, double z0, double z1,
    double c000, double c100, double c110, double c010, double c001, double c101, double c111, double c011){

    double xd = (x - x0)/(x1 - x0);
    double yd = (y - y0)/(y1 - y0);
    double zd = (z - z0)/(z1 - z0);

    double c00 = c000*(1 - xd) + c100*xd;
    double c01 = c001*(1 - xd) + c101*xd;
    double c10 = c010*(1 - xd) + c110*xd;
    double c11 = c011*(1 - xd) + c111*xd;

    double c0 = c00*(1 - yd) + c10*yd;
    double c1 = c01*(1 - yd) + c11*yd;

    return c0*(1 - zd) + c1*zd;
}

// returns the normal vectors Nx, Ny, Nz for the signed distance field
Derivative normal(vector<double> &arr, double dx, double dy, double dz, double M, double N, double P){
    vector<double> Nx;
    vector<double> Ny;
    vector<double> Nz;
    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                if (i==0 || i==(M-1) || j==0 || j==(N-1) || k==0 || k==(P-1)){
                    Nx.push_back(0);
                    Ny.push_back(0);
                    Nz.push_back(0);
                    continue;
                }
                double phix = (arr[(i+1)+j*N+k*P*P] - arr[(i-1)+j*N+k*P*P])/(2*dx);
                if (phix == 0){
                    phix = (arr[(i+1)+j*N+k*P*P] - arr[i+j*N+k*P*P])/(dx);
                }
                double phiy = (arr[i+(j+1)*N+k*P*P] - arr[i+(j-1)*N+k*P*P])/(2*dy);
                if (phiy == 0){
                    phiy = (arr[i+(j+1)*N+k*P*P] - arr[i+j*N+k*P*P])/(dy);
                }
                double phiz = (arr[i+j*N+(k+1)*P*P] - arr[i+j*N+(k-1)*P*P])/(2*dz);
                if (phiz == 0){
                    phiz = (arr[i+j*N+(k+1)*P*P] - arr[i+j*N+k*P*P])/(dz);
                }
                Nx.push_back(phix/abs(phix));
                Ny.push_back(phiy/abs(phiy));
                Nz.push_back(phiz/abs(phiz));
            }
        }
    }
    return Derivative{Nx, Ny, Nz};
}

// saves the coordinates of all the particles to a .txt-files
void plotParticles(string filename, vector<Particle> particles){
    ofstream file;
    file.open(filename);
    if (!file.is_open()){cerr << "could not open file." << endl;}

    int count = 0;
    for (int i = 0; i < particles.size(); ++i){
        file << particles[i].x << "," << particles[i].y << "," << particles[i].z << "," << endl;
    }
    file.close();
}
