#include "initialization.h"

Point Point::operator+(Point const &p){
    Point temp;
    temp.x = x + p.x;
    temp.y = y + p.y;
    temp.z = z + p.z;
    return temp;
}

void Point::operator=(Point const &p){
    x = p.x;
    y = p.y;
    z = p.z;
}

double length(Point const &p0, Point const &p1){
    return sqrt((p0.x - p1.x)*(p0.x - p1.x) + (p0.y - p1.y)*(p0.y - p1.y) + (p0.z - p1.z)*(p0.z - p1.z));
}

bool isInsideSphere(double r, Point c, Point p){
    if ((length(c, p) - r) < 0){
        return true;
    } else {
        return false;
    }
}

double signedDistanceSphere(double r, Point c, Point p){

    if (isInsideSphere(r, c, p)){
        return -(r - length(p, c));
    }
    return length(p, c) - r;
}

void signedDistanceField(vector<double> &arr, vector<double> x, vector<double> y, vector<double> z, double r, Point c, int M, int N, int P){ // Fix
    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                arr.push_back(signedDistanceSphere(r, c, Point(x[i], y[j], z[k])));
            }
        }
    }
}

vector<double> linspace(double start, double end, int n){

    vector<double> vec;

    if (n == 0) { 
        return vec;
    }
    if (n == 1) {
      vec.push_back(start);
      return vec;
    }

    double dx = (end - start)/(n - 1);

    for(int i = 0; i < n - 1; ++i){
        vec.push_back(start + dx * i);
    }
    vec.push_back(end);

    return vec;
}

void saveScalarField(string filename, vector<double> const &arr, vector<double> x, vector<double> y, vector<double> z, int M, int N, int P){
    ofstream file;
    file.open(filename);
    if (!file.is_open()){cerr << "could not open file." << endl;}
        
    file << M << "," << N << "," << P << endl;
    int count = 0;
    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                file << x[i] << "," << y[j] << "," << z[k] << "," << arr[count] << "," << endl;
                ++count;
            }
        }
    }
    file.close();
}
