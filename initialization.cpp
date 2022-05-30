#include "initialization.h"

// addition of two points
Point Point::operator+(Point const &p){
    Point temp;
    temp.x = x + p.x;
    temp.y = y + p.y;
    temp.z = z + p.z;
    return temp;
}

// a point is assigned the same coordinates as another point
void Point::operator=(Point const &p){
    x = p.x;
    y = p.y;
    z = p.z;
}

// returns length between two points
double length(Point const &p0, Point const &p1){
    return sqrt((p0.x - p1.x)*(p0.x - p1.x) + (p0.y - p1.y)*(p0.y - p1.y) + (p0.z - p1.z)*(p0.z - p1.z));
}

// check if a point is within a sphere of center c and radius r
bool isInsideSphere(double r, Point c, Point p){
    if ((length(c, p) - r) < 0){
        return true;
    } else {
        return false;
    }
}

// returns the signed distance from a point to the surface of a sphere of center c and radius r
double signedDistanceSphere(double r, Point c, Point p){

    if (isInsideSphere(r, c, p)){
        return -(r - length(p, c));
    }
    return length(p, c) - r;
}


// generates a signed distance field for all points in [xmin, xmax] * [ymin, ymax] * [zmin, zmax] with reference to a sphere of center c and radius r
void signedDistanceField(vector<double> &arr, vector<double> x, vector<double> y, vector<double> z, double r, Point c, int M, int N, int P){ // Fix
    for (int k = 0; k < P; ++k){
        for (int j = 0; j < N; ++j){
            for (int i = 0; i < M; ++i){
                arr.push_back(signedDistanceSphere(r, c, Point(x[i], y[j], z[k])));
            }
        }
    }
}

// returns a vector of n indexes with equally spaced values from start to end
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

// saves a signed distance field to .txt-file
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
