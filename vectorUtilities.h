#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
using namespace std;

// multiplies the elements of two vectors
vector<double> operator*(vector<double> const &vec1, vector<double> const &vec2);

// divides the elements of one vector with another vector
vector<double> operator/(vector<double> const &vec1, vector<double> const &vec2);

// adds the elements of two vectors
vector<double> operator+(vector<double> const &vec1, vector<double> const &vec2);

// multiplies the elements of a vector with a scalar
vector<double> operator*(double const &scalar, vector<double> const &vec);

// multiplies the elements of a vector with a scalar
vector<double> operator*(vector<double> const &vec, double const &scalar);

// divides the elements of a vector with a scalar
vector<double> operator/(vector<double> const &vec, double const &scalar);

// divides a scalar with the elements of a vector
vector<double> operator/(double const &scalar, vector<double> const &vec);

// subtracts the elements of one vector with another vector
vector<double> operator-(vector<double> const &vec1, vector<double> const &vec2);

// adds the elements of a vector with a scalar
vector<double> operator+(vector<double> const &vec, double const &scalar);

// adds the elements of a vector with a scalar
vector<double> operator+(double const &scalar, vector<double> const &vec);

// subtracts the elements of a vector with a scalar
vector<double> operator-(vector<double> const &vec, double const &scalar);

// subtracts a scalar with the elements of a vector
vector<double> operator-(double const &scalar, vector<double> const &vec);

// takes the absolute value of all elements of a vector
vector<double> vectorAbs(vector<double> const &vec);

// returns the maximum value of all elements in a vector
double vectorMax(vector<double> const &vec);

// takes the square root of all elements of a vector
vector<double> vectorSqrt(vector<double> const &vec);
