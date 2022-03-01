#pragma once

#include <vector>
#include <iostream>
using namespace std;

vector<double> operator*(vector<double> const &vec1, vector<double> const &vec2);

vector<double> operator+(vector<double> const &vec1, vector<double> const &vec2);

vector<double> operator*(double const &scalar, vector<double> const &vec);

vector<double> operator*(vector<double> const &vec, double const &scalar);

vector<double> operator-(vector<double> const &vec1, vector<double> const &vec2);