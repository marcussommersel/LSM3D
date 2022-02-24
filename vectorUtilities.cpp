#include "vectorUtilities.h"

vector<double> operator*(vector<double> const &vec1, vector<double> const &vec2){
    vector<double> res;
    if (vec1.size() != vec2.size()){
        cerr << "Vectors multiplication with different sized vectors.";
    }
    for (int i = 0; i < vec1.size(); ++i){
        res.push_back(vec1[i]*vec2[i]);
    }
    return res;
}

vector<double> operator+(vector<double> const &vec1, vector<double> const &vec2){
    vector<double> res;
    if (vec1.size() != vec2.size()){
        cerr << "Vectors addition with different sized vectors.";
    }
    for (int i = 0; i < vec1.size(); ++i){
        res.push_back(vec1[i] + vec2[i]);
    }
    return res;
}

vector<double> operator*(double const &scalar, vector<double> const &vec){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(scalar*vec[i]);
    }
    return res;
}

vector<double> operator*(vector<double> const &vec, double const &scalar){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(scalar*vec[i]);
    }
    return res;
}

vector<double> operator-(vector<double> const &vec1, vector<double> const &vec2){
    vector<double> res;
    if (vec1.size() != vec2.size()){
        cerr << "Vectors subtraction with different sized vectors.";
    }
    for (int i = 0; i < vec1.size(); ++i){
        res.push_back(vec1[i] - vec2[i]);
    }
    return res;
}
