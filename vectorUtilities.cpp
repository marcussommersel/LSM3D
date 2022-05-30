#include "vectorUtilities.h"

// multiplies the elements of two vectors
vector<double> operator*(vector<double> const &vec1, vector<double> const &vec2){
    vector<double> res;
    if (vec1.size() != vec2.size()){
        cerr << "Vectors multiplication with different sized vectors." << endl;
    }
    for (int i = 0; i < vec1.size(); ++i){
        res.push_back(vec1[i]*vec2[i]);
    }
    return res;
}

// divides the elements of one vector with another vector
vector<double> operator/(vector<double> const &vec1, vector<double> const &vec2){
    vector<double> res;
    if (vec1.size() != vec2.size()){
        cerr << "Vectors division with different sized vectors." << endl;
    }
    for (int i = 0; i < vec1.size(); ++i){
        res.push_back(vec1[i]/vec2[i]);
    }
    return res;
}

// adds the elements of two vectors
vector<double> operator+(vector<double> const &vec1, vector<double> const &vec2){
    vector<double> res;
    if (vec1.size() != vec2.size()){
        cerr << "Vectors addition with different sized vectors." << endl;
    }
    for (int i = 0; i < vec1.size(); ++i){
        res.push_back(vec1[i] + vec2[i]);
    }
    return res;
}

// multiplies the elements of a vector with a scalar
vector<double> operator*(double const &scalar, vector<double> const &vec){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(scalar*vec[i]);
    }
    return res;
}

// multiplies the elements of a vector with a scalar
vector<double> operator*(vector<double> const &vec, double const &scalar){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(vec[i]*scalar);
    }
    return res;
}

// divides the elements of a vector with a scalar
vector<double> operator/(vector<double> const &vec, double const &scalar){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(vec[i]/scalar);
    }
    return res;
}

// divides a scalar with the elements of a vector
vector<double> operator/(double const &scalar, vector<double> const &vec){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(scalar/vec[i]);
    }
    return res;
}

// subtracts the elements of one vector with another vector
vector<double> operator-(vector<double> const &vec1, vector<double> const &vec2){
    vector<double> res;
    if (vec1.size() != vec2.size()){
        cerr << "Vectors subtraction with different sized vectors." << endl;
    }
    for (int i = 0; i < vec1.size(); ++i){
        res.push_back(vec1[i] - vec2[i]);
    }
    return res;
}

// adds the elements of a vector with a scalar
vector<double> operator+(vector<double> const &vec, double const &scalar){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(vec[i]+scalar);
    }
    return res;
}

// adds the elements of a vector with a scalar
vector<double> operator+(double const &scalar, vector<double> const &vec){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(scalar + vec[i]);
    }
    return res;
}

// subtracts the elements of a vector with a scalar
vector<double> operator-(vector<double> const &vec, double const &scalar){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(vec[i]-scalar);
    }
    return res;
}

// subtracts a scalar with the elements of a vector
vector<double> operator-(double const &scalar, vector<double> const &vec){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(scalar - vec[i]);
    }
    return res;
}

// takes the absolute value of all elements of a vector
vector<double> vectorAbs(vector<double> const &vec){
    vector<double> res;
    for (unsigned int i = 0; i < vec.size(); ++i){
        if (vec[i] < 0){
            res.push_back(-vec[i]);
        } else {
            res.push_back(vec[i]);
        }
    }
    return res;
}

// returns the maximum value of all elements in a vector
double vectorMax(vector<double> const &vec){
    double max = numeric_limits<double>::lowest();
    for (unsigned int i = 0; i < vec.size(); ++i){
        if (isfinite(vec[i]) && (vec[i] > max)){
            max = vec[i];
        }
    }
    return max;
}

// takes the square root of all elements of a vector
vector<double> vectorSqrt(vector<double> const &vec){
    vector<double> res;
    for (unsigned int i = 0; i < vec.size(); ++i){
        res.push_back(sqrt(vec[i]));
    }
    return res;
}
