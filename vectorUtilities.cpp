#include "vectorUtilities.h"

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
        res.push_back(vec[i]*scalar);
    }
    return res;
}

vector<double> operator/(vector<double> const &vec, double const &scalar){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(vec[i]/scalar);
    }
    return res;
}

vector<double> operator/(double const &scalar, vector<double> const &vec){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(scalar/vec[i]);
    }
    return res;
}

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

vector<double> operator+(vector<double> const &vec, double const &scalar){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(vec[i]+scalar);
    }
    return res;
}

vector<double> operator+(double const &scalar, vector<double> const &vec){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(scalar + vec[i]);
    }
    return res;
}

vector<double> operator-(vector<double> const &vec, double const &scalar){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(vec[i]-scalar);
    }
    return res;
}

vector<double> operator-(double const &scalar, vector<double> const &vec){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(scalar - vec[i]);
    }
    return res;
}

vector<double> vectorCBRT(vector<double> const &vec){
    vector<double> res;
    for (int i = 0; i < vec.size(); ++i){
        res.push_back(cbrt(vec[i]));
    }
    return res;
}

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

double vectorMax(vector<double> const &vec){
    double max = numeric_limits<double>::lowest();
    for (unsigned int i = 0; i < vec.size(); ++i){
        if (isfinite(vec[i]) && (vec[i] > max)){
            max = vec[i];
        }
    }
    return max;
}

vector<double> vectorSqrt(vector<double> const &vec){
    vector<double> res;
    for (unsigned int i = 0; i < vec.size(); ++i){
        res.push_back(sqrt(vec[i]));
    }
    return res;
}

// vector<double> vectorSign(vector<double> const &vec){
//     vector<double> s;
//     for (unsigned int i = 0; i < vec.size(); ++i){
//         if (vec)
//     }
// }