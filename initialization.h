#include <vector>
using namespace std;

class Point 
{
public:
    double x;
    double y;
};

class Points
{
public:
    vector<Point> points;
    double minDist(Point point, int n)
    {
        double min = INFINITY;
        for (const Point &p: this->points){
            // if ()
        }
    }
};