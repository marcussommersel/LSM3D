#include <vector>
#include <iostream>
#include "initialization2D.h"
#include "initialization2DTests.h"
#include "initialization.h"
#include "initializationTests.h"
#include "schemes.h"
#include "schemesTest.h"
#include "vectorUtilities.h"
using namespace std;

int main(){ 

    bool runTests = true;

    if (runTests)
    {
        // minDistTest();
        // interceptsTest();
        // isInsideTest();
        // signedDistanceTest();
        // simple2DAdvectionTest();
        // cubeTest();
        // isInsideSphereTest();
        // signedDistanceSphereTest();
        // signedDistanceFieldTest();
        euler_upwindTest();
    }


    return 0;
}