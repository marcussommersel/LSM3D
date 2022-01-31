#pragma once

#include "initialization.h"

void cubeTest(){
    double rho = 1;
    double theta;
    double phi;
    double x = rho*cos(theta)*sin(phi);
    double y = rho*sin(theta)*sin(phi);
    double z = rho*cos(phi);
}