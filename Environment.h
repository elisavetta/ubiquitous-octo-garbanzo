#pragma once
#include <cstdio>
#include <iostream>
#include <string>
#include "Defs.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>

using Eigen::Vector3d;
using namespace std;

const double Re = 6378.245e3;

//environment parameters
struct envParam {
//Earth rotation
	double omega_Earth = 7.292115e-5;     //the Earth rotation velocity
//Earth gravity parameters
	double vmu = 3.986e14;     //the Earth gravity constant
	double vJ2 = 1.082626e-3;  // J2 coefficient
	double re = Re;    // mean radius
//Earth magnetic parameters
	double vmue = 7.94e+22;
	double vmu0 = 12.5664e-7;

};