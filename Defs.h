#pragma once
#include <array>

const double PI = 3.14159265358979323846;
const double twoPI = 2 * PI;
const double deg2rad = PI / 180;
const double years2sec = 365.25 * 24 * 3600;
//const double c = 299792458; //speed of light
//const double G = 1361; //Solar constant (irradiance)
const double Psun = 4.54e-6; //Solar radiation pressure G/c
#define SGN(a) (((a)<0) ? -1 : 1)

const std::size_t Dimension = 7;
using Array = std::array<double, Dimension>;

typedef std::array<double, 7> OutArray;
