#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <array>
#include "rk78.h"

const double PI = 3.14159265358979323846;
const double twoPI = 2 * PI;
const double deg2rad = PI / 180;
const double years2sec = 365.25 * 24 * 3600;

const double itm_A = 250;
const double itm_B = 250;
const double itm_C = 100;

const double mu = 3.986e14;     //the Earth gravity constant
const double eccentricity = 0;
const double w0 = 2 * PI / 5400;
const double w02 = w0 * w0;
const double a_major = pow(mu / w02, 1.0 / 3.0);//semi-major axis
const double ecc_param3 = (1 - eccentricity * eccentricity) * (1 - eccentricity * eccentricity) * (1 - eccentricity * eccentricity);

using namespace std;
using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::Vector4d;

Vector4d quatproduct(const Vector4d& q, const Vector4d& m) {
	Vector4d x; // function is used for finding initial quaternion
	x(0) = q(0) * m(0) - q(1) * m(1) - q(2) * m(2) - q(3) * m(3);
	x(1) = q(0) * m(1) + q(1) * m(0) + q(2) * m(3) - q(3) * m(2);
	x(2) = q(0) * m(2) + q(2) * m(0) + q(3) * m(1) - q(1) * m(3);
	x(3) = q(0) * m(3) + q(3) * m(0) + q(1) * m(2) - q(2) * m(1);
	return x;
}

template<class T>
T eqf0(T t, T y[8]) {
	double q_norm = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2] + y[3] * y[3]);
	for (int i = 0; i < 4; i++) y[i] = y[i] / q_norm;
    return 0.5 * (-y[1] * y[4] - y[2] * y[5] - y[3] * y[6]);
}

template<class T>
T eqf1(T t, T y[8]) {
	double q_norm = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2] + y[3] * y[3]);
	for (int i = 0; i < 4; i++) y[i] = y[i] / q_norm;
    return 0.5 * (y[0] * y[5] + y[2] * y[6] - y[3] * y[5]);
}

template<class T>
T eqf2(T t, T y[8]) {
	double q_norm = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2] + y[3] * y[3]);
	for (int i = 0; i < 4; i++) y[i] = y[i] / q_norm;
    return 0.5 * (y[0] * y[5] - y[1] * y[6] + y[3] * y[4]);
}

template<class T>
T eqf3(T t, T y[8]) {
	double q_norm = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2] + y[3] * y[3]);
	for (int i = 0; i < 4; i++) y[i] = y[i] / q_norm;
    return 0.5 * (y[0] * y[6] + y[1] * y[5] - y[2] * y[4]);
}

template<class T>
T eqf4(T t, T y[8]) {
    return (itm_B - itm_C)*y[5]*y[6]/itm_A;
}

template<class T>
T eqf5(T t, T y[8]) {
    return (itm_C - itm_A) * y[4] * y[6] / itm_B;
}

template<class T>
T eqf6(T t, T y[8]) {
    return (itm_A - itm_B) * y[5] * y[4] / itm_C;
}

template<class T>
T eqf7(T t, T y[8]) {
    return w0 * (1 + eccentricity * cos(y[7])) * (1 + eccentricity * cos(y[7])) / (sqrt(ecc_param3));
}

int main(int argc, char* argv[]) {
    const char* datafile = "asteroid.dat";
    RKF78<double, 8> RKF;
    RKF.f[0] = &eqf0<double>;
    RKF.f[1] = &eqf1<double>;
    RKF.f[2] = &eqf2<double>;
    RKF.f[3] = &eqf3<double>;
    RKF.f[4] = &eqf4<double>;
    RKF.f[5] = &eqf5<double>;
    RKF.f[6] = &eqf6<double>;
    RKF.f[7] = &eqf7<double>;

    //Initial conditions
	double sigma_0, ro_0, psi_0, fi_0, tetta_0;
	sigma_0 = 12 * deg2rad;
	ro_0 = 50 * deg2rad;
	psi_0 = 75 * deg2rad;
	tetta_0 = 30 * deg2rad;
	fi_0 = 20 * deg2rad; // if tetta_0 = PI*k, then fi_0 = 0 
	double mod_omega_0 = w0 * 20; //module of omega_0
	double nu_0 = 0; //initial value of true anomaly

	//omega_0 
	Eigen::Matrix3d itm_B_inv;
	itm_B_inv(0, 0) = 1 / itm_A; itm_B_inv(1, 1) = 1 / itm_B; itm_B_inv(2, 2) = 1 / itm_C;  // inertia tensor
	itm_B_inv(0, 1) = itm_B_inv(0, 2) = itm_B_inv(1, 0) = itm_B_inv(1, 2) = itm_B_inv(2, 0) = itm_B_inv(2, 1) = 0;
	Vector3d e_K = { sin(tetta_0) * sin(fi_0), sin(tetta_0) * cos(fi_0), cos(tetta_0) };
	Vector3d e_omega = itm_B_inv * e_K * (1 / mod_omega_0);
	double e_omega_norm = sqrt(e_omega[0] * e_omega[0] + e_omega[1] * e_omega[1] + e_omega[2] * e_omega[2]);
	for (int i = 0; i < 3; i++) e_omega[i] = e_omega[i] / e_omega_norm;

	Vector3d omega_0 = e_omega * mod_omega_0;

	//calculation of the initial quaternion value
	Vector4d q_1, q_2, q_3, q_4, q_5, q_6, q_7, q_8, q_0;
	q_1 = { cos(sigma_0 / 2), 0, 0, sin(sigma_0 / 2) };
	q_2 = { cos(ro_0 / 2), 0, sin(ro_0 / 2), 0 };
	q_3 = quatproduct(q_1, q_2);
	q_4 = { cos(psi_0 / 2), 0, 0, sin(psi_0 / 2) };
	q_5 = quatproduct(q_3, q_4);
	q_6 = { cos(tetta_0 / 2), sin(tetta_0 / 2), 0, 0 };
	q_7 = quatproduct(q_5, q_6);
	q_8 = { cos(fi_0 / 2), 0, 0, sin(fi_0 / 2) };
	q_0 = quatproduct(q_7, q_8); // initial quaternion value
	double q0_norm = sqrt(q_0[0] * q_0[0] + q_0[1] * q_0[1] + q_0[2] * q_0[2] + q_0[3] * q_0[3]);
	for (int i = 0; i < 4; i++) q_0[i] = q_0[i] / q0_norm;

    double rkf[8] = { q_0[0], q_0[1], q_0[2], q_0[3], omega_0[0], omega_0[1], omega_0[2], nu_0}; //final vector of initial conditions
    try {
        RKF.solve(1, 1e-3, rkf, 1e-9, 0.0, 1800, datafile);
    }
    catch (invalid_argument& e) {
        cerr << e.what() << endl;
        return -1;
    }
    return 0;
}
