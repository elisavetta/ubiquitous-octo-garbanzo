#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "Integrator.h"
#include "Defs.h"
#include "Environment.h"
#include "Orbit.h"
#include "Satellite.h"

using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::Vector4d;
using namespace std;

int sign(double a)
{
	return ((a > 0) - (a < 0));
}

Vector3d rotateByQ(const Vector4d& q, const Vector3d& r)
{
	Vector3d x;
	x(0) = (q(0) * q(0) + q(1) * q(1) - q(2) * q(2) - q(3) * q(3)) * r(0) + 2 * (q(1) * q(2) + q(0) * q(3)) * r(1) + 2 * (q(1) * q(3) - q(0) * q(2)) * r(2);
	x(1) = 2 * (q(1) * q(2) - q(0) * q(3)) * r(0) + (q(0) * q(0) + q(2) * q(2) - q(1) * q(1) - q(3) * q(3)) * r(1) + 2 * (q(2) * q(3) + q(0) * q(1)) * r(2);
	x(2) = 2 * (q(1) * q(3) + q(0) * q(2)) * r(0) + 2 * (q(2) * q(3) - q(0) * q(1)) * r(1) + (q(0) * q(0) + q(3) * q(3) - q(1) * q(1) - q(2) * q(2)) * r(2);
	return x;
}

void rhside(const satParam& sat, const envParam& env, orbitParam* orb, const double t, const Array& x, Array* dx)
{
	Vector4d q = { x[0], x[1], x[2], x[3] };
	Vector3d omega = { x[4], x[5], x[6] };
	Vector3d omega_dot;

	//recalculate OY and OZ unit vectors of the orbital reference frame wrt ECI frame
	Vector3d eZorb_B;
	orb->orbitalYZ(t);

	eZorb_B = rotateByQ(q, (*orb).eZorb_ECI);

	//gravity gradient torque
	Vector3d Tgg;
	//Tgg(0) = 0;
	//Tgg(1) = 0;
	//Tgg(2) = 0;
	Tgg = 3 * (*orb).w02 * eZorb_B.cross(sat.itm_B * eZorb_B);

	
	//torque due to eddy currents
	double rrr = (1340e3 + 6378.245e3) * (1340e3 + 6378.245e3) * (1340e3 + 6378.245e3);
	double coeff = 0.25 * env.vmu0 * env.vmue / (rrr *PI); // как тут правильно обращаться?
	Vector3d B;
	Vector3d ke_ECI; //dipole direction
	ke_ECI[0] = 0; // это правильно??? если прямой липоль да 
	ke_ECI[1] = 0;
	ke_ECI[2] = -1;
	Vector3d ke_B;
	ke_B = rotateByQ(q, ke_ECI);
	B = coeff * ((3 * eZorb_B * (ke_B.dot(eZorb_B))) - ke_B);
	Vector3d Tec;
	Tec = -B.cross(sat.mag_B*(omega.cross(B)));

	//gyroscpoic term
	Vector3d wJw;
	wJw = omega.cross(sat.itm_B * omega);

	omega_dot = sat.itm_B.inverse() * (-wJw + Tgg);
	//cout << t << "," << omega_dot[0] * omega_dot[0] + omega_dot[1] * omega_dot[1] + omega_dot[2] * omega_dot[2] << endl;

	(*dx)[0] = 0.5 * (-q[1] * omega[0] - q[2] * omega[1] - q[3] * omega[2]);
	(*dx)[1] = 0.5 * (q[0] * omega[0] + q[2] * omega[2] - q[3] * omega[1]);
	(*dx)[2] = 0.5 * (q[0] * omega[1] - q[1] * omega[2] + q[3] * omega[0]);
	(*dx)[3] = 0.5 * (q[0] * omega[2] + q[1] * omega[1] - q[2] * omega[0]);

	(*dx)[4] = omega_dot[0];
	(*dx)[5] = omega_dot[1];
	(*dx)[6] = omega_dot[2];
}

double YacobiIntegral(const orbitParam& orb, const satParam& sat, double t, const Array& x)
{
	Vector4d q = { x[0], x[1], x[2], x[3] };
	Vector3d omega = { x[4], x[5], x[6] };

	//get OY and OZ unit vectors of the orbital reference frame wrt ECI frame
	Vector3d eY_B, eZ_B;

	eY_B = rotateByQ(q, orb.eYorb_ECI);
	eZ_B = rotateByQ(q, orb.eZorb_ECI);

	Vector3d wt = omega - orb.w0 * eY_B;

	double Tr, Vg, T0;
	Tr = 0.5 * wt.transpose() * (sat.itm_B * wt);

	Vg = 1.5 * orb.w02 * eZ_B.transpose() * (sat.itm_B * eZ_B);

	T0 = 0.5 * orb.w02 * eY_B.transpose() * (sat.itm_B * eY_B);

	return Tr + Vg - T0;

	//return T0;
}

void output(const satParam& sat, const envParam& env, const orbitParam& orb, const double t, const Array& x, OutArray* out)
{
	Vector4d q = { x[0], x[1], x[2], x[3] };
	Vector3d omega = { x[4], x[5], x[6] };

	//get OY and OZ unit vectors of the orbital reference frame wrt ECI frame
	Vector3d eYorb_B, eZorb_B, eXsorb_B, eZsorb_ECI, eZsorb_B, K_B, eK_B;
	orbitParam loc_orb = orb;
	loc_orb.orbitalYZ(t);

	eZsorb_ECI(0) = cos(loc_orb.Omega0 + loc_orb.Omega_p * t);
	eZsorb_ECI(1) = sin(loc_orb.Omega0 + loc_orb.Omega_p * t);
	eZsorb_ECI(2) = 0;

	//eZsorb_ECI(0) = sin(loc_orb.Omega0 + loc_orb.Omega_p * t);
	//eZsorb_ECI(1) = 0;
	//eZsorb_ECI(2) = cos(loc_orb.Omega0 + loc_orb.Omega_p * t);

	eYorb_B = rotateByQ(q, loc_orb.eYorb_ECI);
	eZorb_B = rotateByQ(q, loc_orb.eZorb_ECI);
	eZsorb_B = rotateByQ(q, eZsorb_ECI);
	//eXorb_B = eYorb_B.cross(eZorb_B);
	eXsorb_B = eYorb_B.cross(eZsorb_B);

	K_B = sat.itm_B * omega;
	eK_B = K_B / sqrt(K_B(0) * K_B(0) + K_B(1) * K_B(1) + K_B(2) * K_B(2));

	double rho, sigma_s, nu;
	double sin_rho, sin_sigma_s, cos_sigma_s;

	rho = acos(eK_B.dot(eYorb_B));
	sin_rho = sin(rho);
	nu = acos(eK_B(2));

	if (abs(sin_rho) <= 1e-5)
	{
		sigma_s = 0;
	}
	else
	{

		
		cos_sigma_s = (eK_B.dot(eXsorb_B)) / sin_rho;
		sin_sigma_s = (eK_B.dot(eZsorb_B)) / sin_rho;

		if (abs(cos_sigma_s) > 1)
			sigma_s = 0.5 * PI * (1 - sign(cos_sigma_s));
		else
			sigma_s = acos(cos_sigma_s);

		if (sin_sigma_s < 0)
			sigma_s = 2 * PI - sigma_s;
		

		/*
		//SIDORENKO
		cos_sigma_s = (eK_B.dot(eZsorb_B)) / sin_rho;
		sin_sigma_s = (eK_B.dot(eXsorb_B)) / sin_rho;
		sigma_s = acos(cos_sigma_s);
		if (sin_sigma_s < 0)
			sigma_s = 2 * PI - sigma_s;
		*/

	}


	double abs_omega = sqrt(omega(0) * omega(0) + omega(1) * omega(1) + omega(2) * omega(2));

	//Vector3d el1_ECI, el1_B, el2_ECI, el2_B;
	//double ksi;
	//el1_ECI = { cos(rho) * sin(sigma_s), -sin(rho), cos(rho) * cos(sigma_s) };
	//el1_B = rotateByQ(q, el1_ECI);
	//el2_ECI = { cos(sigma_s), 0, -sin(sigma_s) };
	//el2_B = rotateByQ(q, el2_ECI);

	//if (abs(sin(nu)) <= 1e-5)
		//ksi = acos(el1_B(0));
	//else
		//ksi = atan(-el1_B(2) / eK_B(2));

	//(*out)[0] = el1_B(0)* el1_B(0) + el1_B(1)* el1_B(1) + el1_B(2)* el1_B(2);
	//(*out)[1] = eK_B(0)* eK_B(0)+ eK_B(1)* eK_B(1)+ eK_B(2)* eK_B(2);
	//(*out)[2] = el2_B(0) * el2_B(0) + el2_B(1) * el2_B(1) + el2_B(2) * el2_B(2);
	//(*out)[1] = abs_omega; // orb.w0;
	//(*out)[2] = twoPI / abs_omega;
	//(*out)[2] = atan(-el1_B(2)/ eK_B(2));
	//(*out)[3] = el1_B(2);
	//(*out)[1] = eK_B(2);

	(*out)[0] = t / 86400;
	(*out)[1] = rho / deg2rad; // Кассини +-
	(*out)[2] = sigma_s / deg2rad;
	(*out)[3] = abs_omega;
	//(*out)[3] = sqrt(K_B(0) * K_B(0) + K_B(1) * K_B(1) + K_B(2) * K_B(2)); //сохранение кинмомента при отсутствии грав момента +
	//(*out)[3] = K_B(3); // сохранение проекции кинмомента на ось с наиб моментом инерции (для симметричного спутника) +-
	//(*out)[3] = twoPI / abs_omega;
	//(*out)[4] = sigma_s / deg2rad;

	//(*out)[4] = delta / deg2rad;
	//(*out)[5] = delta_m / deg2rad;
	//(*out)[6] = theta / deg2rad;
	//(*out)[8] = sigma_b / deg2rad;
}

int main()
{
	satParam sat;
	envParam env;
	integParams solver;
	orbitInit orbInit;

	Integrator integrator;

	Array x;

	x[0] = 0.147016;
	x[1] = 0.691655;
	x[2] = 0.147016;
	x[3] = 0.691655;
	x[4] = 0;
	x[5] = 0;
	x[6] = -0.52359;

	OutArray x_out;

	orbitParam orb(&env, orbInit.h0_km, orbInit.i0_deg, orbInit.Omega0_deg, orbInit.u0_deg);
	
	std::cout.precision(6);
	std::cout.setf(std::ios::fixed);
	std::ofstream  resultFile("data.csv");
	
	/*
	std::cout.precision(6);
	std::cout.setf(std::ios::fixed);
	std::ofstream  resultFile("Yacobi.csv");
	*/

	double t = solver.t0;
	double h = solver.h;
	double t_end = solver.t_end * years2sec;
	int save_counter = 0;
	
	////cout << t << "," << YacobiIntegral(orb, sat, t, x) << endl;
	//cout << t << " " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << endl;
	/*
	cout << orbInit.h0_km << endl;
	cout << orbInit.i0_deg << endl;
	cout << orbInit.Omega0_deg << endl;
	cout << orbInit.u0_deg << endl;
	cout << orbInit.start_JD << endl;
	*/
	//resultFile << setprecision(8) << t << "," << YacobiIntegral(orb, sat, t, x) << endl;

	while (t < t_end)
	{
		double q_norm = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
		for (int i = 0; i < 4; i++) x[i] = x[i] / q_norm;
		if (t >= save_counter * solver.t_save)
		{
			output(sat, env, orb, t, x, &x_out);
			resultFile << setprecision(8) << x_out[0] << "," << x_out[1] << "," << x_out[2] << endl;
			//cout << t << " " << x[0] << " " << x[1] << " " << x[2] << endl;
			//resultFile << setprecision(8) << x_out[0] << "," << x_out[3] << endl;

			//cout << x_out[3] << endl;
			//cout << t << sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]) << endl;
			//resultFile << setprecision(8) << t << "," << YacobiIntegral(orb, sat, t, x) << endl;
			//cout << t << "," << YacobiIntegral(orb, sat, t, x) << endl; //сохранение интеграла Якоби для валидации
			//cout << t << "," << (init.h0_km * 1000 + env.re)/1000 << endl; //сохранение радиуса орбиты для валидации
			//cout << t << "," << orb.Omega << endl; //прецессия долготы восходящего узла для валидации
			save_counter += 1;
			//cout << t << endl;
		}
		integrator.rk78(&t, &x, &h, solver.tolerance, solver.hMin, solver.hMax, [&sat, &env, &orb](const double t, const Array& x, Array* dx) { rhside(sat, env, &orb, t, x, dx); });
		//cout << t << " " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << endl;
		//cout << t << "," << YacobiIntegral(orb, sat, t, x) << endl;
		//resultFile << setprecision(8) << t << "," << YacobiIntegral(orb, sat, t, x) << endl;

	}
	
	//double q_norm = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
	//for (int i = 0; i < 4; i++) x[i] = x[i] / q_norm;

	//output(sat, env, orb, t, x, &x_out);
	//resultFile << setprecision(8) << x_out[0] << "," << x_out[1] << "," << x_out[2] << "," << x_out[3] << "," << x_out[4] << endl;

	resultFile.close();
	//cout << YacobiIntegral(orb, sat, t, x) << endl;
	//cout << setprecision(10) << 0.52359*sin(PI / 3) * sin(PI / 6) << " " << 0.52359*sin(PI / 3) * cos(PI / 6) << " " << 0.52359*cos(PI / 3) << endl;
	//cout << setprecision(6) << 0.52359 * sin(PI / 3) << " " << 0.52359 * 0.5 << endl;
	return 0;
}

