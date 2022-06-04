#pragma once
#include "Defs.h"
#include "Environment.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>

using Eigen::Vector3d;


struct orbitInit
{
    double h0_km = 1340.0;
    double i0_deg = 66.0;
    double Omega0_deg = 0.0;
    double u0_deg = 0.0;
};


//orbit paraneters
struct orbitParam
{
private:

public:
    orbitParam(envParam *env, double h0_km = 1340, double i0_deg = 66, double Omega0_deg = 0, double u0_deg = 0)
        : h0_km(h0_km), i0_deg(i0_deg), Omega0_deg(Omega0_deg), u0_deg(u0_deg)
    {
        double incl = i0_deg * deg2rad;
        Omega0 = Omega0_deg * deg2rad;
        r = h0_km * 1000 + env->re;
        r32 = sqrt(r * r * r);
        r72 = r * r * r * sqrt(r);
        cos_incl = cos(incl);
        sin_incl = sin(incl);
        u0 = u0_deg * deg2rad;

        // evolution of the orbit----------------------------
        vmu12 = sqrt(env->vmu);
        w0 = vmu12 / r32;
        w02 = w0 * w0;
        Omega_p = -1.5 * env->vJ2 * vmu12 * env->re * env->re * cos_incl / r72;
        wd = w0 * (1 - 1.5 * env->vJ2 * env->re * env->re * (1 - 4 * cos_incl * cos_incl) / r / r);

		orbitalYZ(0);
    
	}

	void orbitalYZ(const double t)
	{
		double u = u0 + t * wd;
		double cos_u = cos(u);
		double sin_u = sin(u);

		double Omega = Omega0 + Omega_p * t;
        	//double Omega = Omega0;
		double cos_Omega = cos(Omega);
		double sin_Omega = sin(Omega);

		eYorb_ECI(0) = sin_incl * sin_Omega;
		eYorb_ECI(1) = -sin_incl * cos_Omega;
		eYorb_ECI(2) = cos_incl;

		eZorb_ECI(0) = cos_u * cos_Omega - sin_u * cos_incl * sin_Omega;
		eZorb_ECI(1) = cos_u * sin_Omega + sin_u * cos_incl * cos_Omega;
		eZorb_ECI(2) = sin_u * sin_incl;
	}

    //auxiliary variables for orbit properties
    double r, r32, r72, cos_incl, sin_incl, u0, Omega0;

    // auxiliary variables for orbit evolution
    double	vmu12, w0, w02, Omega_p, wd; //1

    //initial values
    double h0_km;       // altitude, km
    double i0_deg;      // inclination, deg
    double Omega0_deg;  // ascending node longitude, deg
    double u0_deg;      // lattitude argument, deg
	double start_JD;    //initial datetime JD
	Vector3d eYorb_ECI;
	Vector3d eZorb_ECI;
};
