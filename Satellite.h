#pragma once
#include <Eigen/Geometry>

//satellite parameters
struct satParam
{
    satParam()
    {
        
        itm_B(0, 0) = 50000; itm_B(1, 1) = 65000; itm_B(2, 2) = 70000;  // inertia tensor
        itm_B(0, 1) = itm_B(0, 2) = itm_B(1, 0) = itm_B(1, 2) = itm_B(2, 0) = itm_B(2, 1) = 0;
        
        
        /*
        itm_B(0, 0) = 50000; itm_B(1, 1) = 50000; itm_B(2, 2) = 70000;  // inertia tensor
        itm_B(0, 1) = itm_B(0, 2) = itm_B(1, 0) = itm_B(1, 2) = itm_B(2, 0) = itm_B(2, 1) = 0;
        */

        mag_B(0, 0) = 50000; mag_B(1, 1) = 65000; mag_B(2, 2) = 70000;  // magnetic tensor
        mag_B(0, 1) = mag_B(0, 2) = mag_B(1, 0) = mag_B(1, 2) = mag_B(2, 0) = itm_B(2, 1) = 0;

    }
	double mass;	
	int maxInertiaIndex;
	int minInertiaIndex;
	Eigen::Matrix3d itm_B;  //inertia tensor matrix
    Eigen::Matrix3d mag_B;  //magnetic tensor matrix
};

