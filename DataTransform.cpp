#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <windows.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <array>

using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::Vector4d;

using namespace std;
typedef std::array<double, 8> Array;

const double PI = 3.14159265358979323846;
const double twoPI = 2 * PI;
const double deg2rad = PI / 180;
const double years2sec = 365.25 * 24 * 3600;
const double w0 = 2 * PI / 5400;
const double w02 = w0 * w0;

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

int main()
{
	Eigen::Matrix3d itm_B;
	itm_B(0, 0) = 250; itm_B(1, 1) = 250; itm_B(2, 2) = 100;  // inertia tensor
	itm_B(0, 1) = itm_B(0, 2) = itm_B(1, 0) = itm_B(1, 2) = itm_B(2, 0) = itm_B(2, 1) = 0;

	SetConsoleCP(1251);//меняем кодировку консоли принудительно
	SetConsoleOutputCP(1251);//меняем кодировку консоли принудительно на вывод				 
	ifstream file("asteroid.dat");//Создаем файловый поток и связываем его с файлом
	if (file.is_open())//Если открытие файла прошло успешно
	{
		//cout << "Файл открыт." << endl;
		string line;//Строчка текста
		while (getline(file, line))
		{
			//cout << line << endl;//Можно посмотреть, что в строчке считалось
			double t, h, x0, x1, x2, x3, x4, x5, x6, x7;
			istringstream iss(line);
			iss >> t >> h >> x0 >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7;
			Array x = { x0, x1, x2, x3, x4, x5, x6, x7 };

			//Calculating the YacobiIntegral
			Vector4d q = { x[0], x[1], x[2], x[3] };
			Vector3d omega = { x[4], x[5], x[6] };
			double nu = x[7];

			//get OX and OZ unit vectors of the orbital reference frame wrt PG frame
			Vector3d eRorb_B, eRorb_PG; //  OX unit vectors of the orbital reference frame wrt PG frame
			eRorb_PG = { cos(x[7]), sin(x[7]), 0 };
			eRorb_B = rotateByQ(q, eRorb_PG);

			Vector3d eZ_PG, eZ_B;
			eZ_PG = { 0, 0, 1 };
			eZ_B = rotateByQ(q, eZ_PG);

			Vector3d wt = omega - w0 * eZ_B;

			double Tr, Vg, T0;
			Tr = 0.5 * wt.transpose() * (itm_B * wt);
			Vg = 1.5 * w02 * eRorb_B.transpose() * (itm_B * eRorb_B);
			T0 = 0.5 * w02 * eZ_B.transpose() * (itm_B * eRorb_B);

			//cout << Tr + Vg - T0 << endl; // check conservation of YacobiIntegral for circular orbit

			Vector3d K_B; //angular momentum 
			K_B = itm_B * omega; 
			cout << K_B(2) << endl;// conservation for dynamically symmetric body 
		}
	}
	return 0;
}