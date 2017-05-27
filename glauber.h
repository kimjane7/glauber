#ifndef __GLAUBER_H__
#define __GLAUBER_H__
#include <cstdio>
#include <boost/random.hpp>

const double pi = 4.0*atan(1.0);
using boost::uniform_01;
using boost::mt19937;

class CGlauber{
public:
	const static double a=0.53, rho0=0.14;
	double R;
	int A;
	CGlauber(int Aset);
	double get_rho(double x, double y, double z);
	double get_T(double x, double y);
	double random_01(mt19937 generator);
	void random_test(int nmax);

};

CGlauber::CGlauber(int A_set){

	A = A_set;
	double A_est, AP_est;
	double R_est = pow(3.0*A/(4.0*pi*rho0),1.0/3.0), R_diff;
	double r, dr = 1E-3*R_est;

	do{
		A_est = 0.0;
		for(r = 0.0; r < R_est; r += dr){
			A_est += 4.0*pi*rho0*r*r*dr/(1.0+exp((r-R_est)/a));
		}
		AP_est = 2.0*pi*rho0*R_est*R_est;
		R_diff = (A-A_est)/AP_est;
		R_est = R_est+R_diff;
	}while(R_diff > dr);

	R = R_est;
}

double CGlauber::get_rho(double x, double y, double z){
	double r = pow(x*x+y*y+z*z,0.5);
	double rho = rho0/(1.0+exp((r-R)/a));
	return rho;
}

double CGlauber::get_T(double x, double y){
	double T=0.0, z, z_bound = pow(R*R-x*x-y*y,0.5), dz = 2.0E-3*z_bound;
	for(z=-z_bound; z<z_bound; z+=dz){
		T+=get_rho(x,y,z)*dz;
	}
	return T;
}

double CGlauber::random_01(mt19937 generator){
	static uniform_01<mt19937> dist(generator);
	return dist();
}

void CGlauber::random_test(int nmax){
	mt19937 generator(time(0));
	printf("A = %i \nR = %lf fm\n", A, R);
	double r, theta, phi, x, y, z;
	for(int n = 0; n < nmax; n+=1){
		r = R*random_01(generator);
		theta = pi*random_01(generator);
		phi = pi*random_01(generator);
		x = r*sin(theta)*cos(phi);
		y = r*sin(theta)*sin(phi);
		z = r*cos(theta);
		printf("rho(%.3lf, %.3lf, %.3lf)=%lf \t T(%.3lf, %.3lf)=%lf\n",
			   x,y,z,get_rho(x,y,z),x,y,get_T(x,y));
	}
}

#endif