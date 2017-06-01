#ifndef GLAUBER_H
#define GLAUBER_H
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <boost/random.hpp>

const double pi = 4.0*atan(1.0);
// uniform rand dist over interval, takes RNG as arg
using boost::uniform_01;
// mersenne twister RNG
using namespace std;
typedef boost::mt19937 RNGType;

class CNucleus {
public:
	
	// declare member variables
	const static double a_ = 0.53, rho0_ = 0.14;
	double R_;
	int A_;

	// declare member functions
	CNucleus(int Aset);
	double get_R();
	double get_rho(double x, double y, double z);
	double get_T(double x, double y);
	double random_01(RNGType generator);
	void random_test(int nmax);
	void T_test(int nmax, string filename);

};

CNucleus::CNucleus(int A_set) {
	
	// set the nucleon number
	A_ = A_set;
	
	// estimate for A(R) and A'(R)
	double A_est, A_prime_est;
	
	// estimate radius from nucleon number
	double R_est = pow(3.0*A_/(4.0*pi*rho0_), 1.0/3.0);
	
	// newton's method increment
	double R_diff; 

	double r;
	double dr = 1E-3*R_est;

	// perform newton's method to find actual radius
	do {
		
		// calculate A from WS dist, density, and radius guess
		A_est = 0.0;
		for (r = 0.0; r < R_est; r += dr) {
			A_est += 4.0*pi*rho0_*r*r*dr/(1.0+exp((r-R_est)/a_));
		}
		
		// update radius guess
		A_prime_est = 2.0*pi*rho0_*R_est*R_est;
		R_diff = (A_-A_est)/A_prime_est;
		R_est = R_est+R_diff;
	
	} while(R_diff > dr);

	R_ = R_est;
}

double CNucleus::get_R() {
	return R_;
}

double CNucleus::get_rho(double x, double y, double z) {
	
	// calculate distance from center
	double r = pow(x*x+y*y+z*z,0.5);
	
	// calculate Woods-Saxon density
	double rho = rho0_/(1.0+exp((r-R_)/a_));

	return rho;
}

double CNucleus::get_T(double x, double y) {

	// initialize variables
	double T=0.0, z, z_bound = pow(R_*R_-x*x-y*y,0.5), dz = 2.0E-3*z_bound;
	
	// integrate the density in the z-direction to get thickness
	for(z=-z_bound; z<z_bound; z+=dz) {
		T+=get_rho(x,y,z)*dz;
	}
	
	return T;
}

double CNucleus::random_01(RNGType generator) {
	static uniform_01<RNGType> dist(generator);
	return dist();
}

void CNucleus::random_test(int nmax) {

	RNGType generator(time(0));
	printf("A = %i \nR = %lf fm\n", A_, R_);
	double r, theta, phi, x, y, z;
	for(int n = 0; n < nmax; n+=1) {

		// pick a random point in the nuclear space
		r = R_*random_01(generator);
		theta = pi*random_01(generator);
		phi = 2.0*pi*random_01(generator);
		
		// transform coords
		x = r*sin(theta)*cos(phi);
		y = r*sin(theta)*sin(phi);
		z = r*cos(theta);
		
		// print test values
		printf("rho(%.3lf, %.3lf, %.3lf)=%lf \t T(%.3lf, %.3lf)=%lf\n",
			   x,y,z,get_rho(x,y,z),x,y,get_T(x,y));
	}

}

void CNucleus::T_test(int nmax, string filename) {

	// open file 
	ofstream file;
	file.open(filename.c_str());

	// write thickness to file for various x and y
	double x, y, dx = 2.5*R_/nmax, dy = 2.5*R_/nmax;
	for(x = -1.25*R_; x < 1.25*R_; x += dx) {
		for(y = -1.25*R_; y < 1.25*R_; y += dy) {

			file << get_T(x,y) << "\t";
		}

		file << "\n";
	}

	// close file
	file.close();
}

class CPairs {
public:
	
	// declare member variables
	const static double sig_nn = 42.0; 
	CNucleus *N1_, *N2_;
	double dEdy_, sig_sat_, fwn_, b_;

	// declare member functions
	CPairs(CNucleus *N1_set, CNucleus *N2_set, double b_set, 
		double dEdy_set, double sig_sat_set, double fwn_set);
	double get_eps_wn(double x, double y);
	double get_eps_sat(double x, double y);
	double get_eps(double x, double y);
	void eps_test(int nmax, string filename);

};

CPairs::CPairs(CNucleus *N1_set, CNucleus *N2_set, double b_set, 
	double dEdy_set, double sig_sat_set, double fwn_set) {

	// set nucleii and impact parameter
	N1_ = N1_set;
	N2_ = N2_set;
	b_ = b_set;

	// set parameters
	dEdy_ = dEdy_set;
	sig_sat_ = sig_sat_set;
	fwn_ = fwn_set; 

}

double CPairs::get_eps_wn(double x, double y) {

	double prefactor, T1, T2, eps_wn;

	// N1 centered at (-b/2,0) and N2 centered at (b/2,0)
	prefactor = 0.5*dEdy_*sig_nn/sig_sat_;
	T1 = N1_->get_T(x+0.5*b_,y);
	T2 = N2_->get_T(x-0.5*b_,y);

	// calculate wounded-nucleon energy density
	eps_wn = prefactor*(T1*(1.0-exp(-T2*sig_sat_))+T2*(1.0-exp(-T1*sig_sat_)));
	
	return eps_wn;

}

double CPairs::get_eps_sat(double x, double y) {

	double prefactor, T1, T2, Tmin, Tmax, eps_sat;

	// N1 centered at (-b/2,0) and N2 centered at (b/2,0)
	prefactor = dEdy_*sig_nn/sig_sat_;
	T1 = N1_->get_T(x+0.5*b_,y);
	T2 = N2_->get_T(x-0.5*b_,y);
	Tmin = 2.0*T1*T2/(T1+T2);
	Tmax = 0.5*(T1+T2);

	// calculate saturation energy density
	eps_sat = prefactor*Tmin*(1.0-exp(-Tmax*sig_sat_));
	
	return eps_sat;

}

double CPairs::get_eps(double x, double y) {

	double eps;

	// calculate weighted average energy density
	eps = fwn_*get_eps_wn(x,y)+(1.0-fwn_)*get_eps_sat(x,y);

	return eps;
}

void CPairs::eps_test(int nmax, string filename) {

	// open file 
	ofstream file;
	file.open(filename.c_str());

	// calculate radii of nuclei
	double R1 = N1_->get_R();
	double R2 = N2_->get_R();

	// find larger radii
	double Rmax;
	if(R1>=R2) Rmax = R1;
	else Rmax = R2;


	// write energy density to file for various x and y
	double x, y, dx = (1.25*(R1+R2)+b_)/nmax, dy = 2.5*Rmax/nmax;
	for(x = -1.25*R1-0.5*b_; x < 1.25*R2+0.5*b_; x += dx) {
		for(y = -1.25*Rmax; y < 1.25*Rmax; y += dy) {

			file << get_eps(x,y) << "\t";
		}

		file << "\n";
	}

	// close file
	file.close();

}



#endif // GLAUBER_H