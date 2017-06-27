#ifndef GLAUBER_H
#define GLAUBER_H
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>


const double pi = 4.0*atan(1.0);

namespace glauber {

class CNucleus {
public:
	
	constexpr static double a_ = 0.546; // only constexpr, const enum, and const int values
	constexpr static double rho0_ = 0.159; // are allowed to be set in class definition
	double R_;
	int A_;

	CNucleus(int Aset);
	double get_R();
	double get_rho(double r);
	double get_rho(double x, double y, double z);
	double get_T(double x, double y);

};

CNucleus::CNucleus(int A_set) {
	
	// set the nucleon number
	A_ = A_set;

	// find radius
	R_ = get_R();
}

double CNucleus::get_R() {

	double r, dr = 1.0E-3;
	double A, A_prime;
	double tolerance = 0.1;

	// estimate R of uniform spherical nucleus
	R_ = pow(3.0*A_/(4.0*pi*rho0_), 1.0/3.0);

	// perform Newton's method to find R s.t. A(R)=A_
	do {

		// calculate A(R)
		A = 0.0;
		for(r = 0.0; r < 2.5*R_; r += dr) {
			A += 4.0*pi*get_rho(r)*r*r*dr;
		}

		// calculate A'(R)
		A_prime = 4.0*pi*rho0_*R_*R_;
		

		// update radius estimate
		R_ += (A_-A)/A_prime;

		//printf("R = %lf\tA = %lf\tA' = %lf\n", R_, A, A_prime);

	} while(fabs(A_-A) > tolerance);

	// check if integrating gives expected nucleon number
	A = 0.0;
	for(r = 0.0; r < 2.5*R_; r+=dr) {
		A += 4.0*pi*get_rho(r)*r*r*dr;
	}
	//printf("R = %lf, A from integrating = %lf\n", R_, A);

	return R_;
}

double CNucleus::get_rho(double r) {

	// calculate Woods-Saxon density
	double rho = rho0_/(1.0+exp((r-R_)/a_));

	return rho;
}

double CNucleus::get_rho(double x, double y, double z) {
	
	// calculate distance from center of nucleus
	double r = pow(x*x+y*y+z*z,0.5);
	
	// calculate Woods-Saxon density
	double rho = get_rho(r);

	return rho;
}

double CNucleus::get_T(double x, double y) {

	// initialize variables
	double z, z_bound = 2.0*pow(R_*R_-x*x-y*y,0.5), dz = 1.0E-3;
	
	// integrate the density in the z-direction to get thickness
	double T=0.0;
	for(z=-z_bound; z<z_bound; z+=dz) {
		T+=get_rho(x,y,z)*dz;
	}
	
	return T;
}

class CPairs {
public:
	
	CNucleus *N1_, *N2_;
	double norm_, sig_sat_, fwn_, b_;

	CPairs(CNucleus *N1_set, CNucleus *N2_set, double b_set,
		double sig_sat_set, double fwn_set);
	double normalize();
	double get_eps_wn(double x, double y);
	double get_eps_sat(double x, double y);
	double get_eps(double x, double y);
	void print_eps(double min, double max, int nmax);

};

CPairs::CPairs(CNucleus *N1_set, CNucleus *N2_set, double b_set,
		double sig_sat_set, double fwn_set) {

	// set nucleii and impact parameter
	N1_ = N1_set;
	N2_ = N2_set;
	b_ = b_set;

	// set parameters
	sig_sat_ = sig_sat_set;
	fwn_ = fwn_set;

	// calculate normalization constant
	norm_ = normalize();
}

double CPairs::normalize() {

	printf("Normalizing energy density...\n");

	double x, xmin, xmax, dx = 0.1, y, ymax, dy = dx;
	double f = 0.005, R1, R2, Rext, eps, sum = 0.0;

	// extension to radius such that rho(R_+Rext)=f*rho0;
	Rext = CNucleus::a_*log((1.0/f)-1.0);

	// find maximum radii to integrate over
	R1 = N1_->get_R()+Rext;
	R2 = N2_->get_R()+Rext;

	// calculate upper y-bound
	ymax = R1*sin(acos((R2*R2-R1*R1-b_*b_)/(2.0*R1*b_)));

	// integrate over half of overlapping region, double the running sum
	for(y = 0; y < ymax; y += dy) {

		// calculate x-bounds
		xmin = 0.5*b_-sqrt(R2*R2-y*y);
		xmax = sqrt(R1*R1-y*y)-0.5*b_;

		// integrate
		for(x = xmin; x < xmax; x += dx){

			eps = fwn_*get_eps_wn(x,y)+(1.0-fwn_)*get_eps_sat(x,y);
			sum += 2.0*eps*dx*dy;
		}	

		//printf("sum = %lf\n", sum);	
	}

	norm_ = 1.0/sum;

	printf("Normalization constant = %lf.\n", norm_);

	return norm_;
}

double CPairs::get_eps_wn(double x, double y) {

	double prefactor, T1, T2, eps_wn;

	// N1 centered at (-b/2,0) and N2 centered at (b/2,0)
	T1 = N1_->get_T(x+0.5*b_,y);
	T2 = N2_->get_T(x-0.5*b_,y);

	// calculate wounded-nucleon energy density
	eps_wn = (0.5/sig_sat_)*(T1*(1.0-exp(-T2*sig_sat_))
		+ T2*(1.0-exp(-T1*sig_sat_)));

	return eps_wn;
}

double CPairs::get_eps_sat(double x, double y) {

	double T1, T2, Tmin, Tmax, eps_sat;

	// N1 centered at (-b/2,0) and N2 centered at (b/2,0)
	T1 = N1_->get_T(x+0.5*b_,y);
	T2 = N2_->get_T(x-0.5*b_,y);

	if(T1+T2==0.0) Tmin = 0.0;
	else Tmin = 2.0*T1*T2/(T1+T2);
	Tmax = 0.5*(T1+T2);

	// calculate saturation energy density
	eps_sat = (1.0/sig_sat_)*Tmin*(1.0-exp(-Tmax*sig_sat_));

	//printf("T1 = %lf\tT2 = %lf\tTmin = %lf\tTmax = %lf\t eps_sat = %lf\n", T2,T2,Tmin,Tmax,eps_sat);
	
	return eps_sat;
}

double CPairs::get_eps(double x, double y) {
	
	// calculate weighted average energy density
	double eps = norm_*(fwn_*get_eps_wn(x,y)+(1.0-fwn_)*get_eps_sat(x,y));

	return eps;
}

void CPairs::print_eps(double min, double max, int nmax) {

	printf("%s\n", "Printing energy densities to 'energy_density.dat'...");

	double x, y, dx = (max-min)/nmax, dy = dx;

	// open file
	std::ofstream file;
	file.open("energy_density.dat");

	// print grid of energy densities
	for(y = min; y < max; y += dy) {
		for(x = min; x < max; x += dx) {
			file << get_eps(x,y) << "\t";
		}
		file << "\n";
	}

	printf("%s\n", "Done.");
}

} // glauber

#endif