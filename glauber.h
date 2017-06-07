#ifndef GLAUBER_H
#define GLAUBER_H
#include <cmath>
#include <cstdlib>
#include <fstream>


const double pi = 4.0*atan(1.0);
using namespace std;

class CNucleus{
public:
	
	const static double a_ = 0.54;
	const static double rho0_ = 0.17;
	double R_;
	int A_;

	CNucleus(int Aset);
	double get_R();
	double get_rho(double r);
	double get_rho(double x, double y, double z);
	double get_T(double x, double y);

};

CNucleus::CNucleus(int A_set){
	
	// set the nucleon number
	A_ = A_set;

	// find radius
	R_ = get_R();
}

double CNucleus::get_R(){

	double r, dr = 1.0E-3;
	double A, A_prime;
	double tolerance = 0.1;

	// estimate R of uniform spherical nucleus
	R_ = pow(3.0*A_/(4.0*pi*rho0_), 1.0/3.0);

	// perform Newton's method to find R s.t. A(R)=A_
	do{

		// calculate A(R)
		A = 0.0;
		for(r = 0.0; r < R_; r += dr){
			A += 4.0*pi*get_rho(r)*r*r*dr;
		}

		// calculate A'(R)
		A_prime = 4.0*pi*rho0_*R_*R_;

		// update radius estimate
		R_ += (A_-A)/A_prime;

		//printf("R = %lf\tA = %lf\tA' = %lf\n", R_, A, A_prime);

	} while(abs(A_-A) > tolerance);

	/* 
	// check if integrating gives expected nucleon number
	A = 0.0;
	for(r = 0.0; r < R_; r+=dr){
		A += 4.0*pi*get_rho(r)*r*r*dr;
	}

	printf("R = %lf, A from integrating = %lf\n", R_, A);
	*/

	return R_;
}

double CNucleus::get_rho(double r){

	// calculate Woods-Saxon density
	double rho = rho0_/(1.0+exp((r-R_)/a_));

	return rho;
}

double CNucleus::get_rho(double x, double y, double z){
	
	// calculate distance from center of nucleus
	double r = pow(x*x+y*y+z*z,0.5);
	
	// calculate Woods-Saxon density
	double rho = get_rho(r);

	return rho;
}

double CNucleus::get_T(double x, double y){

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
	
	const static double sig_nn = 42.0; 
	CNucleus *N1_, *N2_;
	double dEdy_, sig_sat_, fwn_, b_;

	// declare member functions
	CPairs(CNucleus *N1_set, CNucleus *N2_set, double b_set, 
	double dEdy_set, double sig_sat_set, double fwn_set);
	double get_eps_wn(double x, double y);
	double get_eps_sat(double x, double y);
	double get_eps(double x, double y);

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
	
	// calculate weighted average energy density
	double eps = fwn_*get_eps_wn(x,y)+(1.0-fwn_)*get_eps_sat(x,y);

	return eps;
}

#endif