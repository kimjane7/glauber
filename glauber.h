#ifndef GLAUBER_H
#define GLAUBER_H
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

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

	} while(abs(A_-A) > tolerance);

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
	int nmax_;
	double norm_, b_, min_, max_;

	std::vector<double> params_, param_step_, Dchi_;
	std::vector<std::vector<double>> liam_, jane_, H_;

	CPairs(CNucleus *N1_set, CNucleus *N2_set, double b_set,
		double min_set, double max_set, int nmax_set);
	double normalize();
	double get_eps(double x, double y);
	double get_epswn(double x, double y);
	double get_epssat(double x, double y);

	void print_eps();
	void fetch_liam();
	void minimize_chi();
	void fill_H();
	void fill_Dchi();
	void fill_step();

	double get_chi();
	double Dchi_Dfwn();
	double Dchi_Dsigsat();
	double D2chi_Dfwn2();
	double D2chi_Dsigsat2();
	double D2chi_Dsigsat_Dfwn();
	double D2chi_Dfwn_Dsigsat();

	double Deps_Dfwn(double x, double y);
	double Deps_Dsigsat(double x, double y);
	double D2eps_Dsigsat2(double x, double y);
	double D2eps_Dsigsat_Dfwn(double x, double y);

	double Depswn_Dsigsat(double x, double y);
	double D2epswn_Dsigsat2(double x, double y);

	double Depssat_Dsigsat(double x, double y);
	double D2epssat_Dsigsat2(double x, double y);

};

CPairs::CPairs(CNucleus *N1_set, CNucleus *N2_set, double b_set,
	double min_set, double max_set, int nmax_set) {

	// set nucleii and impact parameter
	N1_ = N1_set;
	N2_ = N2_set;
	b_ = b_set;

	// initialize parameter vector (fwn, sig_sat)
	params_.resize(2);
	params_[0] = 0.5;
	params_[1] = 42.0;

	// calculate normalization constant
	norm_ = normalize();

	// set x- & y-bounds and number of points for printing
	min_ = min_set;
	max_ = max_set;
	nmax_ = nmax_set;

	// resize newton's method increment, partial chi's, jacobian matrix
	param_step_.resize(2);
	Dchi_.resize(2);
	H_.resize(2);
	for(int i = 0; i < 2; i++) H_[i].resize(2);

	// resize data matrices to nmax x nmax
	jane_.resize(nmax_);
	liam_.resize(nmax_);
	for(int i = 0; i < nmax_; i++) {
		jane_[i].resize(nmax_);
		liam_[i].resize(nmax_);
	}

	fetch_liam();

}

double CPairs::normalize() {

	printf("Normalizing energy density...\n");

	double x, xmin, xmax, dx = 0.1, y, ymax, dy = dx;
	double f = 0.005, fwn = params_[0], R1, R2, Rext, eps, sum = 0.0;

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

			eps = fwn*get_epswn(x,y)+(1.0-fwn)*get_epssat(x,y);
			sum += 2.0*eps*dx*dy;
		}	

		//printf("sum = %lf\n", sum);	
	}

	norm_ = 1.0/sum;

	printf("Normalization constant = %lf.\n", norm_);

	return norm_;
}

double CPairs::get_epswn(double x, double y) {

	double prefactor, T1, T2, sigsat = params_[1], epswn;

	// N1 centered at (-b/2,0) and N2 centered at (b/2,0)
	T1 = N1_->get_T(x+0.5*b_,y);
	T2 = N2_->get_T(x-0.5*b_,y);

	// calculate wounded-nucleon energy density
	epswn = (0.5/sigsat)*(T1*(1.0-exp(-T2*sigsat))+T2*(1.0-exp(-T1*sigsat)));

	return epswn;
}

double CPairs::get_epssat(double x, double y) {

	double T1, T2, Tmin, Tmax, sigsat = params_[1], epssat;

	// N1 centered at (-b/2,0) and N2 centered at (b/2,0)
	T1 = N1_->get_T(x+0.5*b_,y);
	T2 = N2_->get_T(x-0.5*b_,y);

	if(T1+T2==0.0) Tmin = 0.0;
	else Tmin = 2.0*T1*T2/(T1+T2);
	Tmax = 0.5*(T1+T2);

	// calculate saturation energy density
	epssat = (1.0/sigsat)*Tmin*(1.0-exp(-Tmax*sigsat));
	
	return epssat;
}

double CPairs::get_eps(double x, double y) {
	
	double fwn = params_[0], eps;

	// calculate weighted average energy density
	eps = norm_*(fwn*get_epswn(x,y)+(1.0-fwn)*get_epssat(x,y));

	return eps;
}

// also fetches jane_
void CPairs::print_eps() {

	printf("%s\n", "Printing energy densities to 'energy_density.dat'...");

	double x, y, dx = (max_-min_)/nmax_, dy = dx, eps;
	int i, j;

	// open file
	FILE *fptr;
	fptr = fopen("energy_density.dat", "w");

	// store and print grid of energy densities 
	for(y = min_; y < max_; y += dy) {
		i = trunc((y-min_)/dy);
		for(x = min_; x < max_; x += dx) {
			j = trunc((x-min_)/dx);
			eps = get_eps(x,y);
			jane_[i][j] = eps;
			fprintf(fptr, "%lf\t", eps);
		}
		fprintf(fptr, "\n");
	}

	fclose(fptr);
	printf("Done.\n");
}

void CPairs::fetch_liam() {

	printf("%s\n", "Fetching data from 'average.dat'...");

	int i, j;

	FILE *fptr;
	fptr = fopen("average.dat", "r");

	// store data in liam_ matrix
	for(i = 0; i < nmax_; i++) {
		for(j = 0; j < nmax_; j++) {
			fscanf(fptr, "%lf", &liam_[i][j]);
		}
	}

	fclose(fptr);
	printf("Done.\n");
}

void CPairs::minimize_chi() {

	printf("%s\n", "Minimizing chi...");

	double tolerance = 1.0E-2;
	do {
		print_eps();
		printf("Chi = %lf\n",get_chi());
		fill_step();
		printf("CHECK 1\n");
		params_[0] += param_step_[0];
		params_[1] += param_step_[1];
	} while( (fabs(Dchi_[0]) > tolerance) && (fabs(Dchi_[0]) > tolerance) );
	

	printf("Done.\n");

}


// H is the negative inverse of the Hessian of chi
void CPairs::fill_H() { 

	double a, c, d, det;

	// b = c
	printf("CHECK 8\n");
	a = D2chi_Dfwn2();
	c = D2chi_Dsigsat_Dfwn();
	d = D2chi_Dsigsat2();
	det = a*d-c*c;

	H_[0][0] = -(1.0/det)*d;
	H_[0][1] = H_[1][0] = (1.0/det)*c;
	H_[1][1] = -(1.0/det)*a;
}

void CPairs::fill_Dchi() {

	printf("CHECK 6\n");
	Dchi_[0] = Dchi_Dfwn();
	printf("CHECK 7\n");
	Dchi_[1] = Dchi_Dsigsat();
}

void CPairs::fill_step() {

	printf("CHECK 2\n");
	fill_Dchi();
	printf("CHECK 3\n");
	fill_H();
	printf("CHECK 4\n");
	param_step_[0] = H_[0][0]*Dchi_[0]+H_[0][1]*Dchi_[1];
	param_step_[1] = H_[1][0]*Dchi_[0]+H_[1][1]*Dchi_[1];
	printf("CHECK 5\n");
}

double CPairs::get_chi() {

	double chi = 0.0;
	int i, j;

	// iterate over data points in 'average.dat' and 'energy_density.dat'
	for(i = 0; i < nmax_; i++) {
		for(j = 0; j < nmax_; j++) {
			chi += (liam_[i][j]-jane_[i][j])*(liam_[i][j]-jane_[i][j]);
		}
	}

	return chi;
}



double CPairs::Dchi_Dfwn() {

	double x, y, dx = (max_-min_)/nmax_, dy = dx, Dchi_Dfwn = 0.0;
	int i, j;

	// iterate over data points in 'average.dat' and 'energy_density.dat'
	for(i = 0; i < nmax_; i++) {
		y = min_+i*dy;
		for(j = 0; j < nmax_; j++) {
			x = min_+i*dx;
			Dchi_Dfwn += -2.0*(liam_[i][j]-jane_[i][j])*Deps_Dfwn(x,y);
			// printf("Dchi_Dfwn = %lf\n",Dchi_Dfwn);
		}
	}

	return Dchi_Dfwn;
}

double CPairs::Dchi_Dsigsat() {

	double x, y, dx = (max_-min_)/nmax_, dy = dx, Dchi_Dsigsat = 0.0;
	int i, j;

	// iterate over data points in 'average.dat' and 'energy_density.dat'
	for(i = 0; i < nmax_; i++) {
		y = min_+i*dy;
		for(j = 0; j < nmax_; j++) {
			x = min_+i*dx;
			Dchi_Dsigsat += -2.0*(liam_[i][j]-jane_[i][j])*Deps_Dsigsat(x,y);
			// printf("Dchi_Dsigsat = %lf\n",Dchi_Dsigsat);
		}
	}

	return Dchi_Dsigsat;
}

double CPairs::D2chi_Dfwn2() { 

	double x, y, dx = (max_-min_)/nmax_, dy = dx, D2chi_Dfwn2 = 0.0;

	for(y = min_; y < max_; y += dy) {
		for(x = min_; x < max_; x += dx) {
			D2chi_Dfwn2 += 2.0*pow(Deps_Dfwn(x,y),2.0);
		}
	}

	return D2chi_Dfwn2;
}

double CPairs::D2chi_Dsigsat2() {

	double x, y, dx = (max_-min_)/nmax_, dy = dx, D2chi_Dsigsat2 = 0.0;
	int i, j;

	// iterate over data points in 'average.dat' and 'energy_density.dat'
	for(i = 0; i < nmax_; i++) {
		y = min_+i*dy;
		for(j = 0; j < nmax_; j++) {
			x = min_+i*dx;
			D2chi_Dsigsat2 += 2.0*(pow(Deps_Dsigsat(x,y),2.0)-(liam_[i][j]-jane_[i][j])*D2eps_Dsigsat_Dfwn(x,y));		
		}
	}

	return D2chi_Dsigsat2;
}

// same as D2chi_Dfwn_Dsigsat
double CPairs::D2chi_Dsigsat_Dfwn() { 

	double x, y, dx = (max_-min_)/nmax_, dy = dx, D2chi_Dsigsat_Dfwn = 0.0;
	int i, j;

	// iterate over data points in 'average.dat' and 'energy_density.dat'
	for(i = 0; i < nmax_; i++) {
		y = min_+i*dy;
		for(j = 0; j < nmax_; j++) {
			x = min_+i*dx;
			D2chi_Dsigsat_Dfwn += 2.0*(Deps_Dfwn(x,y)*Deps_Dsigsat(x,y)-(liam_[i][j]-jane_[i][j])*D2eps_Dsigsat_Dfwn(x,y));		
		}
	}

	return D2chi_Dsigsat_Dfwn;
}

double CPairs::Deps_Dfwn(double x, double y) {

	double Deps_Dfwn = get_epswn(x,y)-get_epssat(x,y);

	return Deps_Dfwn;
}

double CPairs::Deps_Dsigsat(double x, double y) {

	double fwn = params_[0], Deps_Dsigsat;

	Deps_Dsigsat = fwn*Depswn_Dsigsat(x,y)+(1-fwn)*Depssat_Dsigsat(x,y);

	return Deps_Dsigsat;
}

double CPairs::D2eps_Dsigsat2(double x, double y) {

	double fwn = params_[0], D2eps_Dsigsat2;

	D2eps_Dsigsat2 = fwn*D2epswn_Dsigsat2(x,y)+(1-fwn)*D2epssat_Dsigsat2(x,y);

	return D2eps_Dsigsat2;
}

// same as D2eps_Dfwn_Dsigsat
double CPairs::D2eps_Dsigsat_Dfwn(double x, double y) {

	double D2eps_Dsigsat_Dfwn = Depswn_Dsigsat(x,y)-Depssat_Dsigsat(x,y);

	return D2eps_Dsigsat_Dfwn;
}

double CPairs::Depswn_Dsigsat(double x, double y) {

	double T1, T2, sigsat = params_[1], Depswn_Dsigsat;

	T1 = N1_->get_T(x+0.5*b_,y);
	T2 = N2_->get_T(x-0.5*b_,y);

	Depswn_Dsigsat = (0.5*norm_/sigsat)*T1*T2*(exp(-T1*sigsat)+exp(-T2*sigsat));
	Depswn_Dsigsat += -(0.5*norm_*T1/(sigsat*sigsat))*(1.0-exp(-T2*sigsat));
	Depswn_Dsigsat += -(0.5*norm_*T2/(sigsat*sigsat))*(1.0-exp(-T1*sigsat));

	return Depswn_Dsigsat;
}

double CPairs::D2epswn_Dsigsat2(double x, double y) {

	double T1, T2, sigsat = params_[1], D2epswn_Dsigsat2;

	T1 = N1_->get_T(x+0.5*b_,y);
	T2 = N2_->get_T(x-0.5*b_,y);

	D2epswn_Dsigsat2 = (norm_*T1/pow(sigsat,3.0))*(1.0-exp(-T2*sigsat));
	D2epswn_Dsigsat2 += (norm_*T2/pow(sigsat,3.0))*(1.0-exp(-T1*sigsat));
	D2epswn_Dsigsat2 += (norm_*T1*T2/(sigsat*sigsat))*(exp(-T1*sigsat)+exp(-T2*sigsat));
	D2epswn_Dsigsat2 += -(0.5*norm_*T1*T2/sigsat)*(T1*exp(-T1*sigsat)+T2*exp(-T2*sigsat));

	return D2epswn_Dsigsat2;
}

double CPairs::Depssat_Dsigsat(double x, double y) {

	double T1, T2, Tmin, Tmax, sigsat = params_[1], Depssat_Dsigsat;

	// N1 centered at (-b/2,0) and N2 centered at (b/2,0)
	T1 = N1_->get_T(x+0.5*b_,y);
	T2 = N2_->get_T(x-0.5*b_,y);

	if(T1+T2==0.0) Tmin = 0.0;
	else Tmin = 2.0*T1*T2/(T1+T2);
	Tmax = 0.5*(T1+T2);

	Depssat_Dsigsat = (norm_/sigsat)*Tmin*Tmax*exp(-Tmax*sigsat);
	Depssat_Dsigsat += -(norm_/(sigsat*sigsat))*Tmin*(1.0-exp(-Tmax*sigsat));

	return Depssat_Dsigsat;
}

double CPairs::D2epssat_Dsigsat2(double x, double y) {
	
	double T1, T2, Tmin, Tmax, sigsat = params_[1], D2epssat_Dsigsat2;

	// N1 centered at (-b/2,0) and N2 centered at (b/2,0)
	T1 = N1_->get_T(x+0.5*b_,y);
	T2 = N2_->get_T(x-0.5*b_,y);

	if(T1+T2==0.0) Tmin = 0.0;
	else Tmin = 2.0*T1*T2/(T1+T2);
	Tmax = 0.5*(T1+T2);

	D2epssat_Dsigsat2 = (2.0*norm_*Tmin/pow(sigsat,3.0))*(1.0-exp(-Tmax*sigsat));
	D2epssat_Dsigsat2 += -(2.0*norm_*Tmin*Tmax/(sigsat*sigsat))*exp(-Tmax*sigsat);
	D2epssat_Dsigsat2 += -(norm_*Tmin*Tmax*Tmax/sigsat)*exp(-Tmax*sigsat);

	return D2epssat_Dsigsat2;
}

} // glauber

#endif