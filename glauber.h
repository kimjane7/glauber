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

} // glauber

#endif // GLAUBER_H