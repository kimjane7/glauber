#ifndef GLAUBER_H
#define GLAUBER_H
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

const double pi = 4.0*atan(1.0);

namespace glauber {

class Nucleus {
public:
	
	constexpr static double a_ = 0.546; // only constexpr, const enum, and const int values
	constexpr static double rho0_ = 0.159; // are allowed to be set in class definition
	double R_;
	int A_;

	Nucleus(int Aset);
	double get_R();
	double get_rho(double r);
	double get_rho(double x, double y, double z);
	double get_T(double x, double y);

};

class NucleusPair {
public:
	
	Nucleus *N1_, *N2_;
	int nmax_;
	double norm_, b_, min_, max_;

	std::vector<double> params_, param_step_, Dchi_;
	std::vector<std::vector<double>> liam_, jane_, H_;

	NucleusPair(Nucleus *N1_set, Nucleus *N2_set, double b_set,
		double min_set, double max_set, int nmax_set);
	
	double normalize();
	double get_epsilon(double x, double y);
	double get_epsilon_wn(double x, double y);
	double get_epsilon_sat(double x, double y);

	void print_epsilon();
	void fetch_liam();
	void minimize_chi();
	void fill_H();
	void fill_Dchi();
	void fill_step();

	double get_chi();
	double Dchi_Df_wn(); //
	double Dchi_Dsigma_sat();
	double D2chi_Df_wn2();
	double D2chi_Dsigma_sat2();
	double D2chi_Dsigma_sat_Df_wn();
	double D2chi_Df_wn_Dsigma_sat();

	double Depsilon_Df_wn(double x, double y);
	double Depsilon_Dsigma_sat(double x, double y);
	double D2epsilon_Dsigma_sat2(double x, double y);
	double D2epsilon_Dsigma_sat_Df_wn(double x, double y);

	double Depsilon_wn_Dsigma_sat(double x, double y);
	double D2epsilon_wn_Dsigma_sat2(double x, double y);

	double Depsilon_sat_Dsigma_sat(double x, double y);
	double D2epsilon_sat_Dsigma_sat2(double x, double y);

};

} // glauber

#endif // GLAUBER_H