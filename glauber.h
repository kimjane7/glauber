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
	int i, j, nmax_;
	double x, y, dx_, min_, max_, b_;
	double N_, fwn_, sigma_sat_;
	double delta_N_, delta_fwn_;
	double Dchi2_DN_, Dchi2_Dfwn_;

	std::vector<std::vector<double>> liam_, jane_;
	std::vector<std::vector<double>> T1_, T2_, H_;
	std::vector<std::vector<double>> epsilon_wn_, epsilon_sat_;

	NucleusPair(Nucleus *N1_set, Nucleus *N2_set, double b_set,
		double min_set, double max_set, int nmax_set);

	void minimize_chi2();
	double get_chi2();

private:
	
	void fetch_jane();
	void fetch_liam();
	void store_thickness();
	void store_epsilons();
	void fill_H();
	void fill_increments();

	double get_epsilon(double x, double y);
	double get_epsilon_wn(double x, double y);
	double get_epsilon_sat(double x, double y);

	void Dchi2_DN();
	void Dchi2_Dfwn();
	double D2chi2_DN2();
	double D2chi2_Dfwn2();
	double D2chi2_DN_Dfwn();


};

} // glauber

#endif // GLAUBER_H