//glauber.cpp
#include "glauber.h"

namespace glauber {

	Nucleus::Nucleus(int A_set) {
	
		// Set the nucleon number
		A_ = A_set;

		// Set the Radius
		R_ = get_R();
	}
	
	double Nucleus::get_R() {

		double r, dr = 1.0E-3;
		double A, A_prime;
		double tolerance = 0.1;

		// Estimate R of uniform spherical nucleus
		R_ = pow(3.0*A_/(4.0*pi*rho0_), 1.0/3.0);

		// Perform Newton's method to find R s.t. A(R)=A_
		do {

			// Calculate A(R)
			A = 0.0;
			for(r = 0.0; r < 2.5*R_; r += dr) {
				A += 4.0*pi*get_rho(r)*r*r*dr;
			}

			// Calculate A'(R)
			A_prime = 4.0*pi*rho0_*R_*R_;
		

			// Update radius estimate
			R_ += (A_-A)/A_prime;

			//printf("R = %lf\tA = %lf\tA' = %lf\n", R_, A, A_prime);

		} while(fabs(A_-A) > tolerance);

		// Check if integrating gives expected nucleon number
		A = 0.0;
		for(r = 0.0; r < 2.5*R_; r+=dr) {
			A += 4.0*pi*get_rho(r)*r*r*dr;
		}
		//printf("R = %lf, A from integrating = %lf\n", R_, A);

		return R_;
	}
	
	double Nucleus::get_rho(double r) {

		// Calculate Woods-Saxon density
		double rho = rho0_/(1.0+exp((r-R_)/a_));

		return rho;
	}

	double Nucleus::get_rho(double x, double y, double z) {
	
		// calculate distance from center of nucleus
		double r = pow(x*x+y*y+z*z,0.5);
	
		// calculate Woods-Saxon density
		double rho = get_rho(r);

		return rho;
	}

	double Nucleus::get_T(double x, double y) {

		// initialize variables
		double z;
		double z_bound = 2.0*pow(R_*R_-x*x-y*y,0.5);
		double dz = 1.0E-3;
	
		// integrate the density in the z-direction to get thickness
		double T=0.0;
		for(z=-z_bound; z<z_bound; z+=dz) {
			T+=get_rho(x,y,z)*dz;
		}
	
		return T;
	}
	
	NucleusPair::NucleusPair(Nucleus *N1_set, Nucleus *N2_set, double b_set,
		double min_set, double max_set, int nmax_set) {

		// set nuclei and impact parameter
		N1_ = N1_set;
		N2_ = N2_set;
		b_ = b_set;

		// fix sigma_sat
		sigma_sat_ = 42.0;

		// initialize parameters and derivatives
		N_ = 1.0;
		fwn_ = 0.5;
		Dchi2_DN_ = 1.0;
		Dchi2_Dfwn_ = 1.0;

		// set x- & y-bounds and number of points for printing
		min_ = min_set;
		max_ = max_set;
		nmax_ = nmax_set;
		dx_ = (max_-min_)/nmax_;

		// resize newton's method increment and hessian matrix
		H_.resize(2);
		for(int i = 0; i < 2; ++i) H_[i].resize(2);

		// resize matrices to nmax x nmax
		jane_.resize(nmax_);
		liam_.resize(nmax_);
		T1_.resize(nmax_);
		T2_.resize(nmax_);
		epsilon_wn_.resize(nmax_);
		epsilon_sat_.resize(nmax_);
		for(int i = 0; i < nmax_; ++i) {
			jane_[i].resize(nmax_);
			liam_[i].resize(nmax_);
			T1_[i].resize(nmax_);
			T2_[i].resize(nmax_);
			epsilon_wn_[i].resize(nmax_);
			epsilon_sat_[i].resize(nmax_);
		}

		store_thickness();
		store_epsilons();
		fetch_liam();
	}

	void NucleusPair::minimize_chi2() {

		printf("%s\n", "Minimizing chi...");

		double tolerance = 1.0E-2;

		while( (fabs(Dchi2_DN_) > tolerance) || (fabs(Dchi2_Dfwn_ > tolerance)) ) {

			// fetch data
			fetch_jane();

			// calculate newton's method increment
			fill_increments();
			
			// apply newton's method
			N_ += delta_N_;
			fwn_ += delta_fwn_;

			printf("Chi2 = %.5e\tN = %.5e\tfwn = %.5e\tDchi2_DN = %.5e\tDchi2_Dfwn = %.5e\n",
				get_chi2(),N_,fwn_,Dchi2_DN_,Dchi2_Dfwn_);
		};

		printf("Done.\n");
	}


	double NucleusPair::get_chi2() {

		double chi = 0.0;
		
		for(i = 0; i < nmax_; ++i) {
			for(j = 0; j < nmax_; ++j) {
				chi += pow(liam_[i][j]-jane_[i][j],2.0);
			}
		}

		return chi;
	}

	void NucleusPair::fetch_jane() {

		double epsilon;

		// store and print grid of energy densities 
		for(y = min_; y < max_; y += dx_) {
			i = trunc((y-min_)/dx_);
			for(x = min_; x < max_; x += dx_) {
				j = trunc((x-min_)/dx_);
				jane_[i][j] = get_epsilon(x,y);
			}
		}
	}

	void NucleusPair::fetch_liam() {

		printf("%s\n", "Fetching data from 'average.dat'...");

		FILE *fptr;
		fptr = fopen("average.dat", "r");

		// store data in liam_ matrix
		for(i = 0; i < nmax_; ++i) {
			for(j = 0; j < nmax_; ++j) {
				fscanf(fptr, "%lf", &liam_[i][j]);
			}
		}

		fclose(fptr);
	}

	void NucleusPair::store_thickness() {

		printf("Storing thicknesses...\n");

		for(y = min_; y < max_; y += dx_) {
			i = trunc((y-min_)/dx_);
			for(x = min_; x < max_; x += dx_) {
				j = trunc((x-min_)/dx_);
				T1_[i][j] = N1_->get_T(x+0.5*b_,y);
				T2_[i][j] = N2_->get_T(x-0.5*b_,y);
			}
		}
	}

	void NucleusPair::store_epsilons() {

		printf("Storing energy densities...\n");

		for(y = min_; y < max_; y += dx_) {
			i = trunc((y-min_)/dx_);
			for(x = min_; x < max_; x += dx_) {
				j = trunc((x-min_)/dx_);
				epsilon_wn_[i][j] = get_epsilon_wn(x,y);
				epsilon_sat_[i][j] = get_epsilon_sat(x,y);
			}
		}
	}

	void NucleusPair::fill_H() {

		double a, b, d, det;

		a = D2chi2_DN2();
		b = D2chi2_DN_Dfwn();
		d = D2chi2_Dfwn2();
		det = a*d-b*b;

		H_[0][0] = -d/det;
		H_[0][1] = H_[1][0] = b/det;
		H_[1][1] = -a/det;

	}

	void NucleusPair::fill_increments() {

		Dchi2_DN();
		Dchi2_Dfwn();
		fill_H();
		delta_N_ = H_[0][0]*Dchi2_DN_+H_[0][1]*Dchi2_Dfwn_;
		delta_fwn_ = H_[1][0]*Dchi2_DN_+H_[1][1]*Dchi2_Dfwn_;
	}

	double NucleusPair::get_epsilon(double x, double y) {
	
		double epsilon;

		// calculate weighted average energy density
		epsilon = N_*(fwn_*get_epsilon_wn(x,y)+(1.0-fwn_)*get_epsilon_sat(x,y));

		return epsilon;
	}

	double NucleusPair::get_epsilon_wn(double x, double y) {

		double T1, T2, epsilon_wn;

		// N1 centered at (-b/2,0) and N2 centered at (b/2,0)
		i = trunc((y-min_)/dx_);
		j = trunc((x-min_)/dx_);

		T1 = T1_[i][j];
		T2 = T2_[i][j];

		// calculate wounded-nucleon energy density
		epsilon_wn = (0.5/sigma_sat_)*(T1*(1.0-exp(-T2*sigma_sat_))+T2*(1.0-exp(-T1*sigma_sat_)));

		return epsilon_wn;
	}

	double NucleusPair::get_epsilon_sat(double x, double y) {

		double T1, T2, Tmin, Tmax, epsilon_sat;

		// N1 centered at (-b/2,0) and N2 centered at (b/2,0)
		i = trunc((y-min_)/dx_);
		j = trunc((x-min_)/dx_);

		T1 = T1_[i][j];
		T2 = T2_[i][j];

		if(T1+T2==0.0) Tmin = 0.0;
		else Tmin = 2.0*T1*T2/(T1+T2);
		Tmax = 0.5*(T1+T2);

		// calculate saturation energy density
		epsilon_sat = (1.0/sigma_sat_)*Tmin*(1.0-exp(-Tmax*sigma_sat_));
	
		return epsilon_sat;
	}

	void NucleusPair::Dchi2_DN() {

		Dchi2_DN_ = 0.0;

		for(i = 0; i < nmax_; ++i) {
			for(j = 0; j < nmax_; ++j) {
				Dchi2_DN_ += 2.0*(jane_[i][j]-liam_[i][j])*
					(fwn_*epsilon_wn_[i][j]+(1-fwn_)*epsilon_sat_[i][j]);			
			}
		}
	}

	void NucleusPair::Dchi2_Dfwn() {

		Dchi2_Dfwn_ = 0.0;
		
		for(i = 0; i < nmax_; ++i) {
			for(j = 0; j < nmax_; ++j) {
				Dchi2_Dfwn_ += 2.0*N_*(jane_[i][j]-liam_[i][j])*
					(epsilon_wn_[i][j]-epsilon_sat_[i][j]);		
			}
		}
	}

	double NucleusPair::D2chi2_DN2() {

		double D2chi2_DN2 = 0.0;

		for(i = 0; i < nmax_; ++i) {
			for(j = 0; j < nmax_; ++j) {
				D2chi2_DN2 += 2.0*pow(fwn_*epsilon_wn_[i][j]+(1-fwn_)*epsilon_sat_[i][j],2.0);
			}
		}

		return D2chi2_DN2;
	}

	double NucleusPair::D2chi2_Dfwn2() {

		double D2chi2_Dfwn2 = 0.0;

		for(i = 0; i < nmax_; ++i) {
			for(j = 0; j < nmax_; ++j) {
				D2chi2_Dfwn2 += 2.0*N_*N_*pow(epsilon_wn_[i][j]-epsilon_sat_[i][j],2.0);
			}
		}

		return D2chi2_Dfwn2;
	}

	// same as D2chi2_Dfwn_DN()
	double NucleusPair::NucleusPair::D2chi2_DN_Dfwn() {

		double epsilon_wn, epsilon_sat, D2chi2_DN_Dfwn = 0.0;
		
		for(i = 0; i < nmax_; ++i) {
			for(j = 0; j < nmax_; ++j) {
				epsilon_wn = epsilon_wn_[i][j];
				epsilon_sat = epsilon_sat_[i][j];
				D2chi2_DN_Dfwn += -2.0*(epsilon_wn-epsilon_sat)*
					(N_*(fwn_*epsilon_wn+(1-fwn_)*epsilon_sat)+jane_[i][j]-liam_[i][j]);
			}
		}

		return D2chi2_DN_Dfwn;
	}



} //glauber
