//glauber.cpp
#include "glauber.h"

namespace glauber {

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
	
	CPairs::CPairs(CNucleus *N1_set, CNucleus *N2_set, double b_set,
		double min_set, double max_set, int nmax_set) {

		// set nuclei and impact parameter
		N1_ = N1_set;
		N2_ = N2_set;
		b_ = b_set;

		// initialize parameter vector (fwn, sigma_sat)
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
		double f = 0.005, fwn = params_[0], R1, R2, Rext, epsilon, sum = 0.0;

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

				epsilon = fwn*get_epsilon_wn(x,y)+(1.0-fwn)*get_epsilon_sat(x,y);
				sum += 2.0*epsilon*dx*dy;
			}	

			//printf("sum = %lf\n", sum);	
		}

		norm_ = 1.0/sum;

		printf("Normalization constant = %lf.\n", norm_);

		return norm_;
	}
	
	double CPairs::get_epsilon_wn(double x, double y) {

		double prefactor, T1, T2, sigma_sat = params_[1], epsilon_wn;

		// N1 centered at (-b/2,0) and N2 centered at (b/2,0)
		T1 = N1_->get_T(x+0.5*b_,y);
		T2 = N2_->get_T(x-0.5*b_,y);

		// calculate wounded-nucleon energy density
		epsilon_wn = (0.5/sigma_sat)*(T1*(1.0-exp(-T2*sigma_sat))+T2*(1.0-exp(-T1*sigma_sat)));

		return epsilon_wn;
	}

	double CPairs::get_epsilon_sat(double x, double y) {

		double T1, T2, Tmin, Tmax, sigma_sat = params_[1], epsilon_sat;

		// N1 centered at (-b/2,0) and N2 centered at (b/2,0)
		T1 = N1_->get_T(x+0.5*b_,y);
		T2 = N2_->get_T(x-0.5*b_,y);

		if(T1+T2==0.0) Tmin = 0.0;
		else Tmin = 2.0*T1*T2/(T1+T2);
		Tmax = 0.5*(T1+T2);

		// calculate saturation energy density
		epsilon_sat = (1.0/sigma_sat)*Tmin*(1.0-exp(-Tmax*sigma_sat));
	
		return epsilon_sat;
	}

	double CPairs::get_epsilon(double x, double y) {
	
		double fwn = params_[0], epsilon;

		// calculate weighted average energy density
		epsilon = norm_*(fwn*get_epsilon_wn(x,y)+(1.0-fwn)*get_epsilon_sat(x,y));

		return epsilon;
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
			print_epsilon();
			printf("Chi = %lf\n",get_chi());
			fill_step();
			printf("CHECK 1\n");
			params_[0] += param_step_[0];
			params_[1] += param_step_[1];
		} while( (fabs(Dchi_[0]) > tolerance) && (fabs(Dchi_[0]) > tolerance) );
	

		printf("Done.\n");

	}
	
	// also fetches jane_
	void CPairs::print_epsilon() {

		printf("%s\n", "Printing energy densities to 'energy_density.dat'...");

		double x, y, dx = (max_-min_)/nmax_, dy = dx, epsilon;
		int i, j;

		// open file
		FILE *fptr;
		fptr = fopen("energy_density.dat", "w");

		// store and print grid of energy densities 
		for(y = min_; y < max_; y += dy) {
			i = trunc((y-min_)/dy);
			for(x = min_; x < max_; x += dx) {
				j = trunc((x-min_)/dx);
				epsilon = get_epsilon(x,y);
				jane_[i][j] = epsilon;
				fprintf(fptr, "%lf\t", epsilon);
			}
			fprintf(fptr, "\n");
		}

		fclose(fptr);
		printf("Done.\n");
	}

	// H is the negative inverse of the Hessian of chi
	void CPairs::fill_H() { 

		double a, c, d, det;

		// b = c
		printf("CHECK 8\n");
		a = D2chi_Df_wn2();
		c = D2chi_Dsigma_sat_Df_wn();
		d = D2chi_Dsigma_sat2();
		det = a*d-c*c;

		H_[0][0] = -(1.0/det)*d;
		H_[0][1] = H_[1][0] = (1.0/det)*c;
		H_[1][1] = -(1.0/det)*a;
	}

	void CPairs::fill_Dchi() {

		printf("CHECK 6\n");
		Dchi_[0] = Dchi_Df_wn();
		printf("CHECK 7\n");
		Dchi_[1] = Dchi_Dsigma_sat();
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
	
	double CPairs::Dchi_Df_wn() {

		double x, y, dx = (max_-min_)/nmax_, dy = dx, Dchi_Df_wn = 0.0;
		int i, j;

		// iterate over data points in 'average.dat' and 'energy_density.dat'
		for(i = 0; i < nmax_; i++) {
			y = min_+i*dy;
			for(j = 0; j < nmax_; j++) {
				x = min_+i*dx;
				Dchi_Df_wn += -2.0*(liam_[i][j]-jane_[i][j])*Depsilon_Df_wn(x,y);
				// printf("Dchi_Df_wn = %lf\n",Dchi_Df_wn);
			}
		}

		return Dchi_Df_wn;
	}

	double CPairs::Dchi_Dsigma_sat() {

		double x, y, dx = (max_-min_)/nmax_, dy = dx, Dchi_Dsigma_sat = 0.0;
		int i, j;

		// iterate over data points in 'average.dat' and 'energy_density.dat'
		for(i = 0; i < nmax_; i++) {
			y = min_+i*dy;
			for(j = 0; j < nmax_; j++) {
				x = min_+i*dx;
				Dchi_Dsigma_sat += -2.0*(liam_[i][j]-jane_[i][j])*Depsilon_Dsigma_sat(x,y);
				// printf("Dchi_Dsigma_sat = %lf\n",Dchi_Dsigma_sat);
			}
		}

		return Dchi_Dsigma_sat;
	}

	double CPairs::D2chi_Df_wn2() { 

		double x, y, dx = (max_-min_)/nmax_, dy = dx, D2chi_Df_wn2 = 0.0;

		for(y = min_; y < max_; y += dy) {
			for(x = min_; x < max_; x += dx) {
				D2chi_Df_wn2 += 2.0*pow(Depsilon_Df_wn(x,y),2.0);
			}
		}

		return D2chi_Df_wn2;
	}

	double CPairs::D2chi_Dsigma_sat2() {

		double x, y, dx = (max_-min_)/nmax_, dy = dx, D2chi_Dsigma_sat2 = 0.0;
		int i, j;

		// iterate over data points in 'average.dat' and 'energy_density.dat'
		for(i = 0; i < nmax_; i++) {
			y = min_+i*dy;
			for(j = 0; j < nmax_; j++) {
				x = min_+i*dx;
				D2chi_Dsigma_sat2 += 2.0*(pow(Depsilon_Dsigma_sat(x,y),2.0) -
					(liam_[i][j]-jane_[i][j])*D2epsilon_Dsigma_sat_Df_wn(x,y));		
			}
		}

		return D2chi_Dsigma_sat2;
	}

	// same as D2chi_Dfwn_Dsigma_sat
	double CPairs::D2chi_Dsigma_sat_Df_wn() { 

		double x, y, dx = (max_-min_)/nmax_, dy = dx, D2chi_Dsigma_sat_Df_wn = 0.0;
		int i, j;

		// iterate over data points in 'average.dat' and 'energy_density.dat'
		for(i = 0; i < nmax_; i++) {
			y = min_+i*dy;
			for(j = 0; j < nmax_; j++) {
				x = min_+i*dx;
				D2chi_Dsigma_sat_Df_wn += 2.0*(Depsilon_Df_wn(x,y)*Depsilon_Dsigma_sat(x,y)-(liam_[i][j]-jane_[i][j])*D2epsilon_Dsigma_sat_Df_wn(x,y));		
			}
		}

		return D2chi_Dsigma_sat_Df_wn;
	}

	double CPairs::Depsilon_Df_wn(double x, double y) {

		double Depsilon_Dfwn = get_epsilon_wn(x,y)-get_epsilon_sat(x,y);

		return Depsilon_Dfwn;
	}

	double CPairs::Depsilon_Dsigma_sat(double x, double y) {

		double fwn = params_[0], Depsilon_Dsigma_sat;

		Depsilon_Dsigma_sat = fwn*Depsilon_wn_Dsigma_sat(x,y)+(1-fwn)*Depsilon_sat_Dsigma_sat(x,y);

		return Depsilon_Dsigma_sat;
	}

	double CPairs::D2epsilon_Dsigma_sat2(double x, double y) {

		double fwn = params_[0], D2epsilon_Dsigma_sat2;

		D2epsilon_Dsigma_sat2 = fwn*D2epsilon_wn_Dsigma_sat2(x,y)+(1-fwn)*D2epsilon_sat_Dsigma_sat2(x,y);

		return D2epsilon_Dsigma_sat2;
	}

	// same as D2epsilon_Df_wn_Dsigma_sat
	double CPairs::D2epsilon_Dsigma_sat_Df_wn(double x, double y) {

		double D2epsilon_Dsigma_sat_Df_wn = Depsilon_wn_Dsigma_sat(x,y)-Depsilon_sat_Dsigma_sat(x,y);

		return D2epsilon_Dsigma_sat_Df_wn;
	}

	double CPairs::Depsilon_wn_Dsigma_sat(double x, double y) {

		double T1, T2, sigma_sat = params_[1], Depsilon_wn_Dsigma_sat;

		T1 = N1_->get_T(x+0.5*b_,y);
		T2 = N2_->get_T(x-0.5*b_,y);

		Depsilon_wn_Dsigma_sat = (0.5*norm_/sigma_sat)*T1*T2*(exp(-T1*sigma_sat)+exp(-T2*sigma_sat));
		Depsilon_wn_Dsigma_sat += -(0.5*norm_*T1/(sigma_sat*sigma_sat))*(1.0-exp(-T2*sigma_sat));
		Depsilon_wn_Dsigma_sat += -(0.5*norm_*T2/(sigma_sat*sigma_sat))*(1.0-exp(-T1*sigma_sat));

		return Depsilon_wn_Dsigma_sat;
	}

	double CPairs::D2epsilon_wn_Dsigma_sat2(double x, double y) {

		double T1, T2, sigma_sat = params_[1], D2epsilon_wn_Dsigma_sat2;

		T1 = N1_->get_T(x+0.5*b_,y);
		T2 = N2_->get_T(x-0.5*b_,y);

		D2epsilon_wn_Dsigma_sat2 = (norm_*T1/pow(sigma_sat,3.0))*(1.0-exp(-T2*sigma_sat));
		D2epsilon_wn_Dsigma_sat2 += (norm_*T2/pow(sigma_sat,3.0))*(1.0-exp(-T1*sigma_sat));
		D2epsilon_wn_Dsigma_sat2 += (norm_*T1*T2/(sigma_sat*sigma_sat))*(exp(-T1*sigma_sat)+exp(-T2*sigma_sat));
		D2epsilon_wn_Dsigma_sat2 += -(0.5*norm_*T1*T2/sigma_sat)*(T1*exp(-T1*sigma_sat)+T2*exp(-T2*sigma_sat));

		return D2epsilon_wn_Dsigma_sat2;
	}

	double CPairs::Depsilon_sat_Dsigma_sat(double x, double y) {

		double T1, T2, Tmin, Tmax, sigma_sat = params_[1], Depsilon_sat_Dsigma_sat;

		// N1 centered at (-b/2,0) and N2 centered at (b/2,0)
		T1 = N1_->get_T(x+0.5*b_,y);
		T2 = N2_->get_T(x-0.5*b_,y);

		if(T1+T2==0.0) Tmin = 0.0;
		else Tmin = 2.0*T1*T2/(T1+T2);
		Tmax = 0.5*(T1+T2);

		Depsilon_sat_Dsigma_sat = (norm_/sigma_sat)*Tmin*Tmax*exp(-Tmax*sigma_sat);
		Depsilon_sat_Dsigma_sat += -(norm_/(sigma_sat*sigma_sat))*Tmin*(1.0-exp(-Tmax*sigma_sat));

		return Depsilon_sat_Dsigma_sat;
	}

	double CPairs::D2epsilon_sat_Dsigma_sat2(double x, double y) {
	
		double T1, T2, Tmin, Tmax, sigma_sat = params_[1], D2epsilon_sat_Dsigma_sat2;

		// N1 centered at (-b/2,0) and N2 centered at (b/2,0)
		T1 = N1_->get_T(x+0.5*b_,y);
		T2 = N2_->get_T(x-0.5*b_,y);

		if(T1+T2==0.0) Tmin = 0.0;
		else Tmin = 2.0*T1*T2/(T1+T2);
		Tmax = 0.5*(T1+T2);

		D2epsilon_sat_Dsigma_sat2 = (2.0*norm_*Tmin/pow(sigma_sat,3.0))*(1.0-exp(-Tmax*sigma_sat));
		D2epsilon_sat_Dsigma_sat2 += -(2.0*norm_*Tmin*Tmax/(sigma_sat*sigma_sat))*exp(-Tmax*sigma_sat);
		D2epsilon_sat_Dsigma_sat2 += -(norm_*Tmin*Tmax*Tmax/sigma_sat)*exp(-Tmax*sigma_sat);

		return D2epsilon_sat_Dsigma_sat2;
	}

} //glauber