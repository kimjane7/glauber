#ifndef GLAUBER_H
#define GLAUBER_H
#include <cstdio>
#include <boost/random.hpp>

const double pi = 4.0*atan(1.0);
// Uniform rand dist over interval, takes RNG as arg
using boost::uniform_01;
// Mersenne Twister RNG
typedef boost::mt19937 RNGType;

class CGlauber {
public:
	// Initialize member variables
	const static double a_=0.53, rho0_=0.14;
	double R_;
	int A_;

	// Declare member functions
	CGlauber(int Aset);
	double get_rho(double x, double y, double z);
	double get_T(double x, double y);
	double random_01(RNGType generator);
	void random_test(int nmax);
};

CGlauber::CGlauber(int A_set) {
	// Set the nucleon number
	A_ = A_set;

	double A_est, A_prime_est;
	
	// Estimate radius from nucleon number
	double R_est = pow(3.0*A_/(4.0*pi*rho0_), 1.0/3.0);
	double R_diff; // Newton's method increment

	double r;
	double dr = 1E-3*R_est;

	// Perform Newton's Method to find actual radius
	do {
		A_est = 0.0;

		// Calculate A from WS dist, density, and radius guess
		for (r = 0.0; r < R_est; r += dr) {
			A_est += 4.0*pi*rho0_*r*r*dr/(1.0+exp((r-R_est)/a_));
		}
		
		// Update radius guess
		A_prime_est = 2.0*pi*rho0_*R_est*R_est;
		R_diff = (A_-A_est)/A_prime_est;
		R_est = R_est+R_diff;
	
	} while(R_diff > dr);

	R_ = R_est;
}

double CGlauber::get_rho(double x, double y, double z) {
	// Calculate distance from center
	double r = pow(x*x+y*y+z*z,0.5);
	
	// Calculate Woods-Saxon density
	double rho = rho0_/(1.0+exp((r-R_)/a_));
	
	return rho;
}

double CGlauber::get_T(double x, double y) {
	// Initialize variables
	double T=0.0, z, z_bound = pow(R_*R_-x*x-y*y,0.5), dz = 2.0E-3*z_bound;
	
	// Integrate the density (rho) in the z direction
	// to get thickness function
	for(z=-z_bound; z<z_bound; z+=dz){
		T+=get_rho(x,y,z)*dz;
	}
	
	return T;
}

double CGlauber::random_01(RNGType generator) {
	static uniform_01<RNGType> dist(generator);
	return dist();
}

void CGlauber::random_test(int nmax) {
	RNGType generator(time(0));
	printf("A = %i \nR = %lf fm\n", A_, R_);
	double r, theta, phi, x, y, z;
	for(int n = 0; n < nmax; n+=1){
		// Pick a random point in the nuclear space
		r = R_*random_01(generator);
		theta = pi*random_01(generator);
		phi = 2*pi*random_01(generator);
		
		// Transform coords
		x = r*sin(theta)*cos(phi);
		y = r*sin(theta)*sin(phi);
		z = r*cos(theta);
		
		// Print test values
		printf("rho(%.3lf, %.3lf, %.3lf)=%lf \t T(%.3lf, %.3lf)=%lf\n",
			   x,y,z,get_rho(x,y,z),x,y,get_T(x,y));
	}
}

#endif // GLAUBER_H