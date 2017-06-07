#include "glauber.h"


int main(int argc, char *argv[]){

	// check for correct number of arguments
	if (argc != 4){
        printf("Input dE/dy, saturation cross section, and weight fwn.\n");
        exit(-1);
	}

	// energy per rapidity, varies from 0.85-1.2
	double dEdy = atof(argv[1]);

	// saturation cross section, varies from 30-50 mb
	double sigma_sat = 100.0*atof(argv[2]);

	// wounded nucleon vs saturation model weight
	double f_wn = atof(argv[3]); // Varies from 0-1
	
	// make new nuclei
	CNucleus *Au = new CNucleus(197);
	CNucleus *Pb = new CNucleus(207);
	CNucleus *U = new CNucleus(235);

	// print radii
	printf("Radius of Au = %lf\n", Au->get_R());
	printf("Radius of Pb = %lf\n", Pb->get_R());
	printf("Radius of U = %lf\n", U->get_R());

	return 0;
}