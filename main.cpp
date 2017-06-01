#include "glauber.h"

using namespace std;

int main(int argc, char *argv[]){

	// Check for correct number of arguments
	if (argc != 4){
        printf("Incorrect number of arguments, try again.\n");
        exit(-1);
	}

	// Energy per rapidity, varies from 0.85-1.2
	double dEdy = atof(argv[1]);
	// Saturation cross section, varies from 30-50
	double sigma_sat = 100.0*atof(argv[2]);
	// Wounded nucleon vs Saturation model weight
	double f_wn = atof(argv[3]); // Varies from 0-1
	
	CNucleus *Au = new CNucleus(197);
	CNucleus *Pb = new CNucleus(207);
	CNucleus *U = new CNucleus(235);

	Au->fixed_test(5,"testing.dat");

	/* WRITE FIXED TEST AND WRITE EPS TO FILE
	printf("\n%s\n","-------- GOLD --------");
	Au->random_test(5);

	printf("\n%s\n","-------- LEAD --------");
	Pb->random_test(5);

	printf("\n%s\n","-------- URANIUM --------");
	U->random_test(5);
	*/

	CPairs *Au_Pb = new CPairs(Au, Pb, 1.0, dEdy, sigma_sat, f_wn);
	printf("Energy density = %lf\n", Au_Pb->get_epsilon(0.5,0.5));



	return 0;
}