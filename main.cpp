#include "glauber.h"

using namespace std;

int main(int argc, char *argv[]){

	if (argc != 4){
        printf("Try again. \n");
        exit(-1);
	}

	double dEdy = atof(argv[1]);
	double sig_sat = 100.0*atof(argv[2]);
	double fwn = atof(argv[3]);
	
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

	CPairs *Au_Pb = new CPairs(Au, Pb, 1.0, dEdy, sig_sat, fwn);
	printf("\neps = %lf\n",Au_Pb->get_eps(0.5,0.5));



	return 0;
}