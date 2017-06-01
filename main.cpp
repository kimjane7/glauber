#include <cmath>
#include "glauber.h"

using namespace std;

int main(int argc, char *argv[]){
	
	CNucleus *Au = new CNucleus(197);
	CNucleus *Pb = new CNucleus(207);
	CNucleus *U = new CNucleus(235);

	/* WRITE FIXED TEST
	printf("\n%s\n","-------- GOLD --------");
	Au->random_test(5);

	printf("\n%s\n","-------- LEAD --------");
	Pb->random_test(5);

	printf("\n%s\n","-------- URANIUM --------");
	U->random_test(5);
	*/

	CPairs *Au_Pb = new CPairs(Au, Pb, 1.0);
	printf("\neps = %lf",Au_Pb->get_eps(0.5,0.5));



	return 0;
}