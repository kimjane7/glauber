#include <cmath>
#include "glauber.h"

int main(int argc, char *argv[]){
	
	CGlauber *Au = new CGlauber(197);
	CGlauber *Pb = new CGlauber(207);
	CGlauber *U = new CGlauber(235);

	printf("\n%s\n","-------- GOLD --------");
	Au->random_test(5);
	//Au->fixed_test(5);
	printf("\n%s\n","-------- LEAD --------");
	Pb->random_test(5);
	//Pb->fixed_test(5);
	printf("\n%s\n","-------- URANIUM --------");
	U->random_test(5);
	//U->fixed_test(5);


	return 0;
}