#include "glauber.h"
#include <iostream>

using namespace glauber;

int main(int argc, char *argv[]){
	
	// make new nuclei
	Nucleus* Pb1 = new Nucleus(207);
	Nucleus* Pb2 = new Nucleus(207);

	// make Pb-Pb nuclei pair with 1.0 fm impact parameter
	NucleusPair* Pb_Pb = new NucleusPair(Pb1, Pb2, 1.0, -10.0, 10.0, 100);

	Pb_Pb->minimize_chi2();

	Pb_Pb->test_chi2();

	return 0;
}