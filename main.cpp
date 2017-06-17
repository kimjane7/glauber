#include "glauber.h"

using namespace glauber;

int main(int argc, char *argv[]){
	
	// make new nuclei
	CNucleus *Pb1 = new CNucleus(207);
	CNucleus *Pb2 = new CNucleus(207);

	// make Pb-Pb nuclei pair with 10.0 fm impact parameter
	CPairs *Pb_Pb = new CPairs(Pb1, Pb2, 10.0, 40.0, 0.5);
	Pb_Pb->print_eps(-5.0, 5.0, 100);

	return 0;
}