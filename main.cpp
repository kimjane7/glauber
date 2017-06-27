#include "glauber.h"
#include "newton.h"
#include <iostream>

using namespace glauber;

int main(int argc, char *argv[]){

	// make new nuclei
	CNucleus *Pb1 = new CNucleus(207);
	CNucleus *Pb2 = new CNucleus(207);

	// make Pb-Pb nuclei pair with 10.0 fm impact parameter
	CPairs *Pb_Pb = new CPairs(Pb1, Pb2, 10.0, -10.0, 10.0, 100);
	Pb_Pb->minimize_chi();


/*
	double a[2][2] = {{1,2},{3,4}};
	double result[2][2];

	result = invert_matrix(a);

	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			std::cout << result[i][j] << std::endl;
		}
	}
	
>>>>>>> f57d3fd3709d8d82371cb4ab4345de1f31eeb555
*/

	return 0;
}