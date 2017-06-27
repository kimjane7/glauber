#include "newton.h"
#include <vector>
#include <iostream>

template<typename T>
using matrix = std::vector<std::vector<T>>;

double* pick_estimate(double f_n, double s_n)
{

}

double* calc_energy(double f_n, double s_n)
{

}

double** jacobian(double f_n, double s_n)
{

}

matrix invert_matrix(matrix M)
{
	double M_inverse[2][2];
	//If you'd like to return a pointer to your
	//object from a function, you must use new

	double a = M[0][0];
	double b = M[0][1];
	double c = M[1][0];
	double d = M[1][1];

	std::printf("%f, %f, %f, %f",a,b,c,d);

	double determinant = a*d-b*c;

	M_inverse[0][0] = d/determinant;
	M_inverse[0][1] = a/determinant;
	M_inverse[1][0] = -b/determinant;
	M_inverse[1][1] = -c/determinant;

	return M_inverse;
}