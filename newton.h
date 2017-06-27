#ifndef NEWTON_H
#define NEWTON_H

#include <cstddef> // for size_t
#include <algorithm> // for transform

double* pick_estimate(double f_n, double s_n);
double* calc_energy(double f_n, double s_n);
double** jacobian(double f_n, double s_n);
double** invert_matrix(double M[2][2], double J_inverse[2][2]);

template <class Type, std::size_t size>
void add_arrays(const Type(&a)[size],
	const Type(&b)[size], Type(&result)[size])
{
	std::transform(a, a+size, b, result, std::plus<Type>());
}

#endif // NEWTON_H