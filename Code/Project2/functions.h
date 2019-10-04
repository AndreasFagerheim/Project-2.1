#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <armadillo>

using namespace std;
using namespace arma;

// makes an tridiagoanl matrix (n x n)
mat tridiagonal(int n, double a, double d);

// returns an array with eiegnevalues given from analytical solutions
vec analytic_eigenvalues(int n, double a, double d);

// Find biggest matrix element in the upper half
double maxoffdiag(mat &A, int * k, int * l, int n);

// find values of sin and cos
void rotate(mat &A, mat &R, int k, int l, int n);

// jacobis method
mat jacobi_method (mat A, mat R, int n);


#endif // FUNCTIONS_H
