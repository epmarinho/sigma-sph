/**@<sph-treecode_eigenvalues.h>::**/

#pragma once

extern double   eigenvectors_tol;

extern int      power_iteration_method(double[DIM], double[DIM][DIM]);

/* define arrays to store both eigenvectors and eigenvalues of the
 * covariance tensor: */
extern double   eigenval[NMAX][DIM],
                eigenvec[NMAX][DIM][DIM]
