/**@<sph-treecode_eigenvalues.c>::**/
#include <sph-treecode_defaults.h>

/**
#include <sph-treecode_eigenvalues.h>
#include <math.h>
#include <stdlib.h>
**/

double          eigenvectors_tol = 1e-22;

int power_iteration_method(double eigenvalues[DIM], double eigenvectors[DIM][DIM])
{
    double          tensor[DIM][DIM];
    int             j,
                    iterations;

    for (j = 0; j < DIM; j++) {
        int             k;

        eigenvalues[j] = 0;
        for (k = 0; k < DIM; k++) {
            tensor[j][k] = eigenvectors[j][k];
        }
        eigenvectors[j][j] = 1;
    }
    for (j = 0; j < DIM; j++) {
        int             k,
                        converged;
        double          a = eigenvalues[j];

        iterations = 0;
        repeat {
            double          u[] = { 0, 0, 0 };
            for (k = 0; k < DIM; k++) {
                int             l;

                for (l = 0; l < DIM; l++) {
                    u[k] += tensor[k][l] * eigenvectors[j][l];
                }
            }
            for (eigenvalues[j] = k = 0; k < DIM; k++) {
                eigenvalues[j] += u[k] * u[k];
            }
            eigenvalues[j] = sqrt(eigenvalues[j]);
            for (k = 0; k < DIM; k++) {
                u[k] /= eigenvalues[j];
                eigenvectors[j][k] = u[k];
            }
            converged = (fabs(eigenvalues[j] - a) < eigenvectors_tol)
                || iterations > 127;
            a = eigenvalues[j];
            iterations++;
        }
        until(converged);

        /* tensor clean-up: dimensionality reduction */
        for (k = 0; k < DIM && j < DIM - 1; k++) {
            int             l;

            for (l = 0; l < DIM; l++) {
                tensor[k][l] -= eigenvectors[j][k] * eigenvectors[j][l] * eigenvalues[j];
            }
        }
    }
    return iterations;
}

#ifdef _EIGENVECTORS_WITH_MAIN_
#    include <stdio.h>

double          eigenvalues[DIM],
                eigenvectors[DIM][DIM];

main(int argc, char const *argv[])
{
    int             i,
                    j;

    if (argc < 10) {
        fprintf(stderr, "deu pau!\n");
        exit(-1);
    }
    for (i = 0; i < DIM; i++) {
        for (j = 0; j < DIM; j++) {
            eigenvectors[i][j] = atof(argv[i * 3 + j + 1]);
        }
    }
    if (argc > 10)
        eigenvectors_tol = atof(argv[10]);
    printf("number of iterations = %d\n", power_iteration_method(eigenvalues, eigenvectors));
    for (i = 0; i < DIM; i++) {
        printf("sigma²[%d] = %lg\ne[%d]=(", i, eigenvalues[i], i);
        for (j = 0; j < DIM - 1; j++) {
            printf("%lg, ", eigenvectors[i][j]);
        }
        printf("%lg)\n", eigenvectors[i][j]);
    }
    return 0;
}
#endif
