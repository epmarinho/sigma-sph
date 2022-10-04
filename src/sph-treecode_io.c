/**@<sph-treecode_io.c>::**/
#include <sph-treecode_defaults.h>

double
/*mass */       m[NMAX],
/*position vector */ x[NMAX][DIM],
/*velocity vector */ v[NMAX][DIM];

double
/*anisotropic treecode -> gravity field */ g[NMAX][DIM],
/*SPH-pressure acceleration */ P_accel[NMAX][DIM],
/*Covariance's eigenvectors */ eigenval[NMAX][DIM],
/*and eigenvalues */ eigenvec[NMAX][DIM][DIM];

double
/*SPH density */ rho[NMAX];

double
/*SPH density rate */ rhodot[NMAX];

double
/*thermal specific energy */ u[NMAX];

double
/*thermal energy rate */ udot[NMAX];

int
/**isgas is a flag to tell the simulation the particle is an SPH particle but not a single n-body particle */ isgas[NMAX];

FILE           *fhandle;

/* getdimensions searches for both line and column number to define the
 * dimension limits for N and dim.
 * if there are more columns than the maximum dimension DIM the remaining
 * columns are neglect. */
size_t getdimensions(FILE * fhandle, int *ncols, double *Mtot)
{
    /** Last revision on Dec 3, 2020 **/
    size_t          i = 0;  /* i denotes line counter */
    int             pos = 0;
    int             c;
    int             dimensions = NCOLS;

    *Mtot = 0;

  blank_or_hashmark:
    while (isspace(c = getc(fhandle)));
    if (c == '#') {
        while ((c = getc(fhandle)) != '\n');
        goto blank_or_hashmark;
    }
    ungetc(c, fhandle);

    while (1) {

        if (pos == dimensions) {
            while ((c = getc(fhandle)) != '\n');
            *ncols = pos;
            pos = 0;
            i++;
        }
        while (isspace(c = getc(fhandle)));
        if (c == EOF) {
            return min(i, NMAX);
        }
        ungetc(c, fhandle);
        if (pos < dimensions) {
            int             nr = -1;

            if (pos == 0) {
                nr = fscanf(fhandle, "%lg", &m[i]); // 1
                *Mtot += m[i];
            } else if (pos <= DIM) {
                nr = fscanf(fhandle, "%lg", &x[i][pos - 1]);    // 2 - 4
            } else if (pos <= 2 * DIM) {
                nr = fscanf(fhandle, "%lg", &v[i][pos - DIM - 1]);  // 5 - 7
            } else if (pos == 2 * DIM + 1) {
                nr = fscanf(fhandle, "%lg", &u[i]); // 8
            } else if (pos == 2 * DIM + 2) {
                nr = fscanf(fhandle, "%lg", &dt_old[i]);    // 9
            } else if (pos == 2 * DIM + 3) {
                nr = fscanf(fhandle, "%lg", &udot[i]);  // 10
            } else if (pos == 2 * DIM + 4) {
                nr = fscanf(fhandle, "%lg", &rhodot[i]);    // 11
            } else if (pos == 2 * DIM + 5) {
                nr = fscanf(fhandle, "%d", &isgas[i]);  // 12
            }
            if (nr == 0) {
                fprintf(stderr, "sph-treecode_io: warning: no characters read\n");
            }
            pos++;
            if (i == 0) {
                if (pos < NCOLS) {
                    while (isspace(c = getc(fhandle))) {
                        if (c == '\n') {
                            *ncols = dimensions = pos;
                            pos = 0;
                            i++;
                        }
                    }
                    ungetc(c, fhandle);
                } else if (pos == NCOLS) {
                    while ((c = getc(fhandle)) != '\n');
                    *ncols = dimensions = pos;
                    pos = 0;
                    i++;
                }
            }
        }
    }

    return 0;
}

/** specific thermal energy first use in show_result() **/
double          u[NMAX];

//flag_t          deleted[NMAX];
extern double   oneoverCv[NMAX];

void show_result(size_t n /*, int ncols */ )
{
    int             i,
                    l;
    char           *format1 = "%20.16le\t";

    extern double   u_unit /*, rho_unit */ ;
    extern double   meanweight,
                    Rgas;

    /*
     *        char *format2 =
     *        "%20.16le\t%20.16le\t%20.16le\t%20.16le\t%20.16le\t%20.16le\t%20.16le";
     *        char *format3 = "%20.16le\n";
     */

    for (i = 0; i < n; i++) {
        if (deleted[i])
            continue;
        printf(format1, m[i]);  // column 1
        for (l = 0; l < DIM; l++) {
            printf(format1, x[i][l]);   // column 2 - 4
        }
        for (l = 0; l < DIM; l++) {
            printf(format1, v[i][l]);   // column 5 - 7
        }
        printf(format1, u[i]);  // column 8
        printf(format1, dt_old[i]); // column 9
        printf(format1, udot[i]);   // column 10
        printf(format1, rhodot[i]); // column 11
        printf("\t%d\t", isgas[i]); // column 12
        printf(format1, 1 + oneoverCv[i]);  // column 13
        printf(format1, rho[i]);    // * rho_unit); // column 14
        printf(format1, (meanweight / Rgas * oneoverCv[i] * u[i] * u_unit) /* ⁰K */ );    // column 15
        for (l = 0; l < DIM; l++) {
            printf(format1, g[i][l]);
        }   // columns 16 - 18
        //         for (l = 0; l < DIM; l++) {
        //             printf(format1, eigenval[i][l]);
        //         }   // columns 16 - 18
        /*WARNING: don't touch next line ever!!! */
        printf("\n");
    }
}
