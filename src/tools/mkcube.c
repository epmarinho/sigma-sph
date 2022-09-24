#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#define DIM 3
int
main(int argc, char const *argv[])
{
    int             n = atoi(argv[1]);
    double          amplitude = .5 / n;
    double          M = 1;
    double          m = M / pow(n, DIM);
    int             i;
    srand48(time(NULL));
    if (argc > 2)
        amplitude = atof(argv[2]);
    for (i = -n / 2; i < n / 2; i++) {
        int             j;
        for (j = -n / 2; j < n / 2; j++) {
            int             k;
            for (k = -n / 2; k < n / 2; k++) {
                printf
                    ("%20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le\n",
                     m,
                     (double) i / (n - 1) + (1 -
                                             2 * drand48()) * amplitude,
                     (double) j / (n - 1) + (1 -
                                             2 * drand48()) * amplitude,
                     (double) k / (n - 1) + (1 -
                                             2 * drand48()) * amplitude, 0,
                     0, 0);
            }
        }
    }
    return 0;
}
