#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int
main(int argc, const char *argv[])
{
    FILE *inputbuffer = stdin;
    if (argc > 1) {
	    inputbuffer = fopen(argv[1], "r");
	    if (inputbuffer == NULL) {
		    fprintf(stderr, "%s: cannot open %s or file not found\n", argv[0], argv[1]);
		    exit(EOF);
	    }
    }
    double          E,
		    E0,
                    Emax,
                    Emin,
		    Emean = 0,
                    E2 = 0,
		    sigma;
    int
                    n = 0;
    int             first = 1;
    while (1) {
        fscanf(inputbuffer, "%lg", &E);
        if (feof(inputbuffer))
            break;
        if (first) {
            Emax = Emin = E0 = E;
            first = 0;
        } else {
            Emax = E > Emax ? E : Emax;
            Emin = E < Emin ? E : Emin;
        }
	Emean += E;
        E2 += E * E;
        n++;
    }
    Emean /= n;
    E2 /= n;
    fprintf(stdout, "mean energy = %lg\n", Emean);
    fprintf(stdout, "RMS energy = %lg\n", sqrt(E2));
    fprintf(stdout, "standard deviation = %lg\n", sigma = sqrt(E2 - Emean * Emean));
    fprintf(stdout, "conservation error = %lg%%\n", sigma / sqrt(E2) * 100);
    fprintf(stdout, "energy has been %lg%% conserved\n", (1 - sigma) * 100);
    fprintf(stdout, "min to max relative amplitude to the RMS = %lg%%\n",
            (Emax - Emin) / sqrt(E2) * 100);
    if (E0 < E) {
    	fprintf(stdout, "energy gain (false heating) = %lg%%\n", (1 - E / E0) * 100);
    } else {
    	fprintf(stdout, "energy loss (false cooling) = %lg%%\n", (1 - E0 / E) * 100);
    }
    return 0;
}
