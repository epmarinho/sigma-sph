#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

main(int argc, char const *argv[])
{
    int             f = (argc == 5);
    FILE           *inbuffer = stdin;
    double          cosx = cos(atof(argv[1 + f]));
    double          cosy = cos(atof(argv[2 + f]));
    double          cosz = cos(atof(argv[3 + f]));
    double          sinx = sin(atof(argv[1 + f]));
    double          siny = sin(atof(argv[2 + f]));
    double          sinz = sin(atof(argv[3 + f]));

    double          A[3][3][3],
                    x[3],
                    xp[3];

    int             i,
                    j,
                    k;

    A[0][0][0] = 1;
    A[0][0][1] = 0;
    A[0][0][2] = 0;
    A[0][1][0] = 0;
    A[0][1][1] = cosx;
    A[0][1][2] = -sinx;
    A[0][2][0] = 0;
    A[0][2][1] = sinx;
    A[0][2][2] = cosx;

    A[1][0][0] = cosy;
    A[1][0][1] = 0;
    A[1][0][2] = siny;
    A[1][1][0] = 0;
    A[1][1][1] = 1;
    A[1][1][2] = 0;
    A[1][2][0] = -siny;
    A[1][2][1] = 0;
    A[1][2][2] = cosy;

    A[2][0][0] = cosz;
    A[2][0][1] = -sinz;
    A[2][0][2] = 0;
    A[2][1][0] = sinz;
    A[2][1][1] = cosz;
    A[2][1][2] = 0;
    A[2][2][0] = 0;
    A[2][2][1] = 0;
    A[2][2][2] = 1;

    if (f)
        inbuffer = fopen(argv[1], "r");

    while (1) {
        double          m;

        int             c =
            fscanf(inbuffer, "%lg %lg %lg %lg", &m, x + 0, x + 1, x + 2);

        if (c == -1)
            break;

        for (k = 0; k < 3; k++) {
            for (i = 0; i < 3; i++) {
                for (xp[i] = j = 0; j < 3; j++) {
                    xp[i] += A[k][i][j] * x[j];
                }
            }
            memcpy(x, xp, 3 * sizeof(double));
        }

        fprintf(stdout, "%lg %lg %lg %lg\n", m, x[0], x[1], x[2]);

        while ((c = getc(inbuffer)) != '\n' && c != -1);

    }

    return 0;
}
