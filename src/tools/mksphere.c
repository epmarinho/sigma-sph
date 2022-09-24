/*
 * mksphere is an accept-reject Monte-Carlo method to generate homogeneous
 * spherical particle distribution
 */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

void
usage(const char *cmdname)
{
    char           *cmdformat =
        "%s n R0 R M V Xoff Yoff Zoff VXoff VYoff VZoff u\n";
    printf(cmdformat, cmdname);
}

unsigned short  Xi[3] = { 0x1FF1, 0xAE0D, 0x3D2F };

int main(int argc, const char *argv[])
{
    if (argc < 13) {
        printf("usage: ");
        usage(argv[0]);
        exit(-1);
    } else {
        int             n = atoi(argv[1]);
        double          R0 = atof(argv[2]);
        double          R = atof(argv[3]);
	double		ratio = pow(R0 / R, 2);
        double          M = atof(argv[4]);
        double          m = M / n;
        double          u;
        double          V = atof(argv[5]);
        double          offset[3];
        double          voffset[3];
        offset[0] = atof(argv[6]);
        offset[1] = atof(argv[7]);
        offset[2] = atof(argv[8]);
        voffset[0] = atof(argv[9]);
        voffset[1] = atof(argv[10]);
        voffset[2] = atof(argv[11]);
        u = atof(argv[12]);
        srand48(time((time_t *) Xi));
        srand48(time((time_t *) & Xi[1]));
        int             i;
        for (i = 0; i < n; i++) {
            int             j;
            double          r,
                            s,
                            x[3],
                            v[3];
            do {
                for (r = j = 0; j < 3; j++) {
                    x[j] = (2 * erand48(Xi) - 1);
                    r += x[j] * x[j];
                }
            } while (ratio > r || r > 1);
            do {
                for (s = j = 0; j < 3; j++) {
                    v[j] = (2 * erand48(Xi) - 1);
                    s += v[j] * v[j];
                }
            } while (s > 1);
            for (j = 0; j < 3; j++) {
                x[j] *= R;
                v[j] *= V * R * R / (r * 128 + R * R);
            }
            /*
             *          for (j = 0; j < 3; j++) {
             *          v[j] -= .25 * x[j];
             }
             */


            printf
                ("%20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le\n",
                 m, x[0] + offset[0], x[1] + offset[1], x[2] + offset[2],
                 v[0] + voffset[0], v[1] + voffset[1], v[2] + voffset[2],
                 u);
        }
        return 0;
    }
}
