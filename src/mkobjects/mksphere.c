/*
 * mksphere is an accept-reject Monte-Carlo method to generate homogeneous spherical particle distribution
 */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
/*
 * Astronomical parameters and constants -- check they out again and again!
 */
double          MSun = 1.9891e+33 /* the sun mass in grams */ ;
double          AU = 1.49597871e+13 /* astronomical unity in cm */ ;
double          Oort = 1e5 * 1.49597871e+13 /* inner radius of the Oort region cm */ ;
double          pc = 3.08567758e+18 /* parsec in centimeters (see IAU definition) */ ;
double          ly = 3.08567758e+18 / 3.2616 /* light-year in cm -- a parser corresponds to 3.2616 light-years */ ;
double          day = 3600 * 24 /* 1 day in seconds */ ;
double          year = 365.25 * 3600 * 24 /* 1 year in seconds */ ;
double          Myr = 1e6 * 365.25 * 3600 * 24 /* 1 milion years in seconds */ ;
double          MGalaxy = 1.5e+12 * 1.9891e+33 /* estimated mass of our galaxy, the Via Lactea, in grams */ ;
double          RGalaxy = 15 * 1000 * 3.08567758e+18 /* galactic radius in cm, see https://en.wikipedia.org/wiki/Milky_Way */ ;
double          zGalaxy = 0.3 * 3.08567758e+18 * 1000   /* estimated Galaxy disc thickness in cetimeters, see
* https://en.wikipedia.org/wiki/Milky_Way */ ;
double          TSun = 240 * 1e6 * 365.25 * 23.9344699 * 3600   /* estimated period of transit of the sun around the Galaxy center
* in seconds */ ;
double          RSun = 8.34 * 3.08567758e+18 * 1000 /* this is the R0 in cm, see https://en.wikipedia.org/wiki/Milky_Way */ ;
/*
 * Fundamental physical constants
 */
double          c = 29979245800 /* the speed of the light in cm s⁻¹ */ ;
double          G = 6.67428e-8 /* the universal gravitation constant in cm³g⁻¹s⁻² */ ;
/*
 * Simulation adimensional unities
 */
double          lunit,
tunit,
munit;
/*
 * normalized speed of the light:
 */
double          c_norm;
void
init_computer_physical_scales(void)
{
    lunit = 39.48 * AU /* Pluto's semi-major axis */ ;
    munit = 1.0014 /* == present-day solar system total mass */  * MSun;
    /*
     * gravity acceleration is given by
     * G M / r² in lunity * tunity⁻²
     * so that
     * since G = 6.67428e-8 cm³g⁻¹s⁻² == 1 lunit³ munit⁻¹ tunit⁻²
     * thus,
     */
    tunit = sqrt(pow(lunit, 3) / (G * munit));
    /*
     * the normalized speed of the light shall be useful in a near future:
     */
    c_norm = c * tunit / lunit;
    printf("# computer-physical units are given by\n");
    printf("# length = %lg cm = %lg au\n", lunit, lunit / AU);
    printf("# time = %lg s = %lg yr\n", tunit, tunit * 1e6 / Myr);
    printf("# mass = %lg g = %lg solar masses\n", munit, munit / MSun);
    printf("# derivative quantities:\n");
    printf("# velocity = %lg cm s⁻¹ = %lg km h⁻¹\n", lunit / tunit, lunit / tunit * 3600 / 1e5);
    double          omegaunit = 2 * 4 * atan(1.0) / tunit;
    printf("# angular velocity = %lg radians s⁻¹ = %lg revolutions per century\n", omegaunit, 100 * year / tunit);
    printf("# acceleration = %lg cm s⁻²\n", lunit / pow(tunit, 2));
    printf("# energy = %lg erg\n", munit * pow(lunit / tunit, 2));
    printf("#\n");
}

unsigned short  Xi[3] = { 0x1FF1, 0xAE0D, 0x3D2F };

void            srand48(long int seedval);
double          erand48(unsigned short xsubi[3]);

void
usage(const char *cmdname)
{
    char           *cmdformat = "%s n R0 R M Vdisp Omega Xoff Yoff Zoff VXoff VYoff VZoff u is_gas epsilon\n";
    printf(cmdformat, cmdname);
    exit(EXIT_FAILURE);
}

int
main(int argc, char const *argv[])
{
    if (argc < 16) {
        printf("%s missing parameters\nusage: ", argv[0]);
        usage(argv[0]);
        exit(EXIT_FAILURE);
    } else {
        int             fdrnd = open("/dev/random", O_RDONLY);
        int             is_gas;
        int             n = atoi(argv[1]);
        double          R0 = atof(argv[2]);
        double          R = atof(argv[3]);
        double          M = atof(argv[4]);
        double          m = M / n;
        double          u;
        double          V = atof(argv[5])/**V is the isotropic, Gaussian velocity dispersion*/;
        double          Omega = atof(argv[6])/**Omega is the z-aligned angular momentum*/;
        double          offset[3];
        double          voffset[3];
        double          epsilon;
        double          r,
        s,
        *x[3],
        *v[3];
        
        offset[0] = atof(argv[7]);
        offset[1] = atof(argv[8]);
        offset[2] = atof(argv[9]);
        voffset[0] = atof(argv[10]);
        voffset[1] = atof(argv[11]);
        voffset[2] = atof(argv[12]);
        u = atof(argv[13]);
        is_gas = atoi(argv[14]);
        epsilon = atof(argv[15]);
        
        if (fdrnd) {
            read(fdrnd, Xi, 6);
        }
        srand48(time((time_t *) Xi));
        srand48(time((time_t *) & Xi[1]));
        int             i;
        x[0] = (double *) malloc(n * sizeof(double));
        x[1] = (double *) malloc(n * sizeof(double));
        x[2] = (double *) malloc(n * sizeof(double));
        v[0] = (double *) malloc(n * sizeof(double));
        v[1] = (double *) malloc(n * sizeof(double));
        v[2] = (double *) malloc(n * sizeof(double));
        
        /*
         * simulate a homogeneous distribution function
         */
        for (i = 0; i < n; i++) {
            int             j;
            
            do {
                for (r = j = 0; j < 3; j++) {
                    x[j][i] = (2 * erand48(Xi) - 1);
                    r += x[j][i] * x[j][i];
                }
            } while (R0 / R > r || r > 1);
            
            do {
                for (s = j = 0; j < 3; j++) {
                    v[j][i] = (2 * erand48(Xi) - 1);
                    s += v[j][i] * v[j][i];
                }
            } while (s > 1);
            
            for (j = 0; j < 3; j++) {
                x[j][i] *= R;
                v[j][i] *= V;
            }
            
            v[1][i] += Omega * x[0][i];
            v[0][i] += -Omega * x[1][i];
        }
        /*
         * take the center of mass and set it as the origin of the reference frame
         */
        double          xc[3],
        vc[3];
        int             j;
        
        for (j = 0; j < 3; j++) {
            for (i = 0; i < n; i++) {
                xc[j] += x[j][i];
                vc[j] += v[j][i];
            }
            xc[j] /= n;
            vc[j] /= n;
            for (i = 0; i < n; i++) {
                x[j][i] -= xc[j];
                v[j][i] -= vc[j];
            }
        }
        /*
         * show result
         */
        for (i = 0; i < n; i++) {
            #ifdef     __IAU_UNITS__
            
            init_computer_physical_scales();
            double          uunit = pow(lunit / tunit, 2) * munit,
            vunit = lunit / tunit;
            printf
            ("%20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le\n",
             m * munit,
             (x[0][i] + offset[0]) * lunit,
             (x[1][i] + offset[1]) * lunit,
             (x[2][i] + offset[2]) * lunit,
             (v[0][i] + voffset[0]) * vunit, (v[1][i] + voffset[1]) * vunit, (v[2][i] + voffset[2]) * vunit, 0, u * uunit, 0., 0.);
            
            #else
            
            printf
            /**   m        x        y        z        vx       vy       vz       u        dt       udot   rhodot is_gas epsilon*/
            ("%20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %d %20.16le\n",
            /**   1        2        3        4        5        6        7        8        9        10       11   12    13 **/
             m,
             x[0][i] + offset[0], x[1][i] + offset[1], x[2][i] + offset[2],
             v[0][i] + voffset[0], v[1][i] + voffset[1], v[2][i] + voffset[2],
             u,
             0., 0., 0., is_gas, epsilon);
            
            #endif
        }
    }
    return 0;
}
