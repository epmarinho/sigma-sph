/**@<sph-treecode_kernel.c>::**/
/*****************************************************************************
 * Smoothing kernel is modeled from the kernel based density estimation      *
 * theory with a compact support set by k-NN method.                         *
 * The present model is the quartic-spline, whose coefficients were computed *
 * using the symbolic processing package Maple 15 (2010)                     *
 *****************************************************************************/

#include <sph-treecode_defaults.h>

/** See Monaghan (1992), page 554, for a 3D cubic B-spline **/
double K3D(double x /** WARNING: 0 <= x <= 2 */ )
{
    double          K;

    x = fabs(x);

    if (x > 2)
        K = 0;
    else if (x > 1)
        /* 1 <= x <= 2 */
        K = .25 * pow(2 - x, 3.);
    else
        /* 0 <= x <= 1 */
        K = 1 + (.75 * x - 1.5) * pow(x, 2.);
    return 0.318309886 * K;
}

double DK3D(double x /** WARNING: 0 <= x <= 2 */ )
{
    double          K;

    x = fabs(x);

    if (x > 2)
        /* if x > 2 */
        K = 0;
    else if (x > 1)
        /* 1 <= x <= 2 */
        K = -.75 * pow(2 - x, 2.);
    else
        /* 0 <= x <= 1 */
        K = (-3 + 2.25 * x) * x;
    return 0.318309886 * K;
}

// #define _HK89_SOFTENING_
/** I have no idea on how to implement quadrupole correction for HK89 softening once the splined kernel is of analytical class C²,
 *  since the kernel's third derivative is discontinuous. Hence, the HK89 idea have to be replaced with a quartic B-spline. **/
#ifdef _HK89_SOFTENING_

/** USAGE: phi = m fHK89(r) **/
double fHK89(double r)
{
    //     if (r < 0) {
    //         fprintf(stderr, "double fHK89(double r): Invalid HK89 spline argument. Must be non-negative. Exiting with error status");
    //         fprintf(stderr, " -1024\n");
    //         exit(-1024);
    //     }

    double          u = r / epsilon;

    if (u <= 1) {
        return -2 / epsilon * (1. / 3 * pow(u, 2) - .15 * pow(u, 4) + .05 * pow(u, 5)) + 1.4 / epsilon;
    }
    if (u <= 2) {
        return -1. / (15 * r) - 1 / epsilon * (4. / 3 * pow(u, 2) - pow(u, 3) + .3 * pow(u, 4) - 1. / 30 * pow(u, 5)) +
            1.6 / epsilon;
    }
    return 1. / r;
}

/** USAGE: \vec{g} = -m \vec{r} gHK89(|\vec{r}|) **/
double gHK89(double r)
{
    //     if (r < 0) {
    //         fprintf(stderr, "double gHK89(double r): Invalid HK89 spline argument. Must be non-negative. Exiting with error status");
    //         fprintf(stderr, " -1024\n");
    //         exit(-1024);
    //     }

    double          u = r / epsilon;

    if (u <= 1) {
        return (4. / 3 - 1.2 * pow(u, 2) + .5 * pow(u, 3)) / pow(epsilon, 3);
    }
    if (u <= 2) {
        return (-1. / 15 + 8. / 3 * pow(u, 3) - 3 * pow(u, 4) + 1.2 * pow(u, 5) - 1. / 6 * pow(u, 6)) / pow(r, 3);
    }
    return 1. / pow(r, 3);
}
#endif
