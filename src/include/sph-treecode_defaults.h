/**@<sph-treecode_defaults.h>::**/
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <libgen.h>
#include <malloc.h>
#include <math.h>

/** pi: */
#define          pi 3.141592654

#ifdef _OPENMP
#    include <omp.h>
void            omp_set_num_threads(int numthreads);
#else
#    define omp_get_thread_num() 0
#endif

#define repeat  do
#define until(x) while(!(x))

double          max(double, double);
double          min(double, double);

#define max(x,y)    ((x)>(y)?(x):(y))
#define min(x,y)    ((x)<(y)?(x):(y))

#define     DIM     3

// #define NCOLS   1+2*DIM+6 /** change to this define as soon as individual epsilon were implemented **/
#define NCOLS   1+2*DIM+5
#define NMAX    131072 /** 2^17 particles **/
#define NCHILDREN  1<<DIM
#define MAXNN   512
#define MAXENN  1024

enum {
    MISSARGS = -666,
    OVERARGS,
    MISSFNAM,
    CANTOPEN,
    EXCEEDSMAXNN,
};

/**
seg jun  3 15:49:10 -03 2019
**/

/** useful to self-explaining index variables **/
typedef int     index_t;
typedef int     descriptor_t;
typedef int     flag_t;
typedef int     counter_t;

/***@<sph-treecode_qnn_binsort.h>::**/

typedef struct _bintree_ BINTREE;

struct _bintree_ {
    double          key;
    int             id;
    BINTREE        *left,
                   *right;
};
extern BINTREE *node;

void            qnn_binsort_put(BINTREE *);

void            qnn_binsort_get(BINTREE *);

BINTREE        *mkbnode(double key, int id);

void            qnn_binsort(size_t n);

extern double   key;

extern int      l;

extern BINTREE *root;

/***@<sph-treecode_eigenvalues.h>::**/

extern double   eigenvectors_tol;

int             power_iteration_method(double[DIM], double[DIM][DIM]);

/***@<sph-treecode_io.h>::**/

extern int      isgas[NMAX];

extern double
/**INPUT::**/
    /*mass */       m[NMAX],
    /*position vector */ x[NMAX][DIM],
    /*velocity vector */ v[NMAX][DIM],
    /*thermal specific energy */ u[NMAX],
    /*specific SPH thermal-energy rate -- use only if energy conservation equation is explicitly time integrated */
                udot[NMAX],
/**OUTPUT::**/
    /*anisotropic treecode -> gravity field */ g[NMAX][DIM],
    /*SPH density */ rho[NMAX],
                rhodot[NMAX],
    /*SPH-pressure acceleration */ P_accel[NMAX][DIM],
    /*SPH smoothing tensor eigenvalues */ eigenval[NMAX][DIM],
    /*SPH smoothing tensor eigenvectors */ eigenvec[NMAX][DIM][DIM],
    /*the smoothing lengths are the smoothing tensor eigenvalues */ h[NMAX][DIM];

extern double   u0 /* thermal specific energy of H2 molecule at lab temperature */ ;
extern double   LabTemp /* Lab temperature */ ;
extern double   vol[NMAX];

extern flag_t   deleted[NMAX];

// extern int                                      *nn[NMAX];

extern FILE    *fhandle;

/* getdimensions searches for both line and column number to define the
 * dimension limits for N and dim.
 * if there are more columns than the maximum dimension DIM the remaining
 * columns are neglect. */

size_t          getdimensions(FILE * fhandle, int *ncols, double *M);

void            show_result(size_t n /*, int ncols */ );

extern FILE    *fhandle;

/***@<sph-treecode_kernel.h>::**/

double          K3D(double);

double          DK3D(double x);

/***@<sph-treecode_sigmaoctree.h>::**/

/*
 * header file for octree definitions
 */
typedef struct _octant_ OCTANT;

extern OCTANT  *sigma_octree_root;

struct _octant_ {
    /* particles number */
    int             n;
    /* particle list */
    int            *plist;
    /* distance from the outermost particle to the center of mass */
    double          topdistance;
    /* octant's total mass */
    double          mass;
    /* octant's center of mass */
    double          x[DIM];
    /* octant's covariance matrix/eigenvectors */
    double          eigenvec[DIM][DIM];
    /* octant's covariance eigenvalues */
    double          eigenval[DIM];
    /* parental link */
    OCTANT         *parent;
    /* childhood link */
    OCTANT         *child[NCHILDREN];
};

/* procedures and functions */
void            mkoctree(OCTANT * p);

OCTANT         *mknode(int n, int *plist, OCTANT * parent);

void            mkchildhood(OCTANT * p);

int            *checkoctants(OCTANT * p);

void            getoctants(OCTANT * p, int *qcount, int *qlist[NCHILDREN]);

double          calc_center_of_mass(int n, int const *plist, double xc[]);

void            calc_covariance(int n, int *plist, double xc[], double sigma2[DIM][DIM]);

extern size_t   totalnodes;

/**Synopsis:
 * double p2p_d(double x[DIM], double y[DIM], double eigenval[DIM], double eigenvec[DIM][DIM]);
 * Computes the squared Mahalanobis (1936) distance. The covariance tensor, Sigma, is
 * represented in diagonal form, having eigenval[DIM] and eigenvec[DIM][DIM] the
 * eigenvalues and eigenvectors of Sigma, respectively. This function is defined in
 * <sph-treecode_sigmaoctree.c> */
double          p2p_d(double x[], double y[], double eigenval[], double eigenvec[DIM][DIM]);

OCTANT         *mkroot(size_t n);

/***@<sph-treecode_knn.h>::**/

typedef struct _neighboring_list_ NNLIST;
struct _neighboring_list_ {
    int             nn,
                    list[MAXNN];
    double          distance[MAXNN],
                    topdistance;
};
extern NNLIST   qnnlist;

extern int      K;

extern int      query;

extern void     aniso_k_NN_search(double, OCTANT const *);

#ifdef _RANGE_SEARCH_
extern double   R;

void            range_search(double, OCTANT *);

void            pump_neighbor(int query, int i);
#endif

/**Synopsis:
 * double p2p_distance(int i, int query);
 * Computes the squared Mahalanobis (1936) distance. The covariance tensor, Sigma, is
 * represented in diagonal form, having the global arrays
 * double eigenval[NMAX][DIM] and double eigenvec[NMAX][DIM][DIM]
 * the eigenvalues and eigenvectors of Sigma, respectively. This function is defined in <sph-treecode_knn.c> */
double          p2p_distance(int i, int query);

double          query2node_distance(double d, int query, OCTANT const *p);

void            pump_neighbor(int query, int i);

void            collect_neighbor(int query, int i);

/* This is the troublessome part of the method.
 * After the kNN it is required that the kNN-list include also the particles
 * that give non-zero contribution to the kernel, aka effective neighbors.
 * To avoid replicating kNN search and then saving time, we can reuse the
 * kNN list to add only the particles which are out from the previous kNN-list
 * but their specific moothing radius include the current query.
 */

typedef struct _ennlist_ ENNLIST;

/** The effective neighbors list: */
struct _ennlist_ {
    int             nn,
                    list[MAXENN];   // adopted MAXENN = MAXNN + sqrt(MAXX)
    double          distance[MAXENN];
};

extern ENNLIST  effneighb[NMAX];

double          p2node_outer_distance(double d, double x[], OCTANT * p);

void            add_neighbor(int query, int i);

void            effective_neighbors_search(double d, OCTANT * p);

extern size_t   leaffound,
                testednodes,
                visitednodes,
                pumped;
extern  /* symmetric_closure introduced in Aug 18, 2014, from original code knn.c at repository path: sources/SPH/sph\ 2.6/src/ */
void            symmetric_closure(size_t n);

void            sph_density(index_t i);

/** Synopsis:
 * void sph_gradients(index_t i);
 * A procedure to compute mostly SPH quantities which depend on the vector components of the smoothing-kernel gradients.
 * The procedure provides for the i-particle:
 * -> P_accel[i]: pressure acceleration;
 * -> rhodot[i]: density derivative;
 * -> udot[i]: adiabatic-heating rate.
 */
void            sph_gradients(index_t i);

double          scalarprod(double a[], double b[]);

/** binary_leapfrog(descriptor_t time_depth): performs the first half timestep at consecutive time-depth, then performs the full timestep in the root time-depth, and next performs again half timestep in the consecutive time-depth. Recently, I have adopted in-order binary integration rather than the prefix binary order of my Ph.D. Thesis (Marinho, 1997) **/
void binary_leapfrog(descriptor_t time_depth);

void            integrate(descriptor_t tdepth);

void            sph_init(size_t);

void            timebin_share(size_t n);

extern double   theta2;

extern double   epsilon2,
                epsilon;

void            init_computer_physical_scales(void);

void            find_selfregulating_kNN_cluster(size_t n);

void            smooth_velocities(size_t n);

void            sigmatreegravity(index_t particle);

double          dotproduct(double x[], double y[]);

extern double   DT;

extern double   udot_old[NMAX],
                rhodot_old[NMAX],
                dt[NMAX],
                dt_old[NMAX];

extern double   Mtot;
extern double   U,
                W,
                T,
                E;
void            integral_quantities(size_t n);

/** USAGE: \vec{g} = -m \vec{r} gHK89(|\vec{r}|) **/
double          gHK89(double);
double          fHK89(double);

extern size_t   numthreads;

extern double psize[NMAX];
