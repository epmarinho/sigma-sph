/**@<sph-treecode_sigmaoctree.c>::**/
#include <sph-treecode_defaults.h>

/*** WARNING: omp parallelism did not work in this module! ***/

double psize[NMAX];

OCTANT         *sigma_octree_root;

OCTANT         *mkroot(size_t n)
{
    /** create the covariance octree root **/
    index_t        *plist = malloc(n * sizeof(int)),
        i;

    for (i = 0; i < n; i++)
        plist[i] = i;

    return mknode(n, plist, NULL);
}

void mkoctree(OCTANT * p)
{
    if (p->n > 1) {
        int             q;

        mkchildhood(p);
        for (q = 0; q < NCHILDREN; q++) {
            OCTANT         *child = p->child[q];

            if (child) {
                mkoctree(child);
            }
        }
    }
}

void mkchildhood(OCTANT * p)
{
    int            *qlist[NCHILDREN];
    int            *qcount = checkoctants(p);
    int             q;

    for (q = 0; q < NCHILDREN; q++) {
        if (qcount[q])
            qlist[q] = malloc(qcount[q] * sizeof(int));
    }
    getoctants(p, qcount, qlist);
    for (q = 0; q < NCHILDREN; q++) {
        if (qcount[q])
            p->child[q] = mknode(qcount[q], qlist[q], p);
    }
}

#include <strings.h>    //required by explicit_bzero
void            explicit_bzero(void *s, size_t n);

int            *checkoctants(OCTANT * p)
{
    int             i;
    int            *qcount = malloc(sizeof(int) * NCHILDREN);

    explicit_bzero(qcount, sizeof(int) * NCHILDREN);
    for (i = 0; i < p->n; i++) {
        int             query = p->plist[i],
            q;
        int             k;

        for (q = k = 0; k < DIM; k++) {
            double          a;
            int             j;

            for (a = j = 0; j < DIM; j++) {
                a += (x[query][j] - p->x[j]) * p->eigenvec[k][j];
            }
            q = (q << 1) + (a > 0);
        }
        qcount[q]++;
    }
    return qcount;
}

void getoctants(OCTANT * p, int *qcount, int *qlist[NCHILDREN])
{
    int             i;

    explicit_bzero(qcount, sizeof(int) * NCHILDREN);

    for (i = 0; i < p->n; i++) {
        int             query = p->plist[i],
            q = 0,
            k;

        for (k = 0; k < DIM; k++) {
            double          a;
            int             j;

            for (a = j = 0; j < DIM; j++) {
                a += (x[query][j] - p->x[j]) * p->eigenvec[k][j];
            }
            q = (q << 1) + (a > 0);
        }
        qlist[q][qcount[q]++] = query;
    }

}

double          maxfromlist(OCTANT * p);

OCTANT         *mknode(int n, int *plist, OCTANT * parent)
{
    int             i /*, iterations */ ;
    OCTANT         *p = malloc(sizeof(OCTANT));

    if (n == 1) {
        psize[plist[0]] = parent->topdistance;
    }
    p->n = n;
    p->plist = plist;
    p->mass = calc_center_of_mass(n, p->plist, p->x);
    calc_covariance(n, p->plist, p->x, p->eigenvec);
    //     iterations =
    power_iteration_method(p->eigenval, p->eigenvec);
    p->topdistance = maxfromlist(p);
    for (i = 0; i < NCHILDREN; i++) {
        p->child[i] = NULL;
    }
    p->parent = parent;
    return p;
}

/*NOTE: compute point-to-point Mahalanobis distance under eigenval[] and eigenvec[] */
double p2p_d(double x[], double y[], double eigenval[], double eigenvec[DIM][DIM])
{
    double          xi2;
    int             k;

    for (xi2 = k = 0; k < DIM; k++) {
        int             j;
        double          a;

        for (a = j = 0; j < DIM; j++) {
            a += (x[j] - y[j]) * eigenvec[k][j];
        }
        xi2 += a * a / eigenval[k];
    }
    /* return the point-to-point, squared Mahalanobis distance: */
    return xi2;
}

double maxfromlist(OCTANT * p)
{
  /** this method performs a rough estimation of the octant size **/
    int             i;
    double          topdistance = 0, distance;

    for (i = 0; i < p->n; i++) {
        // topdistance = max(topdistance, p2p_d(p->x, x[p->plist[i]], p->eigenval, p->eigenvec));
        distance = 0;
        for (int j = 0; j < DIM; j++) {
            distance += pow(p->x[j] - x[p->plist[i]][j], 2);
        }
        distance = sqrt(distance);
        topdistance = max(topdistance, distance);
    }
    return topdistance;
}

double calc_center_of_mass(int n, int const *plist, double xc[])
{
    double          totmass = 0;

    if (n == 1) {
        memcpy(xc, x[plist[0]], sizeof(*xc) * DIM);
        totmass = m[plist[0]];
    } else {

        explicit_bzero(xc, sizeof(*xc) * DIM);

        for (int i = 0; i < n; i++) {
            int             q = plist[i];

            for (int j = 0; j < DIM; j++) {
                xc[j] += m[q] * x[q][j];
            }
            totmass += m[q];
        }

        for (int i = 0; i < DIM; i++)
            xc[i] /= totmass;

    }

    return totmass;
}

void calc_covariance(int n, int *plist, double xc[], double sigma2[DIM][DIM])
{
    double          totmass = 0;

    for (int i = 0; i < DIM; i++) {
        explicit_bzero(sigma2[i], sizeof(*sigma2[i]) * DIM);
    }

    if (n == 1) {
        return;
    }

    for (int i = 0; i < n; i++) {
        int             query = plist[i];

        for (int j = 0; j < DIM; j++) {
            for (int k = j; k < DIM; k++) {
                sigma2[j][k] += m[query] * x[query][j] * x[query][k];
            }
        }
        totmass += m[query];
    }

    /** normalize the partial covariance matrix */
    for (int i = 0; i < DIM; i++) {
        for (int j = i; j < DIM; j++) {
            sigma2[i][j] /= totmass;
        }
    }

    /** conclude the covariance matrix */
    for (int i = 0; i < DIM; i++) {
        /** process diagonal first **/
        sigma2[i][i] -= xc[i] * xc[i];
        /** upper triangle must be equal to lower trinagle: **/
        for (int j = i + 1; j < DIM; j++) {
            sigma2[j][i] = (sigma2[i][j] -= xc[i] * xc[j]);
        }
    }
}

/****** Under construction *****/

extern void release_sigma_octree_node(OCTANT *p);

void release_sigma_octree_node ( OCTANT* p )
{
    free(p->plist);
    p->parent = NULL;
    for (descriptor_t octant = 0; octant < NCHILDREN; octant++) {
        if (p->child[octant]) {
            free(p->child[octant]);
            p->child[octant] = NULL;
        }
    }
    free(p);
    p = NULL;
}

void sigma_octree_destroy(OCTANT *t)
{
    for (descriptor_t octant = 0; octant < NCHILDREN; octant++) {
        OCTANT *p = t->child[octant];
        if (p != NULL) {
            sigma_octree_destroy(p);
            release_sigma_octree_node(p);
        }
    }
    release_sigma_octree_node(t);
}
