/**@<sph-treecode_knn.c>::**/
/* This is a modification to Marinho (1991, 1997) recursive method for
 * K-NN search.
 * This is the last working version since Mar 27 2012. */
#include <sph-treecode_defaults.h>

NNLIST          qnnlist;

index_t         query,
                K;

void aniso_k_NN_search(double d, OCTANT const *p)
{
    /* The pointer p abstracts a 3D convex region bounded by the 3 principal planes,
     * coming out from the principal component analysis done over the parent cluster */
    if (p->n > 1) {
        /* the current cluster is not a singleton */
        int             q;
        OCTANT         *child;

        /* NCHILDREN equals 2^3 in 3D spaces */
        for (q = 0; q < NCHILDREN; q++) {
            /* perform the search through the next non-empty child */
            if ((child = p->child[q]) != NULL) {

                /* take the query-to-child anisotropic distance */
                double          r = query2node_distance(d, query, child);

                /* check whether the outermost K is away from the new neighbor found or
                 * if the K list is not-yet saturated */
                if (r < qnnlist.topdistance || qnnlist.nn < K) {
                    aniso_k_NN_search(r, child);
                }
            }
        }   // end for
    } else {
        /* insert the newly found neighbor by pumping out the outer neighbors; the outermost neighbor
         * drops out from the list if the population has reached the K upper limit */
        if (isgas[p->plist[0]]) {
            pump_neighbor(query, p->plist[0]);
        }
    }
}

double query2node_distance(double d, int query, OCTANT const *p)
{
  /** d is the inherited distance **/
  /** D is the face-on distance **/
    double          D = 0;

    int             k;

    for (k = 0; k < DIM; k++) {
        double          a = 0,
            b = 0;

    /** the parent's eigenvectors form a set of projectors onto the local
     * principal directions **/
        int             j;

        for (j = 0; j < DIM; j++) {

            a += (x[query][j] - p->parent->x[j]) * p->parent->eigenvec[k][j];

            b += (p->parent->x[j] - p->x[j]) * p->parent->eigenvec[k][j];

        }

        /*
         * if the product a*b is positive then the octant's lower boundary
         * is facing the query point
         */
        if ((a * b) > 0) {
     /** point-to-plane Mahalanobis distance **/
            D = max(D, a * a / eigenval[query][k]);
        }
    }

    return max(D, d);
}

/* compute point-to-point Mahalanobis distance:
 *
 * Distances are calculated according to the Mahalanobis metric. However, the original
 * context is slightly distracted, so that the peer-to-peer distance is induced by a
 * metric tensor, which is precisely the covariance tensor evaluated on the K cluster
 * centered in the query particle.
 *
 * NOTE:
 * Both eigenval[] and eigenvec[] are training aspect vectors.
 * The information retrieved is the anisotropic K-list under the morphological
 * information given by eigenval[] and eigenvec[].
 */
double p2p_distance(int i, int query)
{

    int             j,
                    k;

    double          sqrd_x_iq,
                    a_k;

    for (sqrd_x_iq = k = 0; k < DIM; k++) {

        for (a_k = j = 0; j < DIM; j++) {

            double          x_iq_j = (x[query][j] - x[i][j]);

            a_k += x_iq_j * eigenvec[query][k][j];
        }

        sqrd_x_iq += (a_k * a_k) / eigenval[query][k];
    }

    return sqrd_x_iq;
}

void pump_neighbor(int query, int i)
{
    double          r = p2p_distance(i, query);

    if (r == 0)
        /* of course i is the query itself */
        return;

    if (qnnlist.nn < K) {
        /* in this case, the knn-list is not yet saturated */

        /* search for the first K neighbors whatever they are */
        int /*j, */     last = qnnlist.nn;

        qnnlist.list[last] = i;
        qnnlist.distance[last] = r;
        qnnlist.nn++;

        if (qnnlist.nn == K) {
            //             qnn_binsort(qnnlist.nn);
            //             qnn_sort_bruteforce();
//             qnn_quicksort(0, qnnlist.nn - 1);
            qnn_quickSort(0, qnnlist.nn - 1);
            qnnlist.topdistance = qnnlist.distance[last];

        }

    } else if (r < qnnlist.topdistance) {

        /* search for the K nearest neighbors along the remaining nodes to
         * descend through */
        int             j,
                        k,
                        nn,
                        last = (nn = qnnlist.nn) - 1;

        for (j = 1; j < nn; j++) {
            if (r < qnnlist.distance[j])
                break;
        }

        if (i == qnnlist.list[j - 1])
            return;

        for (k = last; k > j; k--) {
            qnnlist.list[k] = qnnlist.list[k - 1];
            qnnlist.distance[k] = qnnlist.distance[k - 1];
        }

        qnnlist.list[k] = i;
        qnnlist.distance[k] = r;
        qnnlist.topdistance = qnnlist.distance[last];

    }
}

/* This is the troublesome part of the present anisotropic SPH method.
 * After the K it is required that the K-list includes also the particles
 * that give non-zero contribution to the kernel, aka effective neighbors.
 * To avoid replicating K search and then saving time, we can reuse the
 * K list to add only the particles which are out from the previous K-list
 * but their specific smoothing radius must include the current query.
 */

ENNLIST         effneighb[NMAX];

#ifdef _USE_EFFECTIVE_NEIGHBORS_SEARCH_ // deprecated use
void add_neighbor(int query, int i)
{
    if (*eigenval[query] >= *eigenval[i]) {
        /* the query radius is either already covering the
         * effective neighbor candidate, or the candidate has a small
         * radius and it is away from the query range.
         */
        return;
    } else {
        double          r = p2p_distance(i, query);

        if (r < *eigenval[i]) {
            /* at this point, i-particle is an effective neighbor since its
             * smoothing radius h[i] encompasses the query particle.
             */
            effneighb[query].list[effneighb[query].nn++] = i;
        }
    }
}

void effective_neighbors_search(double d, OCTANT * p)
{
    if (p->n > 1) {
        int             q;
        OCTANT         *child;

        for (q = 0; q < NCHILDREN; q++) {
            if (child = p->child[q]) {
                double          r = query2node_distance(d, query, child);

                if (child->n == 1 || 1 <= r && r <= child->topdistance) {
                    effective_neighbors_search(r, child);
                }
            }
        }
    } else {
        add_neighbor(query, p->plist[0]);
    }
}
#else
#    define _USE_SYMMETRIC_CLOSURE_
                                // for posterior usage
/* symmetric_closure introduced in Aug 18, 2014, after the original code knn.c at repository path:
 * sources/SPH/sph\ 2.6/src/ */
void symmetric_closure(size_t n)
{
    descriptor_t    k;

    /** for each i-particle, do the search: **/
    for (index_t i = 0; i < n; i++) {

        /** visit all the effneighb[i].nn particles pointed by effneighb[i].list[j] observing 
         *  that the list is growing along the search: **/
        for (descriptor_t j = 1; j < effneighb[i].nn; j++) {

            index_t         l = effneighb[i].list[j];

            /** search for the first match effneighb[l].list[k] == i or exhaust the search 
             *  until k == effneighb[l].nn **/
            for (k = 1; k < effneighb[l].nn && effneighb[l].list[k] != i; k++);

            /** if the search above was unmatched, k == effneighb[l].nn, then add i to the l's kNN list **/
            if (k == effneighb[l].nn && effneighb[l].nn < MAXENN) {

                /** add the query i as the new symmetric neighbor of particle l: **/
                effneighb[l].list[effneighb[l].nn] = i;
                /** i-to-l distance must be computed from the l's POV: **/
                effneighb[l].distance[effneighb[l].nn] = sqrt(p2p_d(x[l], x[i], eigenval[l], eigenvec[l]));
                /** count the new appended neighbor **/
                effneighb[l].nn++;

            }

        }
    }
    /** buble-sort the just formed list since new entries may not be ordered by anisotropic distances: **/
    //     for (index_t i = 0; i < n; i++) {
    //         size_t nn = effneighb[i].nn;
    //         for (descriptor_t j = 1; j < nn - 1; j++) {
    //             double dj = effneighb[i].distance[j];
    //             index_t lj = effneighb[i].list[j];
    //             for (descriptor_t k = j + 1; k < nn; k++) {
    //                 double dk = effneighb[i].distance[k];
    //                 if (dj > dk) {
    //                     effneighb[i].distance[j] = dk;
    //                     effneighb[i].distance[k] = dj;
    //                     effneighb[i].list[j] = effneighb[i].list[k];
    //                     effneighb[i].list[k] = lj;
    //                 }
    //             }
    //         }
    //     }
}
#endif
