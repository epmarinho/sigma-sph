/**@<sph-treecode_main.c>::**/

/********************************************************
 *  Author: Prof. Dr. Eraldo Pereira Marinho            *
 *                                                      *
 *  First working version on Mar 27 2012.               *
 ********************************************************/

#ifdef _FIND_SELFCONSISTENT_KNN_CLUSTER_
#    undef _ISOTROPIC_SPH_
#else
#    define _ISOTROPIC_SPH_
#endif

#ifdef _ISOTROPIC_SPH_
#    define _ISOTROPIC_KNN_INITIALIZATION_
#endif

/* overrule above definition if you do insist in having isotropic k-NN in initialization regardless anisotropic simulations: */
// #    define _ISOTROPIC_KNN_INITIALIZATION_

/* project headers * */
#include <sph-treecode_defaults.h>

/** define default number of omp threads (check the CPU specifications before using too many threads - I have head problems when setting high number of threads): **/
size_t          numthreads = 4;

double calc_knn_center_of_mass(const NNLIST * qnnlist, double xc[])
{
    double          totmass = 0;

    bzero(xc, sizeof(*xc) * DIM);
    int             i,
                    j;

#pragma omp parallel num_threads(numthreads)
    {
        double          aux_totmass = 0,
            aux_xc[DIM];

        bzero(aux_xc, sizeof(*xc) * DIM);

#pragma omp for
        for (i = 0; i < qnnlist->nn; i++) {
            index_t         query = qnnlist->list[i];

            aux_totmass += m[query];

            for (j = 0; j < DIM; j++) {
                aux_xc[j] += m[query] * x[query][j];
            }
        }   // end for

#pragma omp critical
        {
            totmass += aux_totmass;

            for (j = 0; j < DIM; j++) {
                xc[j] += aux_xc[j];
            }
        }
    }   // end parallel

    for (j = 0; j < DIM; j++) {
        xc[j] /= totmass;
    }

    return totmass;
}

void calc_knn_covariance(const NNLIST * qnnlist, const double xc[], double sigma2[DIM][DIM])
{
    double          totmass = 0;

    for (int i = 0; i < DIM; i++) {
        /** clear the entire row i **/
        bzero(sigma2[i], sizeof(*sigma2[i]) * DIM);
    }

#pragma omp parallel num_threads(numthreads)
    {
        double          aux_mass = 0;

        /** aux_sigma2 initialization **/
        double          aux_sigma2[DIM][DIM];

        for (int i = 0; i < DIM; i++) {
            /** clear the entire row i **/
            bzero(aux_sigma2[i], sizeof(*sigma2[i]) * DIM);
        }

#pragma omp for
        for (descriptor_t i = 0; i < qnnlist->nn; i++) {
            index_t         query = qnnlist->list[i];

            for (int j = 0; j < DIM; j++) {
                for (int k = j; k < DIM; k++) {
                    aux_sigma2[j][k] += m[query] * x[query][j] * x[query][k];
                }
            }

            aux_mass += m[query];

        }   // end for

#pragma omp critical
        {
            totmass += aux_mass;

            for (int j = 0; j < DIM; j++) {
                for (int k = j; k < DIM; k++) {
                    sigma2[j][k] += aux_sigma2[j][k];
                }
            }

        }   // end critical

    }   // end parallel

    for (int i = 0; i < DIM; i++) {
        for (int j = i; j < DIM; j++) {
            sigma2[i][j] /= totmass;
        }
    }

    for (int i = 0; i < DIM; i++) {
        /** compute diagonal **/
        sigma2[i][i] -= xc[i] * xc[i];
        for (int j = i + 1; j < DIM; j++) {
            /** compute upper and lower triangles **/
            sigma2[j][i] = (sigma2[i][j] -= xc[i] * xc[j]);
        }
    }
}

/** query is a shared variable to index SPH particles **/
/** WARNING: already defined in sph-treecode_knn.c */
// index_t         query;

/* Convergence-control parameters for the find_selfregulating_kNN_cluster method:
 * MAX_ITERATIONS_ALLOWED is a parameter that expresses the calm with which the learning
 * scheme awaits the presumed but not-always certain convergence */
#ifdef _FIND_SELFCONSISTENT_KNN_CLUSTER_
#    ifndef MAX_ITERATIONS_ALLOWED
#        define     MAX_ITERATIONS_ALLOWED      16
#    endif
#    define     MAXHISTORYSIZE              256
#else
    /** MAX_ITERATIONS_ALLOWED cannot be smaller than 1 **/
#    define     MAX_ITERATIONS_ALLOWED      1
#    define     MAXHISTORYSIZE              1
#endif

double          convtol,
                history[MAXHISTORYSIZE],
                previous_diff,
                oldprincipal;

/** the smoothing lengths are the smoothing tensor eigenvalues: **/
double          h[NMAX][DIM];
double          vol[NMAX];
flag_t          find_selfregulating_kNN_cluster_has_started = 0;

#define not !

/* convergence statistics variables */
double          mean_error = 0,
    mean_iterations = 0 /*,
                           mean_fault_iterations = 0 */ ;

counter_t       nfaults = 0;

/** Self-regulated, anisotropic kNN cluster: a fixed point problem */
void find_selfregulating_kNN_cluster(size_t n)
{
    /* The present is a machine-learning that searches for the optimal anisotropic, covariance based,
     * kNN cluster as follows:
     *           1. Initialize the metric tensor as the identity matrix.
     *           2. Perform the anisotropic kNN search, returning a cluster of the kNN particles from the query one.
     *           3. Get the principal components assuming that the kNN cluster's center of mass is the query itself.
     *           4. Repeat Step 2 until the the yielded kNN cluster becomes unchanged (fixed-point). 
     * WARNING: don't try to omp parallelize the present code neither aniso_k_NN_search */

    convtol = pow(K, (double) -2 / 3);  // / (double) (1 << 4);

    for (query = 0; query < n; query++) {

        if (isgas[query]) {

            flag_t          itsover /*, not_fault_yet = 1 */ ;

            descriptor_t    i,
                            j,
                            jj = -1;

            counter_t       hits = 0,
                stored = 0,
                step = 1,
                iterations = 0;

            flag_t          histoverflow = 0,
                ask_for_period_match = 0;

            double          knnx[DIM],
                            diff = 0,
                mindiff = INFINITY;

            if (not find_selfregulating_kNN_cluster_has_started) {
#ifdef _ISOTROPIC_KNN_INITIALIZATION_
                for (i = 0; i < DIM; i++) {
                    eigenval[query][i] = eigenvec[query][i][i] = 1;
                    for (j = i + 1; j < DIM; j++) {
                        eigenvec[query][i][j] = eigenvec[query][j][i] = 0;
                    }
                }
#else
                for (i = 0; i < DIM; i++) {
                    eigenval[query][i] = sigma_octree_root->eigenval[i];
                    memcpy(eigenvec[query][i], sigma_octree_root->eigenvec[i], sizeof(double) * DIM);
                }
#endif
            }

            for (oldprincipal = 0; 1;) {

                /* trying to reuse the k-nn list (zero initially) */
                for (i = 1; i < qnnlist.nn; i++) {

                    if (qnnlist.list[i] == query) {
                        index_t         tmp = qnnlist.list[0];

                        qnnlist.list[0] = query;
                        qnnlist.list[i] = tmp;

                        for (j = 1; j < qnnlist.nn; j++)
                            qnnlist.distance[j] = p2p_d(x[qnnlist.list[j]], x[qnnlist.list[0]], eigenval[query], eigenvec[query]);

                        /* binary-sort the just modified k-NN list: */
                        qnn_binsort(qnnlist.nn);

                        /* get the distance from the query to the outermost k-NN particle: */
                        qnnlist.topdistance = qnnlist.distance[j - 1];

                        break;

                    }   // end if (qnnlist.list[i] == query)

                }   //end for (i = 1; i < qnnlist.nn; i++)

                /* if it was no longer possible to reuse, reset the k-nn list */
                if (query == 0 || i == qnnlist.nn) {
                    qnnlist.list[0] = query;
                    qnnlist.nn = 1;
                    qnnlist.topdistance = 0;
                }   //end if (query == 0 || i == qnnlist.nn)

                /* find the anisotropic k-nn region given the present inverse-covariance eigenvalues/eigenvectors */
                aniso_k_NN_search(0, sigma_octree_root);
                /* bail out the learning-iteration if the anisotropic pattern remains (nearly) unchanged */
                iterations++;
                previous_diff = diff;
                diff = fabs(oldprincipal - eigenval[query][0]);
                double          noise = convtol * oldprincipal;

/** check if the iteration time is long enough or convergence tolerance has been reached **/
#define _TIGHT_CONVERGENCE_CRITERIA
#ifdef _TIGHT_CONVERGENCE_CRITERIA
                itsover =
            /** inserted on Fri Oct 8, 2021 to avoid filamentar ellipsoids **/
                    eigenval[query][0] > 3 * eigenval[query][1] ||
            /** **/
                    (iterations == MAX_ITERATIONS_ALLOWED) ||
                    (iterations > 2 &&
                     diff < noise && diff < previous_diff && fabs(previous_diff - diff) < convtol * (previous_diff + diff)
                    );
#else
                itsover =
                    iterations >= MAX_ITERATIONS_ALLOWED ||
                    (iterations > 2 && (diff < noise || fabs(previous_diff - diff) < convtol * (previous_diff + diff)));
#endif

                /*NOTE:******************************************************************
                 *  Aug 19, 2014, from 9:00 am through 11:00 pm.                        *
                 *  This is the new heuristic approach to try convergence               *
                 *  by means of a circularity-detection machine                         *
                 ***********************************************************************/
                /* check against convergence fault */
                flag_t          convergence_fault = query && not itsover && iterations > ceil(mean_iterations);

                if (convergence_fault) {
                    int             prevstored = stored;

                    nfaults++ /** = not_fault_yet; not_fault_yet = 0 **/ ;

                    for (i = 0; i < stored; i++) {
                        if (diff == history[i]) {
                            if (jj < 0 /* is it the first entry? */ ) {
                                jj = i;
                                mindiff = diff;
                            } else if (i != jj) {
                                mindiff = min(mindiff, diff);
                            } else {
                                ask_for_period_match = 1;
                            }
                            break;
                        }
                    }   //end for (i = 0; i < stored; i++)

                    if (i == stored) {
                        if (stored < MAXHISTORYSIZE) {
                            history[i] = diff;
                            stored++;
                        } else {
                            prevstored = 0;
                            if (histoverflow >= 0 && histoverflow <= MAXHISTORYSIZE) {
                                history[i - histoverflow] = diff;
                                histoverflow += step;
                            } else {
                                if (histoverflow > MAXHISTORYSIZE) {
                                    histoverflow--;
                                    step = -1;
                                } else {
                                    histoverflow++;
                                    step = 1;
                                }
                                hits++;
                            }
                        }
                    }   // end if (i == stored)

                    if (ask_for_period_match && prevstored == stored)
                        if (diff == mindiff) {
                            /* matched a period */ itsover = 1;
                        }
                }   //end if (!itsover && query && iterations > mean_iterations)

                if (itsover) {

                    /* renormalization of the eigenvalues to have the minimal ellipsoid that covers the k-NN cluster
                       -- this is the ellipsoid whose semi-major axes measure as the principal smoothing lengths */
                    for (vol[query] = 1, i = 0; i < DIM; i++) {
                        /** previously, eigenval[query][0..2] were the eigenvalues (sigma²) of the covariance tensor;
                            it is now the eigevalues of tensor associated to the minimal ellipsoid that includes all
                            the k-NN list **/
                        eigenval[query][i] *= qnnlist.topdistance; /** distances are square defined **/
                        /** now, the top non-square distance equals to 1 **/
                        /** Since May 4 2020, the hull ellipsoid is no longer the smoothing ellipsoid. Now the
                            latter is half the former to comply the commonly used in literature: */
                        eigenval[query][i] *= .25;
                        /** now, the top non-square distance equals to 2 if using eigenvector²/eigenval as in Mahalanobis metric **/
                        /** h[NMAX][DMI] are the eigenvalues of H: */
                        h[query][i] = sqrt(eigenval[query][i]);
                        /** h's are now half the semimajor axes of the kNN hull ellipsoid **/

                        /** vol[query] is the volume of the circumscribed parallelepiped to the hull ellipsoid, necessary to
                            compute the smoothing kernel in the module <sph-treecode_quantities.c> **/
                        vol[query] *= h[query][i];
                    }

#ifdef _FIND_SELFCONSISTENT_KNN_CLUSTER_
                    mean_error = (query * mean_error + fabs(diff) / oldprincipal) / (query + 1);
                    mean_iterations = (mean_iterations * query + iterations) / (query + 1);
#endif
                    break /** bail-out the for-loop **/ ;

                }   //end if (istover)

                /* redo the covariance analysis for the newly-found neighborhood (not-yet converged) */
                calc_knn_center_of_mass(&qnnlist, knnx);
                calc_knn_covariance(&qnnlist, knnx, eigenvec[query]);
                oldprincipal = eigenval[query][0];
                power_iteration_method(eigenval[query], eigenvec[query]);

            }   //end for (oldprincipal = 0; 1;)

            effneighb[query].nn = qnnlist.nn;
            memcpy(effneighb[query].list, qnnlist.list, qnnlist.nn * sizeof query);

            for (i = 0; i < qnnlist.nn; i++) {
                /** The distances gotta be doubled as the smoothing ellipsoid has been
                    reduced to half the k-NN ellipsoid hull (Monaghan's cubic B-spline kernel
                    vanishes for r = 2h which means that whenever the normalized distance equals
                    to one, it corresponds to twice the eigenvalues of the smoothing tensor
                    -- I don't like this kernel model) **/
                effneighb[query].distance[i] = 2 * sqrt(qnnlist.distance[i] / qnnlist.topdistance);
            }

        }   // end if (isgas[query])...

        //         fprintf(stderr, "thread %d down for query %d\n", thread_ID, query);

    }   //end parallel for (query = 0; query < n; query++)

    /** informe the algorithm to not reinitialize the exterior eigenvector/eigenvalue fields **/
    find_selfregulating_kNN_cluster_has_started = 1;

}   //end void find_selfregulating_kNN_cluster(size_t n)

/* system total mass */
double          Mtot;

double          xmean[DIM],
                R2 = 0,
    R;
flag_t          already_initiated_find_selfregulating_kNN_cluster;

void usage(char *command)
{
    fprintf(stderr, "Usage:\n%s <inputdata> <K> <DT> <theta> <epsilon> <numthreads>\n", command);
    fprintf(stderr, "where left to right occurrence must be preserved and <inputdata> is mandatory.\n");
}

int main(int argc, char *argv[])
{
#define MAXARGC 7
#ifdef _MAXARGC_REQUIRED_
    if (argc < MAXARGC) {
        fprintf(stderr, "%s: missing parameters; argc = %d whereas it should be %d\n", argv[0], argc, MAXARGC);
        exit(MISSARGS);
    }
#endif
    if (argc > MAXARGC) {
        fprintf(stderr, "%s: too many parameters!\n", argv[0]);
        usage(argv[0]);
        exit(OVERARGS);
    }
    if (argc < 2) {
        usage(argv[0]);
        exit(MISSFNAM);
    } else if ((fhandle = fopen(argv[1], "r")) == NULL) {
        fprintf(stderr, "%s: cannot read %s... exiting with error status\n", argv[0], argv[1]);
        usage(argv[0]);
        exit(MISSARGS);
    } else {

        int             ncols /** NOTE: ncols is by now obsolete **/ ,
                        n = getdimensions(fhandle, &ncols, &Mtot);

        fprintf(stderr, "\n\ninput file: %s\n", basename(argv[1])) /**1*/ ;
#ifndef __CALC_INTEGRALS_ONLY__
        double start_t, end_t;
#endif
        extern double omp_get_wtime(void);
        /* create an octree for the n-particle data set */
        //         fprintf(stderr, "mkoctree(sigma_octree_root = mkroot(n));\n");
        //         start_t = omp_get_wtime();
        sigma_octree_root = mkroot(n);
        mkoctree(sigma_octree_root);
        //         fprintf(stderr, "done in %lg sec\n", (end_t = omp_get_wtime()) - start_t);

    /********************************* simulation-parameters default: ****************************************/

        for (index_t particle = 0; particle < n; particle++) {
            for (descriptor_t j = 0; j < DIM; j++) {
                xmean[j] += m[particle] * x[particle][j];
            }
            R2 += m[particle] * dotproduct(x[particle], x[particle]);
        }
        for (descriptor_t j = 0; j < DIM; j++) {
            xmean[j] /= Mtot;
        }
        R2 /= Mtot;
        R2 -= dotproduct(xmean, xmean);

        /** the old Aarseth's softening length *//*epsilon2 = max(.5 * R2 / pow((double) n, (double) 2 / (double) 3), 1e-6); */

            /**********************************************************
             * GM = R³/DT² => DT = sqrt(GM R³)                        *
             *      DT = sqrt (Mtot/n * pow(epsilon2, (double) 3/2)); *
             **********************************************************/

        /** major time-step: *//*DT =
           sqrt((Mtot / (double) n) * pow(epsilon2, 1.5));
           double          DT_sav = DT;

           for (DT = 1; DT > DT_sav; DT *= .5) {
           ;
           } */

        /** suggested number k of nearest neighbors (K) based upon Poisson counting-error: */
        K = .5 * ceil(sqrt((double) n));
        /** suggested root time-step **/
        DT = 1. / 1024.;
        /** suggested theta **/
        theta2 = .25;
        theta2 *= theta2;
        /** suggested epsilon **/
        epsilon = 0.066;
        epsilon2 = epsilon * epsilon;

    /************* set-up the simulation parameters to override the previously suggested ones ****************/
        double          daux;
        int             iaux;

        switch (argc) {
        case 7:
            iaux = atoi(argv[6]);
            if (iaux > 0) {
                numthreads = iaux;
            } else {
                fprintf(stderr, "using default number of threads\n");
            }
        case 6:
            daux = atof(argv[5]);
            if (daux > 0) {
                epsilon = daux;
                epsilon2 = epsilon * epsilon;
            }
        case 5:
           /** squared tolerance parameter **/
            daux = atof(argv[4]);
            if (daux >= 0) {
                theta2 = daux;
                theta2 *= theta2;
            }
        case 4:
           /** root time-step **/
            daux = atof(argv[3]);
            if (daux > 0) {
                DT = daux;
            }
        case 3:
           /** number K of nearest neighbors **/
            if ((iaux = atoi(argv[2])) > 2) {
                if (iaux > MAXNN - 1) {
                    fprintf(stderr,
                            "%s: number of nearest neighbors (%d) exceeds maximum number of particles minus 1 (%d)\n", argv[0], K,
                            MAXNN - 1);
                    exit(-3);
                }
                K = iaux;
            } else {
                fprintf(stderr, "%s: number of nearest neighbors must be greater than or equal to 3. Using default K = %d\n",
                        argv[0], K);
            }
        }

#define _VERBOSE_INPUT_
#ifdef _VERBOSE_INPUT_
        fprintf(stderr, "number of threads = %zu\n", numthreads);
        fprintf(stderr, "n = %d\n", n) /**2*/ ;
        fprintf(stderr, "K = %d\n", K) /**3*/ ;
        fprintf(stderr, "DT = %lg = 1/%lg\n", DT, 1 / DT) /**4*/ ;
        fprintf(stderr, "theta = %lg\n", sqrt(theta2)) /**5*/ ;
        fprintf(stderr, "epsilon = %lg\n", epsilon) /**6*/ ;
        fprintf(stderr, "Total mass = %lg\n", Mtot) /**7*/ ;
        fprintf(stderr, "mean radius (rms) = %lg\n", R = sqrt(R2)) /**8*/ ;
        fprintf(stderr, "sqrt(G M / R) = %lg\n", sqrt(Mtot / R)) /**9*/ ;
#endif
#ifdef _OPENMP
        /** set default number of threads (redundant) **/
        omp_set_num_threads(numthreads);
#endif

#ifndef __CALC_INTEGRALS_ONLY__
        /***********************************************************************************
         * NOTE: therein is my not-yet published anisotropic-kNN cluster detection scheme, *
         * which I have temporarily named self-consistent kNN cluster                      *
         ***********************************************************************************/

        fprintf(stderr, "find_selfregulating_kNN_cluster(n);\n");
        start_t = omp_get_wtime();
        find_selfregulating_kNN_cluster(n);
        fprintf(stderr, "done in %lg sec\n", (end_t = omp_get_wtime()) - start_t);

#    if defined(_VERBOSE_INPUT_) && defined(_FIND_SELFCONSISTENT_KNN_CLUSTER_)
        fprintf(stderr, "\nmaximum iterations allowed = %d\n", MAX_ITERATIONS_ALLOWED);
        fprintf(stderr, "mean iterations = %lg\n", mean_iterations);
        fprintf(stderr, "total convergence faults = %d\n", nfaults);
        fprintf(stderr, "fraction of particles in convergence fault = %lg\n", nfaults / (double) n);
        fprintf(stderr, "mean convergence error = %lg %%\n\n", mean_error * 100);
#    endif
        //         fprintf(stderr, "init_computer_physical_scales();\n");
        //         start_t = omp_get_wtime();
        init_computer_physical_scales();
        //         fprintf(stderr, "done in %lg sec\n", (end_t = omp_get_wtime()) - start_t);
        already_initiated_find_selfregulating_kNN_cluster = 1;

        /*********************************************************************************
         * WARNING:                                                                      *
         * The following procedure has a number of issues therein, which may compromise  *
         * the overall simulation results                                                *
         *********************************************************************************/
        //         fprintf(stderr, "sph_init(n);\n");
        //         start_t = omp_get_wtime();
        sph_init(n);
        //         fprintf(stderr, "done in %lg sec\n", (end_t = omp_get_wtime()) - start_t);
        //         fprintf(stderr, "timebin_share(n);\n");
        //         start_t = omp_get_wtime();
        timebin_share(n);
        //         fprintf(stderr, "done in %lg sec\n", (end_t = omp_get_wtime()) - start_t);
        fprintf(stderr, "binary_leapfrog(0);\n");
        start_t = omp_get_wtime();
        binary_leapfrog(0);
        fprintf(stderr, "done in %lg sec\n", (end_t = omp_get_wtime()) - start_t);
        //         fprintf(stderr, "show_result(n /*, ncols */ );\n");
        //         start_t = omp_get_wtime();
        show_result(n /*, ncols */ );
        //         fprintf(stderr, "done in %lg sec\n", (end_t = omp_get_wtime()) - start_t);
        //         fprintf(stderr, "integral_quantities(n);\n");
        //         start_t = omp_get_wtime();
#else
        extern void     find_well_distant_nodes(int particle, const OCTANT * cluster);

        for (index_t particle = 0; particle < n; particle++)
            find_well_distant_nodes(particle, sigma_octree_root);
#endif
        integral_quantities(n);
        //         fprintf(stderr, "done in %lg sec\n", (end_t = omp_get_wtime()) - start_t);

        return 0;
    }
}

/* Last changes:

   Thu Mar 15 14:20:31 GMT -03, 2018
   Fri Mar 16 15:10:29 GMT -03, 2018
   Tue Mar 27 16:19:00 GMT -03, 2018
   Mon Mar 23 16:43:11 GMT -03, 2020
   Thu 16 Apr 2020 11:22:13 AM -03
   Fri Apr 23 12:17 GMT -03, 2021
   Thu Apr 29 09:53:01 AM -03 2021

 */
