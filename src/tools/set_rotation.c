#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <ctype.h>
#include <octree.h>


double          mass[N];
double          x[N][DIM];
double          v[N][DIM];
double          g[N][DIM];

int             disallow[MAXNODES];

double
distance(double x[DIM], double y[DIM])
{
    int             j;
    double          r;
    for (r = j = 0; j < DIM; j++) {
        double          a = x[j] - y[j];
        r += a * a;
    }
    return r;
}

double          Omega = .25;

int
faraway(int i, int node)
{
    if (Octree[node].identifier) {
        return Octree[node].identifier - 1 != i;
    } else {
        double          r = distance(x[i], Octree[node].vertex);
        return Omega * r >= Octree[node].size;
    }
}

int             welldistant[N];
int             accepted = 0;

void
Octree_descent(int i)
{
  /**
   * As an application, we are performing a topdown search for
   * the nearest neighbors close enough to be encompassed by
   * a solid angle Omega. The main interpretation for this kind
   * of NN search is the fast calculation of the gravity forces
   * on N-body simulations.
   */
    int             level;

    bzero(disallow, MAXNODES * sizeof(int));

    bzero(welldistant, N * sizeof(int));

    accepted = 0;

    for (level = 1; level < depth; level++) {
        int             node;

    /** forall node in this level, do: **/
        for (node = levelbegin[level]; node < levelend[level]; node++) {
            disallow[node] = disallow[node] ||
                (faraway(i, node) && (welldistant[accepted++] = node));
        }                       // end-forall

    /** forall node in this level, do: **/
        for (node = levelbegin[level]; node < levelend[level]; node++) {
            int             q;

      /** forward node disallowance to its children **/
            for (q = 0; q < GRD; q++) {
                if (Octree[node].child[q]) {
                    disallow[Octree[node].child[q]] = disallow[node];
                }
            }

        }                       // end-forall

    }                           // end-for level
}

double          epsilon = .04;

void
octree_gravity(int n)
{
    int             i,
                    k;

    for (i = 0; i < n; i++) {
        int             k;
        double          r,
                        vc;

        printf("%24.16le", mass[i]);

        for (k = 0; k < DIM; k++) {
            printf("%24.16le", x[i][k]);
        }

    /** perform the topdown search for the nearest neighbors seen inside
     * the solid angle Omega */
        Octree_descent(i);

        for (bzero(g[i], sizeof(double) * DIM), k = 0; k < accepted; k++) {
            int             nu = welldistant[k],
                j;
            double          r =
                distance(Octree[nu].vertex, x[i]) + epsilon;

            for (j = 0; j < DIM; j++) {
                double          delta = Octree[nu].vertex[j] - x[i][j];
                g[i][j] += Octree[nu].mass * delta / (r * sqrt(r));
            }

        }

        r = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        vc = sqrt(g[i][0] * g[i][0] + g[i][1] * g[i][1]);
        vc = sqrt(vc / r);

        v[i][0] += vc * (-x[i][1]);
        v[i][1] += vc * (x[i][0]);

        for (k = 0; k < DIM; k++) {
            printf("%24.16le", v[i][k]);
        }
#ifdef _SAVE_GRAVITY_FIELD_
        for (k = 0; k < DIM; k++) {
            printf("%24.16le", g[i][k]);
        }
#endif

        printf("\n");
    }
}

int
skipspaces(FILE * input)
{
    int             c;
    while (isspace(c = getc(input)));
  /** check eof **/
    if (c == -1)
        return -1;
    ungetc(c, input);
    return 0;
}

int
getdata(int n, FILE * input)
{
    int             j,
                    i;
    for (i = 0; 1; i++) {
        int             c;

        skipspaces(input);
        fscanf(input, "%le", &mass[i]);

        for (j = 0; j < DIM; j++) {
            skipspaces(input);
            fscanf(input, "%le", &x[i][j]);
        }

        for (j = 0; j < DIM; j++) {
            skipspaces(input);
            fscanf(input, "%le", &v[i][j]);
        }

    /** seek eol **/
        while ((c = getc(input)) != '\n');

        if (skipspaces(input))
            return i + 1;
    }
}

main(int argc, char const *argv[])
{
    int             n = N,
        i;
    FILE           *input;

    if (argc > 1) {
        n = atoi(argv[1]);
        input = fopen(argv[1], "r");
    } else
        input = stdin;
    if (argc > 2) {
        Omega = atof(argv[2]);
    }
    if (argc > 3) {
        epsilon = atof(argv[3]);
    }

    Octree_build(n = getdata(n, input));

    octree_gravity(n);

    exit(0);
}
