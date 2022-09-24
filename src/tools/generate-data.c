#include <stdio.h>
#include <stdlib.h>

void
generatedata(int n)
{
    int             i;
    for (i = 0; i < n; i++) {
        int             j;
        printf("%24.16le", (1 - .5 * drand48()) / n);
        for (j = 0; j < 3; j++) {
            printf("%24.16le", drand48() - .5);
        }
        for (j = 0; j < 3; j++) {
            printf("%24.16le", .1 * (drand48() - .5));
        }
        printf("\n");
    }
}

main(int argc, char const *argv[])
{
    generatedata(atoi(argv[1]));
    exit(0);
}
