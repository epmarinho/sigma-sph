#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

int
main(int argc, char const *argv[])
{
    if (argc < 8) {
        fprintf(stderr, "missing arguments\n");
        exit(-1);
    } else {
        FILE           *data = fopen(argv[1], "r");
        if (data == NULL) {
            fprintf(stderr, "cannot open %s\n", argv[1]);
            exit(-2);
        } else {
            double
                            X = atof(argv[2]),
                Y = atof(argv[3]),
                Z = atof(argv[4]),
                Vx = atof(argv[5]),
                Vy = atof(argv[6]),
                Vz = atof(argv[7]),
                xm = 0,
                ym = 0,
                zm = 0,
                vxm = 0,
                vym = 0,
                vzm = 0,
                M = 0;
            int             n = 0;

            while (1) {
                int             c;
                double          m,
                                x,
                                y,
                                z,
                                vx,
                                vy,
                                vz;
                while (isspace(c = getc(data)));
                if (c == -1)
                    break;
                ungetc(c, data);
                fscanf(data, "%le %le %le %le %le %le %le",
                       &m, &x, &y, &z, &vx, &vy, &vz);
                xm += x * m;
                ym += y * m;
                zm += z * m;
                vxm += vx * m;
                vym += vy * m;
                vzm += vz * m;
                M += m;
                while ((c = getc(data)) != '\n');
                n++;
            }
            xm /= M;
            ym /= M;
            zm /= M;
            vxm /= M;
            vym /= M;
            vzm /= M;
            rewind(data);
            for (; n; n--) {
                int             c;
                double          m,
                                x,
                                y,
                                z,
                                vx,
                                vy,
                                vz;
                while (isspace(c = getc(data)));
                if (c == -1)
                    break;
                ungetc(c, data);
		double dt, u, udot, rhodot, epsilon; int isgas;
                /**           m   x   y   z   vx  vy  vz  u   dt  udot rhodot isgas epsilon**/
                fscanf(data, "%le %le %le %le %le %le %le %le %le %le %le %d %le", &m, &x, &y, &z, &vx, &vy, &vz, &u, &dt, &udot, &rhodot, &isgas, &epsilon);
                /**                m         x        y        z       vx       vy       vz       dt       u        udot    rhodot**/
                fprintf(stdout, "%20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %d %20.16le\n",
                        m, x - xm + X, y - ym + Y, z - zm + Z, vx - vxm + Vx, vy - vym + Vy, vz - vzm + Vz,  u, dt, udot, rhodot, isgas, epsilon);
                while ((c = getc(data)) != '\n');
            }
        }
    }
    return 0;
}
