/** @<sph-treecode_quantities.c>::**/

#define not !

#include <sph-treecode_defaults.h>

extern void     explicit_bzero(void *s, size_t n);

flag_t          deleted[NMAX];

double dotproduct(double x[], double y[])
{
    double          s = 0;

    for (descriptor_t j = 0; j < DIM; j++) {
        s += x[j] * y[j];
    }
    return s;
}

/** Astronomical parameters and constants -- check them out again and again! */
double          MSun = 1.9891e+33 /** the sun mass in grams */ ;
double          RSun = 6.9634e5 /** solar radius in  cm */ ;
double          Mearth = 5.972e+27 /** the earth mass in grams */ ;
double          AU = 1.49597871e+13 /** astronomical unity in cm */ ;
double          Oort = 1e5 * 1.49597871e+13 /** inner radius of the Oort region cm */ ;
double          pc = 3.08567758e+18 /** parsec in centimeters (see IAU definition) */ ;
double          ly = 3.08567758e+18 / 3.2616
/** a light-year in cm -- a parser corresponds to 3.2616 light-years */
    ;
double          Myr = 1e6 * 365.25 * 3600 * 24 /** 1 million years in seconds */ ;
double          MGalaxy = 1.5e+12 * 1.9891e+33 /** estimated mass of our galaxy, the Via Lactea, in grams */ ;
double          RGalaxy = 15 * 1000 * 3.08567758e+18
/** galactic radius in cm, see https://en.wikipedia.org/wiki/Milky_Way */ ;
double          zGalaxy = 0.3 * 3.08567758e+18 * 1000
/** estimated Galaxy disc thickness in cetimeters, see https://en.wikipedia.org/wiki/Milky_Way */
    ;
double          TSun = 240 * 1e6 * 365.25 * 23.9344699 * 3600
/** estimated period of transit of the sun around the Galaxy center in seconds */
    ;
double          R0 = /* 8.34 kpc = */ 8.34 * 1000 * 3.08567758e+18
/** this is the solar orbit radius in cm, see https://en.wikipedia.org/wiki/Milky_Way */ ;

/** Fundamental physical constants: */
double          c = 29979245800 /* the speed of the light in cm s⁻¹ */ ;
double          G = 6.67428e-8 /* the universal gravitation constant in cm³g⁻¹s⁻² */ ;
double          Rgas = 8.314e+07 /* erg mole⁻¹ Kelvin⁻¹ */ ;
double          NA = 6.0221415e+23 /* Avogadros number in particles per mole */ ;

#ifndef H_ATOM
#    define H2_MOLECULE
#endif

/** standard density at normal conditions: */
//double          rho0 = 2 /*g mole⁻¹ for H2 */  / 22.414e+03 /*cm³ mole⁻¹ at 1 atm 0⁰C */ ;
/** Oxygen molecule mixture */
// rho0 = 32 / 22.414e+03 /* g cm⁻³ for O2; */
double          P0; // = /* 1 atm = */ 1013250 /* g cm⁻¹ s⁻² */ ;

/** Typical ISM off-core density */
#ifdef H2_MOLECULE
double          rho0 = 1 /*H2 molecules per cm³ */  * 2 /*g mole⁻¹ */  / 6.0221415e+23 /*molecules mole⁻¹ */ ;
#else
double          rho0 = 1 /*H atoms per cm³ */  * 1 /*g mole⁻¹ */  / 6.0221415e+23 /*molecules mole⁻¹ */ ;
#endif

/** computational speed of the light: */
double          c_norm;

double          T_cmb = 2.725 /* Kelvin */ ;

/** adiabatic index **/
double          gamma_index;

/** The reciprocal heat capacity at constant volume for an ideal diatomic gas, 1/Cv = 2/phi: */
#define POINTFOUR	0.4
#define TWOTHIRDS	0.6666666666666667
double          oneoverCv[NMAX];

#ifdef H2_MOLECULE
double          Oneovercv = POINTFOUR,
    Cv = 2.5;
#else
double          Oneovercv = TWOTHIRDS,
    // presumed H-atom gas
    Cv = 1.5;
#endif

/** NOTE: Cp = 1 + Cv, gamma_index = Cp / Cv = 1 + 1 / Cv, Cv = 5/2 for diatomic. **/

double          H2 = 2 / 6.0221415e+23,
    H = 1 / 6.0221415e+23,
    O2 = 32 / 6.0221415e+23;

/** Simulation adimensional unities: */
double          l_unit;
double          t_unit;
double          v_unit;
double          m_unit;
double          vol_unit;
double          rho_unit;
double          P_unit;
double          u_unit;
double          T_unit;

double          u0,
                LabTemp,
                u_ism,
                u_cmb;

double          meanweight =
#ifdef H2_MOLECULE
    2;
#else
    1;
#endif

    /** do not uncomment **/
double          P00 /* = P0 / P_unit */ ;

    /** do not uncomment **/
double          rho00 /* = rho0 / rho_unit */ ;

double          u00 /* = u0 / u_unit */ ;

void init_computer_physical_scales(void)
{
    //     l_unit = 39.48 * AU /** Pluto's semi-major axis */ ;
    //     l_unit = 1.0 * AU;
    //     l_unit = 1.0 * Rsun;
    //     l_unit = 1000 * AU;
    //     l_unit = R0; /* this to perform Galaxy simulation having the unit length of the Sun-to-centre distance */
    l_unit = RGalaxy;   /* this to perform Galaxy simulation */
    //     l_unit = 1000 * pc /* near-galactic length scale */;
    //     l_unit = 25000 * pc /* galactic length scale */;
    vol_unit = pow(l_unit, 3);

    //     m_unit = 1.0014 /* == present-day solar system total mass */  * MSun;
    //     m_unit = 1.0 * MSun;
    //     m_unit = 222.22 * MSun;
    m_unit = MGalaxy;   /* this to perform Galaxy simulation */

    /** gravity acceleration is given by G M / r² in l_unity * t_unity⁻² so that
        G = 6.67428e-8 cm³g⁻¹s⁻² == 1 l_unit³ m_unit⁻¹ t_unit⁻² thus, */

    //     t_unit = TSun;

    //     t_unit = 100 * Myr /* this is a galactic time scale */;
    //     m_unit = 1.3 * MSun;

    //     vol_unit = G * m_unit * pow(t_unit, 2);
    //     l_unit = pow(vol_unit, 1.0 / 3.0);

    //     m_unit = vol_unit / (G * pow(t_unit, 2));   /* mass unit is estimated from gravity free-fall scale */
    t_unit = sqrt(vol_unit / (G * m_unit)); /* time unit is estimated from gravity free-fall scale */

    /** velocity unit should be at least 1 km s⁻¹ for solar system simulation **/
    v_unit = l_unit / t_unit;

    /** One important derivative unit is the density unit: */
    rho_unit = m_unit / vol_unit;

    /** pressure unit: */
    //    P_unit = m_unit / (pow(t_unit, 2) * l_unit);
    //    P_factor = P_unit / P0;
    //    u0 /** erg g⁻¹ */  = Cv * Rgas    /** erg mole⁻¹ Kelvin⁻¹ */
    //    * 273.15 /** Kelvin */  / 2 /** g mole⁻¹ for H2 */ ;
    /** We redefine u0 to fit the commonly measured ISM temperature ~10K: */
    u0 /** erg g⁻¹ */  = Cv * Rgas /** erg mole⁻¹ Kelvin⁻¹ */  * 10 /** Kelvin */  / meanweight;
    /** defining P0 as the typical H2 pressure of the ISM **/
    P0 = rho0 / Cv * u0;
    /** P0 = rho0 / Cv * Cv * Rgas * 10 / meanweight; **/
    /** P0 = rho0 * Rgas * 10 / meanweight; **/

    /** specific thermal energy unit: */
    u_unit = pow(v_unit, 2);
    P_unit = u_unit * rho_unit;
    // u_unit = P_unit / rho_unit;
    /** temperature in Kelvin given u = 1 [u] */
    T_unit = meanweight * u_unit / (Cv * Rgas);
    //
    /** the computational speed of the light shall be useful in a near future: */
    c_norm = c / v_unit;
    u_ism /** erg g⁻¹ */  =
        Cv * Rgas /** erg mole⁻¹ Kelvin⁻¹ */  * 10 /** Kelvin */  / meanweight /** g mole⁻¹ for H2 */ ;
    /** Cosmic background radiation thermalized to H2 molecules: */
    u_cmb /** erg g⁻¹ */  = Cv * Rgas * T_cmb / meanweight;

    P00 = P0 / P_unit;
    rho00 = rho0 / rho_unit;
    u00 = u0 / u_unit;

    // #define _VERBOSE_PHYS_UNITS_
#ifdef _VERBOSE_PHYS_UNITS_
    extern double   Mtot,
                    R;
    double          Tvirial = .5 * G * 2 /*g mole⁻¹ */  * (Mtot * m_unit) / (R * l_unit) / (Cv * Rgas);
    double          radperyear = 3600 * 24 * 365.25 / t_unit,
        rotperyear = radperyear / (2 * pi);

    fprintf(stderr, "\n[m] = %lg g = %lg solar masses\n", m_unit, m_unit / MSun);
    fprintf(stderr, "[l] = %lg cm = %lg pc = %lg AU\n", l_unit, l_unit / pc, l_unit / AU);
    fprintf(stderr, "[t] = %lg s = %lg Myr = %lg yr\n", t_unit, t_unit / Myr, t_unit / Myr * 1e6);
    fprintf(stderr, "[v] = %lg cm s⁻¹ = %lg km s⁻¹\n", v_unit, v_unit * 1e-5);
    fprintf(stderr, "[u] = %lg erg g⁻¹\n", u_unit);
    fprintf(stderr, "[P] = %lg erg cm⁻³\n", P_unit);
    fprintf(stderr, "[rho] = %lg g cm⁻³ = %lg molecules cm⁻³\n", rho_unit, rho_unit
#    ifdef H2_MOLECULE
            / 2
#    endif
            * NA);
    fprintf(stderr, "temperature under 1 [u] = %lg Kelvin\n", T_unit);
    fprintf(stderr, "rho0 = %lg g cm⁻³ = %lg particles cm⁻³ = %lg [rho] \n", rho0, rho0
#    ifdef H2_MOLECULE
            / 2
#    endif
            * NA, rho0 / rho_unit);
    fprintf(stderr, "u0 = %lg erg g⁻¹ = %lg [u]\n", u0, u0 / u_unit);
    fprintf(stderr, "P0 = %lg erg cm⁻³ = %lg [P]\n", P0, P0 / P_unit);
    fprintf(stderr, "[Omega] = %lg rad yr⁻¹ = %lg rot yr⁻¹\n", radperyear, rotperyear);
    fprintf(stderr, "computational speed of the light, c_norm = %lg [l][t]⁻¹\n", c_norm);
    fprintf(stderr, "temperature on virial = %lg K = %lg ⁰C\n", Tvirial, Tvirial - 273.15);
    fprintf(stderr, "softening length = %lg cm = %lg AU\n", epsilon * l_unit, epsilon * l_unit / AU);
#endif

}   // end void init_computer_physical_scales(void)

double          cs[NMAX];

double          mu_max[NMAX];

double          rhocrit;

void sph_init(size_t n)
{
    size_t          n_gas = 0;

#pragma omp parallel for reduction(+:n_gas)
    for (descriptor_t i = 0; i < n; i++)
        n_gas += (isgas[i] != 0);

    //fprintf(stderr, "entering sph_init\n");
    /** search for effective neighbors; namely, those ones that return nonzero contribution for the symmetrizing kernel*/
    symmetric_closure(n_gas);

    /** do the anisotropic density estimation: */
    for (index_t particle = 0; particle < n; particle++)
        if (isgas[particle]) {

            rho[particle] = 0;
            for (descriptor_t l = 0; l < effneighb[particle].nn; l++) {

                /*j is the effective-neighbor identifier */
                index_t         j = effneighb[particle].list[l];

                double /*Anisotropic gather distance */ xi_ij_i =
                    effneighb[particle].distance[l];
                double /*anisotropic-gather kernel */ W_ij_gather =
                    K3D(xi_ij_i) / vol[particle];

                double /*Anisotropic scatter distance */ xi_ij_j;
                descriptor_t    ll;

                /** Find the distance measured from particle j: */
                for (ll = 0; ll < effneighb[j].nn; ll++) {
                    if (effneighb[j].list[ll] == particle)
                        break;
                }
                xi_ij_j = effneighb[j].distance[ll];
                double /*anisotropic-scatter kernel */ W_ij_scatter =
                    K3D(xi_ij_j) / vol[j];

                /** Gather-scatter SPH density is explicitly computed in the following summation: **/
                rho[particle] += .5 * (W_ij_gather + W_ij_scatter) * m[j];
            }   // end for

#ifndef H_ATOM
// #    define _MULTIPHASE_MODEL_2_
#    ifdef _MULTIPHASE_MODEL_2_
#        undef _MULTIPHASE_MODEL_
#        undef _ADIABATIC_
#    endif
#else
#    undef H2_MOLECULE
#endif

            // #define _ADIABATIC_
#ifdef _ADIABATIC_
#    undef _MULTIPHASE_MODEL_
#    undef _MULTIPHASE_MODEL_2_
#endif

#ifdef _MULTIPHASE_MODEL_
            /** the multiphase adiabatic indices are computed here: **/

            /** rho critical is the density at which the fluid starts to smoothly reduce the number 
            of degrees of freedom down to 1/3, which is the approach proposed by Monaghan (1994) **/
            rhocrit = 100;
            double          a = 10,
                zeta = rho[particle] / rhocrit;

            oneoverCv[particle] = isgas[particle] * (2. / (((2 * Cv + 1 / 3.) / (exp(a * (zeta - 1)) + 1)) +
                                                           1. / 3. * (exp(a * (zeta - 1)) - 1) / (exp(a * (zeta - 1)) + 1)));
#elif defined(_MULTIPHASE_MODEL_2_)
            rhocrit = 22.1;
            //             rhocrit = 100;
            double          zeta = rho[particle] / rhocrit;

            //         if (zeta > 5) {
            //             oneoverCv[particle] = 1 / Cv;
            //         } else {
            oneoverCv[particle] = isgas[particle] * (1 / Cv * (1 - exp(-zeta)));
            //         }
#else
            /** otherwise adopt the conventional adiabatic model, either for H atom or for H2 molecule **/
            oneoverCv[particle] = isgas[particle] * (1 / Cv);
#endif

            /** Compute the adiabatic speed of the sound: */
            gamma_index = 1 + oneoverCv[particle] /** gamma_index = cp / cv = 1 + 1/cv = 1 + 0.4 for H_2 **/ ;
            cs[particle] = sqrt(gamma_index * oneoverCv[particle] * u[particle]);

        }   // end if

}   // end-void sph_init...

/** Synopsis:
 * void sph_gradients(index_t i);
 * A procedure to compute mostly SPH quantities which depend on the vector components of the smoothing-kernel gradients.
 * The procedure provides the following div-v dependent quantities for i-particle:
 * -> P_accel[i]: pressure acceleration;
 * -> rhodot[i]: density derivative;
 * -> udot[i]: adiabatic-heating rate. */
void sph_gradients(index_t i)
{
    descriptor_t    l;

    explicit_bzero((void *) P_accel[i], sizeof(double) * DIM);
    mu_max[i] = 0;
    udot_old[i] = udot[i];
    udot[i] = 0;
    /** This is necessary since the variable udot[i] denotes a summation. 
        To make sense it is necessary to save the just read/computed version of udot, 
        from the input file or from the previous integration scheme, in udot_old[i] 
        before cleaning up, and then applying the Marinho's integration scheme
        (Marinho 1997, Ph.D. Thesis, c.f. Sec. 8.4) **/

    /** Similarly to udot and udot_old **/
    rhodot_old[i] = rhodot[i];
    rhodot[i] = 0;

    for (l = 1; l < effneighb[i].nn; l++) {

        index_t         j = effneighb[i].list[l];

        descriptor_t    k,
                        r;

        double          grad_xi_ij_i[DIM];
        double          grad_xi_ij_j[DIM];
        double          x_ij[DIM];

        /** Get the Anisotropic gather distance: */
        double          xi_ij_i = effneighb[i].distance[l];

        /** and its square: */
        double          xi2_ij_i = pow(xi_ij_i, 2);

        /** Anisotropic scatter distance is not necessarily symmetric if switching the gather particle with the
            scatter one due to different local covariance tensors. Thus, search the ll-descriptor that points back to 
            the query i-particle among the effective neighbors of j-particle: */
        descriptor_t    ll;

        for (ll = 0; ll < effneighb[j].nn; ll++) {
            if (effneighb[j].list[ll] == i)
                break;
        }
        /** Now we get the anisotropic scatter distance: */
        double          xi_ij_j = effneighb[j].distance[ll];

        /** and also its square: */
        double          xi2_ij_j = pow(xi_ij_j, 2);

        /** Gather critical direction (namely, xi_ij_i-gradient): i's p.o.v. */
        for (k = 0; k < DIM; k++)
            grad_xi_ij_i[k] = x[i][k] - x[j][k];
        for (r = 0; r < DIM; r++)
            x_ij[r] = dotproduct(eigenvec[i][r], grad_xi_ij_i) / eigenval[i][r];
        for (k = 0; k < DIM; k++) {
            for (grad_xi_ij_i[k] = r = 0; r < DIM; r++) {
                grad_xi_ij_i[k] += x_ij[r] * eigenvec[i][r][k];
            }
            grad_xi_ij_i[k] /= xi_ij_i;
        }

        /** Scatter critical direction (namely, xi_ij_i-gradient): j's p.o.v. */
        for (k = 0; k < DIM; k++)
            grad_xi_ij_j[k] = x[i][k] - x[j][k];
        for (r = 0; r < DIM; r++)
            x_ij[r] = dotproduct(eigenvec[j][r], grad_xi_ij_j) / eigenval[j][r];
        for (k = 0; k < DIM; k++) {
            for (grad_xi_ij_j[k] = r = 0; r < DIM; r++) {
                grad_xi_ij_j[k] += x_ij[r] * eigenvec[j][r][k];
            }
            grad_xi_ij_j[k] /= xi_ij_j;
        }

        /** The kernel gradient: */
        double          grad_W_ij[DIM];

        /** Compute the gather-scatter, Anisotropic simmetrizing kernel-gradient: */
        for (k = 0; k < DIM; k++) {
            grad_W_ij[k] = .5 * (DK3D(xi_ij_i) / vol[i] * grad_xi_ij_i[k] + DK3D(xi_ij_j) / vol[j] * grad_xi_ij_j[k]);
        }

#define __USE_EXPLICIT_THERMAL_ENERGY_EQUATION__
#ifdef __USE_EXPLICIT_THERMAL_ENERGY_EQUATION__
        /** Pressure model works for a diatomic, adiabatic, perfect gas: P = u * rho / Cv.
            So, P / rho^2 = 1 / Cv * u / rho -- I think this is the right approach since ite relates specific thermal energy with
            the mean number of degrees of freedon. However, in the present varying polytropic index, it may conduct to zero pressure,
            which is desired to simulate isothermal processes -- the big problem here is 1 / Cv becoming zero, in the adopted model, for
            low densities */
        double          P_ij = .5 * (oneoverCv[i] * u[i] / rho[i] + oneoverCv[j] * u[j] / rho[j]);
#else
        /** Another approach is considering gamma = 1 + 1/n, where n is the polytropic index, then
            P = P0 * (rho / rho0) ^ gamma_index, so that P/rho² = P0 * (rho / rho0) ^ gamma_index / rho² -- it should be the right approach
            for adiabatic transformation, but varying polytropic index conducts the result to pressure proportional to densities, which
            allows a residual heating for high polytropic indices -- the big problem is presetting P0 **/
        /** P0 = rho0 * Rgas * 10 / meanweight; **/
        double          P00 = rho0 * Rgas * 10 / meanweight / P_unit;
        double          P_ij =
            .5 * P00 * (pow(rho[i] / rho00, 1 + oneoverCv[i]) / pow(rho[i], 2) +
                        pow(rho[j] / rho00, 1 + oneoverCv[j]) / pow(rho[j], 2));
#endif

        /** Artificial-viscosity parameters: */
#define _ARTIFICIAL_VISC_PARM_4
#ifdef _ARTIFICIAL_VISC_PARM_1
        double          alpha = .5;
        double          beta = 2;
#elif defined(_ARTIFICIAL_VISC_PARM_2)
        double          alpha = 1;
        double          beta = 10;
#elif defined(_ARTIFICIAL_VISC_PARM_3)
        double          alpha = 1;
        double          beta = 20;
#elif defined(_ARTIFICIAL_VISC_PARM_4)  /*best until now */
        double          alpha = 5;
        double          beta = 20;
#elif defined(_ARTIFICIAL_VISC_PARM_5)
        double          alpha = 5;
        double          beta = 10;
#elif defined(_ARTIFICIAL_VISC_PARM_6)
        double          alpha = 10;
        double          beta = 40;
#elif defined(_ARTIFICIAL_VISC_PARM_7)
        double          alpha = 1;
        double          beta = 30;
#elif defined(_ARTIFICIAL_VISC_PARM_8)
        double          alpha = 4;
        double          beta = 20;
#elif defined(_ARTIFICIAL_VISC_PARM_9)
        double          alpha = 1;
        double          beta = 2;
#else
        double          alpha = .5;
        double          beta = 1;
#endif
        double          eta2 = .01;
        double          mu_ij;

        /** end artificial viscosity setting **/

        /** create dx_ij[DIM] and dv_ij[DIM] to easy readability: */
        double          dx_ij[DIM],
                        dv_ij[DIM];

        for (k = 0; k < DIM; k++) {
            /** i-to-j relative position vector: */
            dx_ij[k] = x[i][k] - x[j][k];
            /** i-to-j relative velocity vector: */
            dv_ij[k] = v[i][k] - v[j][k];
        }

/** STILL UNDER TEST **/

        /********************************************************************************************************
         * Project: AN UNDER CONSTRUCTION ANISOTROPIC ARTIFICIAL VISCOSITY MODEL FOR ANISOTROPIC SPH SIMULATIONS
         * Author: Eraldo Pereira Marinho, Ph.D.
         *
         * This is my tensor model for artificial viscosity. Firstly I have defined the smoothing tensor as a
         * quadratic form rather than the Shapiro etal model, which is linear in coordinates. In this case,
         * the anisotropic distance, here, is defined as
         *      D(\vec a, \vec b) = \sqrt{\transpose(\vec a - \vec b) \tensor{H}^{-1} (\vec a - \vec b)}.
         * By the way, the smoothing tensor is defined here as proportional to the covariance tensor, namely
         *      \tensor{H} = \xi^2_\mathrm{top} \tensor{\sigma},
         * where \tensor{\Sigma} is the self-regulating covariance tensor for the self-regulating,
         * anisotropic, k-NN cluster; \xi^2_\mathrm{top} is the top distance (Anisotropic measured) from the
         * query particle to its k-NN. Thus, after re-normalization, distances \delta^2 having \tensor{H} as the
         * metric tensor are converted from Anisotropic distance \xi^2 as \delta^2 = \xi^2 / \xi^2_\mathrm{top}.
         * Therefore, \tensor{H} is the metric tensor for which the unit sphere in the H-space correspond to
         * the elliposoid in the Euclidean space with semi-major axis h_1, h_2, h_3, where
         * h^2_j = \xi^2_\mathrm{top}\sigma^2_j, j = 1, 2, 3, and sigma^2_j being the j-th principal component
         * of the self-regulating covariance tensor. Such ellipsoid is named the smoothing ellipsoid.
         *
         *      In order to have the analogous, linear scale, smoothing tensor, I have introduced the
         * square root of the smoothing tensor, namely, \tensor{H}^{-1/2}.
         */
        /** In this case, velocities are left-transformed by the square root of the smoothing tensor, whereas
            the relative positions are right trasformed, yielding the following scalar: */
        for (mu_ij = r = 0; r < DIM; r++) {
            /***********************************************************************************************************************
             * The gather, diagonal form would be:
             * \mu_{ij}=\transpose(\vec v_i-\vec v_j)\frac{\hat{\vec e_r}_i\transpose{\hat{\vec e_r}_i}}{h_{i_r}}(\vec x_i-\vec x_j)
             **********************************************************************************************************************/
            /** The gather-scatter form: */
            mu_ij += .5 * (dotproduct(dv_ij, eigenvec[i][r]) * dotproduct(eigenvec[i][r], dx_ij) / h[i][r]
                           + dotproduct(dv_ij, eigenvec[j][r]) * dotproduct(eigenvec[j][r], dx_ij) / h[j][r]);
        }/** [mu_ij] = [v] */

        /** Particles i and j are approaching if and only if the mu_ij velocity projection is negative: */
        if (mu_ij < 0) {
            /** [mu_ij] = [v] */
            mu_ij = -mu_ij / (.5 * (xi2_ij_i + xi2_ij_j) + eta2);

            mu_max[i] = max(mu_max[i], mu_ij);
            /** midpoint adiabatic speed of sound: cs = sqrt(gamma_index P / rho) = sqrt(gamma_index/cv u) **/
            double          cs_ij = .5 * (cs[i] + cs[j]);

            /** midpoint density **/
            double          rho_ij = .5 * (rho[i] + rho[j]);

            /** Thus, the Monaghan (1992, p. 550) artificial viscosity pressure formula is invariant (either isotropic or
                anisotropic): **/
            P_ij += (alpha * cs_ij + beta * mu_ij) * mu_ij / rho_ij;
            /** Do I have to find an artificial viscosity in the form of the stress tensor? The answer is surely no! **/
        }
        /**
                (1) density derivative;
                (2) pressure acceleration;
                (3) adiabatic-heating rate. **/

        for (k = 0; k < DIM; k++) {
            double          m_grad_W_ij_k = m[j] * grad_W_ij[k],
                m_grad_dv_ij_k = m_grad_W_ij_k * (v[i][k] - v[j][k]);

            /** (1)*/ rhodot[i] += m_grad_dv_ij_k /** Checked! */ ;
            /** (2)*/ P_accel[i][k] -= m_grad_W_ij_k * P_ij /** Checked! */ ;
            /** (3)*/ udot[i] += .5 * m_grad_dv_ij_k * P_ij /** Checked! */ ;

            /** for debug purposes only **/
            if (isnan(udot[i])) {
                fprintf(stderr, "thermal energy rate goes to NaN for particle %d and neighbor %d\n", i, j);
                exit(-5);
            }

        }   // end-for (k = 0; k < DIM; k++)

    }   // end-for (l = 1; l < effneighb[i].nn; l++)

}   // end void sph_gradients(index_t i)

/** performs the cross product, aka vector product, between two 3D vectors a and b yielding the result in the vector c*/
void crossproduct(double const *a, double const *b, double *c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

/******************************************************************************/
/*************** Aperture Criterion is computed hereafter *********************/
/** simulate the line of sight solid angle encompassing the selected octant: */
/** WARNING: it requires urgently a careful revision **/
double p2node_solid_angle(int particle, OCTANT const *p)
{
    double          x1,
                    x2,
                    x3,
                    rx,
                    ry,
                    rz,
                    S1,
                    S2,
                    S3,
                    r,
                    r2,
                    p1,
                    p2,
                    p3,
                    L1,
                    L2,
                    L3;

    /** NOTE:
     *       1. p is a pointer to the the cluster member of a covariance octree partition (node, or somtimes
     *       called cell);
     *       2. x[particle] abstracts the observer or measure particle. Thus, the (rx, ry, rz)-array represents the
     *       particle-to-cluster position vector. 2. x1, x2, x3 are, respectively, the p->x - x[particle] projection along
     *       the e_1, e_2, e_3 principal directions. **/

    /** line of sight particle-to-node vector in explicit form (rx, ry, rz) **/
    rx = p->x[0] - x[particle][0];
    ry = p->x[1] - x[particle][1];
    rz = p->x[2] - x[particle][2];

    /** particle-to-node Euclidean distance **/
    r = sqrt(r2 = rx * rx + ry * ry + rz * rz);

    /** line of sight particle-to-node projection onto the first principal direction **/
    x1 = fabs(p->eigenvec[0][0] * rx + p->eigenvec[0][1] * ry + p->eigenvec[0][2] * rz) / r;
    /** line of sight particle-to-node projection onto the second principal direction **/
    x2 = fabs(p->eigenvec[1][0] * rx + p->eigenvec[1][1] * ry + p->eigenvec[1][2] * rz) / r;
    /** line of sight particle-to-node projection onto the third principal direction **/
    x3 = fabs(p->eigenvec[2][0] * rx + p->eigenvec[2][1] * ry + p->eigenvec[2][2] * rz) / r;

    /** p1 is the principal component **/
    p1 = sqrt(p->eigenval[0]);
    /** sqrt(1 - x1 * x1) is the sine of the formed angle between (rx,ry,rz) and the principal direction p1 **/
    L1 = p1 * p1 * (1 - x1 * x1);
    /** p2 is the second principal component, etc. **/
    p2 = sqrt(p->eigenval[1]);
    L2 = p2 * p2 * (1 - x2 * x2);
    /** p3 is the third principal component, etc. **/
    p3 = sqrt(p->eigenval[2]);
    L3 = p3 * p3 * (1 - x3 * x3);
    /** Thus, L1, L2, L3 are aparent sizes of the 3 principal components respectively **/

    /** S1, S2 and S3 are, respectively, the apparent areas of the reduced ellipses against the line of sight **/
    S1 = 4 * pi * p2 * p3 * x1;
    S2 = 4 * pi * p3 * p1 * x2;
    S3 = 4 * pi * p1 * p2 * x3;
    /** The 4-factor comes from considering the ellipsoid twice the size of the confidence ellipsoid - I've not optimized these expressions **/

    /** This is an estimation of solid angle inscribing the apparent ellipsoid **/
    double          Omega = max(S1, max(S2, S3)) / r2;

    /** L1 / r is a rough estimation of the apparent angle aperture to cover L1 etc. **/
    return /*this is a calibration factor to match BH (1986) scheme */ 5. * max(Omega, max(L1, max(L2, L3)) / r2);
}

// double p2node_solid_angle(int particle, OCTANT const *p)
// {
//     double Omega; double a[DIM], b[DIM], L2, d; double d2, r[DIM];
//     index_t key = p->plist[0];
//     d2 = 0;
//     for (descriptor_t j = 0; j <DIM; j++) {
//         a[j] = b[j] = x[key][j];
//         r[j] = p->x[j] - x[particle][j];
//         d2 += r[j] * r[j];
//     }
//     d = sqrt(d2);
//     for (descriptor_t i = 0; i < p->n; i++) {
//         key = p->plist[i];
//         for (descriptor_t j = 0; j < DIM; j++) {
//             a[j] = min(a[j], x[key][j]);
//             b[j] = max(b[j], x[key][j]);
//         }
//     }
//     Omega = 0;
//     for (descriptor_t j = 0; j <DIM; j++) {
//         L2 = pow(b[j] - a[j], 2) * (1 - pow(r[j] / d, 2));
//         Omega = max(Omega, L2);
//     }
//     Omega = Omega / d2;
//     return Omega;
// }

/*******************************************************************************/
/************ Covariance-octree based gravity acceleration computation *********/

double          theta2;

double          epsilon,
                epsilon2;

double          Phi[NMAX];

/** void find_well_distant_nodes(int particle, const OCTANT * cluster); is a recursive method of finding
 *  sigma-octree nodes that are distant enough to have their projected solid angle encompassed by theta².
 *  The method is efficient but it is a very irregular program to parallel using omp tasks. Several attempts
 *  were made but either the speedup were negative or the calculation became inaccurate. I am still studying
 *  a way of parallelizing this code. **/

// #define _HK89_SOFTENING_
/** I have no idea on how to implement quadrupole correction for HK89 softening once the splined kernel is of analytical class C²,
 *  since the kernel's third derivative is discontinuous. Hence, the HK89 idea have to be replaced with a quartic B-spline. **/
#ifdef _HK89_SOFTENING_

void find_well_distant_nodes(int particle, const OCTANT * cluster)
{
    for (descriptor_t octant = 0; octant < NCHILDREN; octant++) {
        OCTANT         *child = cluster->child[octant];

        /** refuse empty child octant */
        if (child == NULL)
            continue;

        /** recursively descend each child to find either a well-distant cluster or a singleton one */
        if (( /** is the octant-cluster a proper singleton[1]? */ child->n == 1 && child->plist[0] != particle) ||
            ( /** or is the octant-cluster well distant? */ p2node_solid_angle(particle, child) < theta2)) {
            /** [1]: NOTE: by proper singleton we mean child->plist[0] != particle, say, not gravity singularity **/

            double          rx,
                            ry,
                            rz,
                            r,
                            r2;

            rx = (child->x[0] - x[particle][0]);
            ry = (child->x[1] - x[particle][1]);
            rz = (child->x[2] - x[particle][2]);

            /** take the square euclidean distance */
            r2 = pow(rx, 2) + pow(ry, 2) + pow(rz, 2);
            /** and then the conventional one: **/
            r = sqrt(r2);

            /** sum up gravity potential (positive definite) contribution from the octant **/
            Phi[particle] += child->mass * fHK89(r);

            /** sum up the gravity acceleration contribution from the child octant */
            double          HK89_g = child->mass * gHK89(r);

            g[particle][0] += HK89_g * rx;
            g[particle][1] += HK89_g * ry;
            g[particle][2] += HK89_g * rz;

        } else {

            /** if the cluster child is not distant enough from the particle, proceed the descent task */
            find_well_distant_nodes(particle, child);

        }   // end if

    }   // end for

}   // end void find_well_distant_nodes

#else

/** Arseth model of softened potential is adopted here **/
void find_well_distant_nodes(int particle, const OCTANT * cluster)
{
    for (descriptor_t octant = 0; octant < NCHILDREN; octant++) {
        OCTANT         *child = cluster->child[octant];

        /** refuse empty child octant -- it happens! **/
        if (child == NULL)
            continue;

        /** recursively descend each child to find either a well-distant cluster or a singleton one */
        if (( /** is the octant-cluster a proper singleton¹? */ child->n == 1 && child->plist[0] != particle) ||
            ( /** or is the octant-cluster well distant? */ p2node_solid_angle(particle, child) < theta2)) {
            /** ¹: NOTE: by proper singleton we mean child->plist[0] != particle, say, not gravity singularity **/

            double          r[DIM],
                            r1,
                            r2,
                            r3;

            for (descriptor_t j = 0; j < DIM; j++) {
                r[j] = child->x[j] - x[particle][j];
            }

            /** take the square softened distance */
            r2 = pow(r[0], 2) + pow(r[1], 2) + pow(r[2], 2) + epsilon2;
            /** and then the softened distance: **/
            r1 = sqrt(r2);
            /** the cubic softened distance: **/
            r3 = r1 * r2;

            /** sum up gravity potential (positive definite) contribution from the octant **/
            Phi[particle] += child->mass / r1;

            /** sum up the gravity acceleration contribution from the child octant */
            double          src = child->mass / r3;

            for (descriptor_t j = 0; j < DIM; j++) {
                g[particle][j] += src * r[j];
            }

#    define _QDP_CONTRIB_
#    ifdef _QDP_CONTRIB_

            if (child->n > 2) {

                /** fifth power of the softened distance **/
                double          r5 = r2 * r3;

                /*** softened quadrupole components ***/
                double          /** product of inertia **/PP[DIM][DIM],
                                /** polar momentum of inertia **/QQ = 0;

                for (descriptor_t j = 0; j < DIM; j++)
                    explicit_bzero(PP[j], sizeof(**PP) * DIM);

                for (descriptor_t l = 0; l < DIM; l++) {
                    QQ += child->mass * child->eigenval[l];
                    for (descriptor_t j = 0; j < DIM; j++) {
                        /*** diagonal contribution: ***/
                        PP[j][j] += child->eigenvec[l][j] * child->eigenval[l] * child->eigenvec[l][j];
                        for (descriptor_t k = j + 1; k < DIM; k++) {
                            /*** upper triangle summation: ***/
                            PP[j][k] += child->eigenvec[l][j] * child->eigenval[l] * child->eigenvec[l][k];
                        }
                    }
                }

                for (descriptor_t j = 0; j < DIM; j++) {
                    /*** diagonal contribution: ***/
                    PP[j][j] *= child->mass;
                    for (descriptor_t k = j + 1; k < DIM; k++) {
                        /*** upper triangle contribution: ***/
                        PP[j][k] *= child->mass;
                        /*** symmetric closure: do the lower triangle: ***/
                        PP[k][j] = PP[j][k];
                    }
                }

                double          tensorcontrib = 0,
                    veccontrib[DIM];

                for (descriptor_t j = 0; j < DIM; j++) {
                    veccontrib[j] = 0;
                    for (descriptor_t k = 0; k < DIM; k++) {
                        tensorcontrib += r[j] * PP[j][k] * r[k];
                        veccontrib[j] += PP[j][k] * r[k];
                    }
                }

                for (descriptor_t j = 0; j < DIM; j++) {
                    g[particle][j] += 1.5 * ((5 * tensorcontrib / r2 - QQ) * r[j] - 2 * veccontrib[j]) / r5;
                }

                Phi[particle] += .5 * (3 * tensorcontrib / r2 - QQ) / r3;

            }   // end if (child->n > 2)

#    endif

        } else {

            /** if the cluster child is not distant enough from the particle, proceed the descent task in the search for smaller partitions */
            find_well_distant_nodes(particle, child);

        }   // end if

    }   // end for

}   // end void find_well_distant_nodes
#endif

void sigmatreegravity(index_t particle)
{
    /** Clean-up the gravity components for the current particle: */
    explicit_bzero(g[particle], sizeof(*g[particle]) * DIM);
    Phi[particle] = 0;
    /** Compute the gravity acceleration for the current particle: */
    find_well_distant_nodes(particle, sigma_octree_root);
}

/****************************** Time integration (the present modified leapfrog is an ad hoc scheme) ****************************/

#define MAXTIMEDEPTH 12
double          DT /* master time-step (user defined or default set-up ) */ ;
double          dt[NMAX],
                dt_old[NMAX];

/** Hilbert-Courant stability conditions: */
size_t          dtfault = 0;
double Hilbert_Courant(index_t particle)
{
    //fprintf(stderr, "entering Hilbert_Courant\n");
    double          dtmin = DT / (1 << MAXTIMEDEPTH);

    //     double          dtmax;
    double          dt_estimate = DT,
        veff = sqrt(dotproduct(v[particle], v[particle])),
        geff;

    // #define COURANT_FACTOR .0150625
    // #define COURANT_FACTOR .03125
#define COURANT_FACTOR .0625
    // #define COURANT_FACTOR .125

    /** gravity is ubiquitous: **/
    sigmatreegravity(particle);
    /** default: softening-length based times-cales: **/
    geff = sqrt(dotproduct(g[particle], g[particle]));
    // dt_estimate = min(dt_estimate, COURANT_FACTOR * 1 * sqrt(epsilon / geff));
    // dt_estimate = min(dt_estimate, COURANT_FACTOR * 1 * epsilon / veff);
    dt_estimate = min(dt_estimate, COURANT_FACTOR * sqrt(psize[particle] / geff));
    dt_estimate = min(dt_estimate, COURANT_FACTOR * .5 * psize[particle] / veff);

    /** gas (SPH) specific time scales: **/
    if (isgas[particle]) {
        sph_gradients(particle);
        /** artificial-viscosity time scale: **/
        dt_estimate = min(dt_estimate, COURANT_FACTOR * h[particle][2] / max(cs[particle], max(mu_max[particle], veff)));

        /** Jeans time scale: **/
        dt_estimate = min(dt_estimate, COURANT_FACTOR / sqrt(4 * pi /*G*/ * rho[particle]));
        dt_estimate = min(dt_estimate, COURANT_FACTOR * h[particle][2] / geff);

        /** This is the critical part of the integration scheme: energies and densities **/
        dt_estimate = min(dt_estimate, COURANT_FACTOR * .125 * u[particle] / fabs(udot[particle]) + 1);
        dt_estimate = min(dt_estimate, COURANT_FACTOR * .125 * rho[particle] / fabs(rhodot[particle]) + 1);
    }

                                        /** count time-faults: **/
    if (dt_estimate <= dtmin) {
        dtfault++;
        return dtmin;
    }
    /***    
    double          dt_cut = min(dt_estimate, DT);

    for (dtmax = 1; dtmax > dt_cut; dtmax *= .5);

    return dtmax;*/

    return dt_estimate;
}

index_t         timebin[MAXTIMEDEPTH + 1][NMAX];
size_t          timebin_population[MAXTIMEDEPTH + 1];
size_t          max_time_depth = 0;

void timebin_share(size_t n)
{
    explicit_bzero(timebin_population, (MAXTIMEDEPTH + 1) * sizeof(size_t));

#pragma omp parallel for
    for (index_t i = 0; i < n; i++) {
        if (deleted[i] == 0) {
#pragma omp atomic write
            dt[i] = Hilbert_Courant(i);
        }
    }   // end omp parallel for

#pragma omp parallel for reduction(max:max_time_depth)
    for (index_t i = 0; i < n; i++) {
        if (deleted[i])
            continue;

        double          dtaux = DT;
        size_t          time_depth;

        for (time_depth = 0; dtaux > dt[i]; time_depth++) {
            dtaux *= .5;
        }   // end for

        max_time_depth = max(time_depth, max_time_depth);

#pragma omp critical
        {
            timebin[time_depth][timebin_population[time_depth]] = i;
            timebin_population[time_depth]++;
        }   // end critical

    }   // end omp parallel for

#define VERBOSE_TIME_BINS
#ifdef VERBOSE_TIME_BINS
    size_t          ntot = 0;

    fprintf(stderr, "\ntime-bins:\n");
    for (descriptor_t i = 0; i <= max_time_depth; i++) {
        fprintf(stderr, "[%i] -> %lu\n", i, timebin_population[i]);
        ntot += timebin_population[i];
    }
    fprintf(stderr, "chksum = %lu\n", ntot);
    fprintf(stderr, "time faults = %lu\n\n", dtfault);
#endif
}

void binary_leapfrog(descriptor_t time_depth) // this is what has been adopted in the last three years
{
    if (time_depth < max_time_depth)
        binary_leapfrog(time_depth + 1);

    integrate(time_depth);

    if (time_depth < max_time_depth)
        binary_leapfrog(time_depth + 1);
}

// void binary_leapfrog(descriptor_t time_depth) // this is the scheme of Marinho (1997, Ph.D. Thesis) and Marinho-Lépine (2000)
// {
//     integrate(time_depth);
//
//     if (time_depth < max_time_depth)
//         binary_leapfrog(time_depth + 1);
//
//     if (time_depth < max_time_depth)
//         binary_leapfrog(time_depth + 1);
// }

double          udot_old[NMAX],
                rhodot_old[NMAX],
                v0[NMAX][DIM],
                u_0[NMAX],
                rho_0[NMAX],
                u_sav[NMAX],
                dtau_m[NMAX],
                dtau_d[NMAX];

extern flag_t   already_initiated_find_selfregulating_kNN_cluster;

// #define _USE_T_CUT_
#ifdef _USE_T_CUT_
extern double   T_mean;
#endif

size_t          nterm = 0;

void integrate(descriptor_t time_depth)
{

    /** n is the time-bin population: **/
    size_t          n = timebin_population[time_depth];

    if (n < 1)
        return;

    // #define VERBOSE_TIME_DEPTH
#ifdef VERBOSE_TIME_DEPTH
    double          dtt = DT / (1 << time_depth);

    fprintf(stderr, "time_depth = %d; dt = %lg = DT / %d for %lu particle", time_depth, dtt, (1 << time_depth), n);
    if (n > 1) {
        fprintf(stderr, "s");
    }
    fprintf(stderr, "\n");
#endif

#pragma omp parallel for
    for (descriptor_t i = 0; i < n; i++) {
        index_t         particle = timebin[time_depth][i];

        if (deleted[particle])
            continue;

        if (dt_old[particle] == 0) {
            dt_old[particle] = dt[particle];
        }

        dtau_m[particle] = .5 * (dt[particle] + dt_old[particle]);
        dtau_d[particle] = .5 * (dt[particle] - dt_old[particle]);

        //fprintf(stderr, "dt[%d]=%lg\n", particle, dt_old[particle]);

    }   // end omp parallel for

    /** predict positions for time-level n+1/2 **/
#pragma omp parallel for
    for (descriptor_t i = 0; i < n; i++) {
        index_t         particle = timebin[time_depth][i];

        if (deleted[particle])
            continue;

        /** Gravity acceleration has been computed for time-level n-1/2 **/
        if (isgas[particle]) {
            for (descriptor_t j = 0; j < DIM; j++) {
                x[particle][j] +=
                    (v[particle][j] + .5 * (g[particle][j] + P_accel[particle][j]) * dtau_d[particle]) * dtau_m[particle];
            }
        } else {
            for (descriptor_t j = 0; j < DIM; j++) {
                x[particle][j] += (v[particle][j] + .5 * g[particle][j] * dtau_d[particle]) * dtau_m[particle];
            }
        }
    }   // end omp parallel for

    /** Gravity acceleration is now predicted to time-level n+1/2 **/
#pragma omp parallel for
    for (descriptor_t i = 0; i < n; i++) {
        index_t         particle = timebin[time_depth][i];

        if (deleted[particle] == 0) {
            sigmatreegravity(particle);
        }

    }   // end omp parallel for

#pragma omp parallel for
    for (descriptor_t i = 0; i < n; i++) {
        index_t         particle = timebin[time_depth][i];

        if (deleted[particle])
            continue;

        if (isgas[particle]) {
#ifdef _USE_T_CUT_
/** Using cut-off temperature in the simulation is a rudimentary way of simulating cooling effects at high densities **/
#    define T_cut 1000
                    /*Kelvin */
            /*** inserted in May 26 2022 ***/
            /*** inserted in June 1st 2022 ***/
            double          T_particle = meanweight * oneoverCv[particle] * u_unit * u[particle] / Rgas;

            if (T_particle > T_cut) {
                //                 u[particle] = T_cut / (meanweight * oneoverCv[i] * u_unit / Rgas);
                u[particle] = T_cut / T_particle * u[particle];
                T_particle = T_cut;
            }
            /** inserted in June the First, 2022: **/
            T_mean += T_particle;
            nterm++;
#endif
            /***/
            rho_0[particle] = rho[particle];
            u_0[particle] = u[particle];
            u_sav[particle] = u[particle];
        }

        memcpy(v0[particle], v[particle], sizeof(double) * DIM);

    }   // end omp parallel for

    counter_t       iter = 0;
    flag_t          not_converged;

    /** convergence check against what? 
     *  - Answer: specific thermal energies **/
    for (not_converged = 1; not_converged;) {
        iter++;

#pragma omp parallel for
        for (descriptor_t i = 0; i < n; i++) {
            index_t         particle = timebin[time_depth][i];

            if (deleted[particle])
                continue;

            if (isgas[particle] == 0)
                continue;

            double          du;

#define MAXCOUNT    127
            counter_t       count = 0;

            /** entering this loop means the timestep is not sufficiently small **/
            while ((du = .5 * dtau_m[particle] * (udot[particle] + udot_old[particle]) / u_0[particle]) <= -1) {
                if (count++ > MAXCOUNT)
                    break;
                udot[particle] *= .5;
                udot_old[particle] *= .5;
            }
            if (du <= -1 && count > MAXCOUNT) {
                fprintf(stderr, "udot reduction reached MAXCOUNT=%d! Try smaller timesteps\n", MAXCOUNT);
                //                 fprintf(stderr, "deleting particle %d at depth %d\n", particle, time_depth);
                //                 deleted[particle] = 1;
                //                 continue;
            }
            u[particle] = u_0[particle] * (1 + du);
            /*** inserted in May 26 2022 ***/
            /*** inserted in June 1st 2022 ***/
#ifdef _USE_T_CUT_
            double          T_particle = meanweight * oneoverCv[particle] * u_unit * u[particle] / Rgas;

            if (T_particle > T_cut) {
                //                 u[particle] = T_cut / (meanweight * oneoverCv[i] * u_unit / Rgas);
                u[particle] = T_cut / T_particle * u[particle];
                T_particle = T_cut;
            }
            /** inserted in June the First, 2022: **/
            T_mean += T_particle;
            nterm++;
#endif
            /***/

            if (u[particle] <= 0) {
                fprintf(stderr, "Cannot deal with non-positive thermal energy u[%d] = %lg. ", particle, u[particle]);
                u[particle] = u_0[particle];
                //                 fprintf(stderr, "deleting particle %d at depth %d\n", particle, time_depth);
                //                 deleted[particle] = 1;
                continue;
            }
        }   // end omp parallel for

#pragma omp parallel for
        for (descriptor_t i = 0; i < n; i++) {
            index_t         particle = timebin[time_depth][i];

            if (deleted[particle])
                continue;

            if (isgas[particle] == 0)
                continue;

            double          drho;

            counter_t       count = 0;

            while ((drho = .5 * dtau_m[particle] * (rhodot[particle] + rhodot_old[particle]) / rho_0[particle]) <= -1) {
                if (count++ > MAXCOUNT)
                    break;
                rhodot[particle] *= .5;
                rhodot_old[particle] *= .5;
            }
            if (drho <= -1 && count > MAXCOUNT) {
                fprintf(stderr, "rhodot reduction reached MAXCOUNT=%d! Try smaller timesteps\n", MAXCOUNT);
                //                 fprintf(stderr, "deleting particle %d at depth %d\n", particle, time_depth);
                //                 deleted[particle] = 1;
                //                 continue;
            }
            rho[particle] = rho_0[particle] * (1 + drho);

            if (rho[particle] <= 0) {
                fprintf(stderr, "Cannot deal with non-positive density rho[%d] = %lg. ", particle, rho[particle]);
                rho[particle] = rho_0[particle];
                //                 fprintf(stderr, "deleting particle %d at depth %d\n", particle, time_depth);
                //                 deleted[particle] = 1;
                //                 continue;
            }

#ifdef _MULTIPHASE_MODEL_
            double          a = 10,
                zeta = rho[particle] / rhocrit;

            oneoverCv[particle] =
                isgas[particle] * (2. / (((2 * Cv + 1 / 3.) / (exp(a * (zeta - 1)) + 1)) +
                                         1. / 3. * (exp(a * (zeta - 1)) - 1) / (exp(a * (zeta - 1)) + 1)));
#elif defined(_MULTIPHASE_MODEL_2_)
            double          zeta = rho[particle] / rhocrit;

            //             if (zeta > 5) {
            //                 oneoverCv[particle] = 1 / Cv;
            //             } else {
            oneoverCv[particle] = isgas[particle] * (1 / Cv * (1 - exp(-zeta)));
            //             }

#else
            oneoverCv[particle] = isgas[particle] * (1 / Cv);
#endif
            gamma_index = 1 + oneoverCv[particle] /** gamma_index = cp / cv = 1 + 1/cv = 1 + 0.4 in H2 model **/ ;
            cs[particle] = sqrt(gamma_index * oneoverCv[particle] * u[particle]);
        }   // end parallel for

        /** get the midpoint velocities, namely, predicted to time-level n+1/2 */
        if (iter > 1) {

#pragma omp parallel for
            for (descriptor_t i = 0; i < n; i++) {
                index_t         particle = timebin[time_depth][i];

                if (deleted[particle])
                    continue;

                if (isgas[particle] == 0)
                    continue;

                for (descriptor_t j = 0; j < DIM; j++) {
                    v[particle][j] = .5 * (v[particle][j] + v0[particle][j]);
                }
            }   // end omp parallel for

        }
        /** compute SPH gradients, yielding pressure acceleration etc. for midpoint velocities, n+1/2 */
#pragma omp parallel for
        for (descriptor_t i = 0; i < n; i++) {
            index_t         particle = timebin[time_depth][i];

            if (deleted[particle])
                continue;

            if (isgas[particle] == 0)
                continue;

            sph_gradients(particle);
        }   //end omp parallel for

        /** integrate velocities covering the enire time interval n -> n+1 */
#pragma omp parallel for
        for (descriptor_t i = 0; i < n; i++) {
            index_t         particle = timebin[time_depth][i];

            if (deleted[particle])
                continue;

            if (isgas[particle]) {
                for (descriptor_t j = 0; j < DIM; j++) {
                    v[particle][j] = v0[particle][j] + (g[particle][j] + P_accel[particle][j]) * dt[particle];
                }
            } else {
                for (descriptor_t j = 0; j < DIM; j++) {
                    v[particle][j] = v0[particle][j] + g[particle][j] * dt[particle];
                }
            }

            double          v2;

            if ((v2 = dotproduct(v[particle], v[particle])) >= c_norm * c_norm) {
                deleted[particle] = 1;
                fprintf(stderr, "particle %d reached light speed v²/c² = %lg, deleting it\n", particle, v2 / (c_norm * c_norm));
            }
        }   //end parallel for

        /** convergence check in function of thermal specific energies: */
        double          diff = 0,
            U = 0;

#pragma omp parallel for
        for (descriptor_t i = 0; i < n; i++) {
            index_t         particle = timebin[time_depth][i];

            if (deleted[particle])
                continue;

            if (isgas[particle] == 0)
                continue;

            U += m[particle] * u[particle];

            diff += m[particle] * fabs(u[particle] - u_sav[particle]);

        }   //end parallel for

#define CONV_TOL 5.0e-12
        not_converged = (diff > CONV_TOL * U);
        if (not not_converged) {
            break;
        }
#define MAXITER 63
        if (iter > MAXITER) {
            fprintf(stderr,
                    "Specif thermal energies not converging, Udiff/U = %lg/%lg = %lg. Breaking iteration %i\n", diff,
                    U, diff / U, iter);
            break;
        }

#pragma omp parallel for
        for (descriptor_t i = 0; i < n; i++) {
            index_t         particle = timebin[time_depth][i];

            if (deleted[particle])
                continue;

            if (isgas[particle] == 0)
                continue;

            u_sav[particle] = u[particle];
        }   //end parallel for

    }   //end for (not_converged = 1; not_converged;)

    //fprintf(stderr, "t = %lg -> %d iterations\n", t / dt, iter);
#pragma omp parallel for
    for (descriptor_t i = 0; i < n; i++) {
        index_t         particle = timebin[time_depth][i];

        if (deleted[particle])
            continue;

        dt_old[particle] = dt[particle];
    }   //end parallel for

}

double
                M = 0,
    E = 0,
    W = 0,
    T = 0,
    U = 0,
    X[DIM],
    V[DIM],
    L[DIM]
#ifdef _USE_T_CUT_
    ,
    T_mean = 0
#endif
    ;

void integral_quantities(size_t n)
{
    //     fprintf(stderr, "\nentering integral_quantities\n");

    M = E = W = T = U = 0;

    explicit_bzero(X, DIM * sizeof(*X));
    explicit_bzero(V, DIM * sizeof(*V));
    explicit_bzero(L, DIM * sizeof(*L));

#pragma omp parallel shared(M, W, T, U, X, V, L) num_threads(numthreads)
    {
        double /** thread specific (cache) variables: **/ thread_M = 0,
            thread_T = 0,
            thread_U = 0,
            thread_W = 0,
            thread_L[DIM],
            thread_X[DIM],
            thread_V[DIM];

        /** temporary variable to store x cross v: **/
        double          xcrossv[DIM];

        explicit_bzero(thread_L, sizeof(*L) * DIM);
        explicit_bzero(thread_X, sizeof(*X) * DIM);
        explicit_bzero(thread_V, sizeof(*V) * DIM);

#pragma omp for
        for (index_t i = 0; i < n; i++) {
            if (deleted[i])
                continue;

            thread_M += m[i];

            thread_T += .5 * m[i] * dotproduct(v[i], v[i]);

            thread_W += .5 * m[i] * Phi[i];

            if (isgas[i]) {
                thread_U += m[i] * u[i];
            }

            crossproduct(x[i], v[i], xcrossv);

            for (descriptor_t j = 0; j < DIM; j++) {
                thread_L[j] += m[i] * xcrossv[j];
                thread_X[j] += m[i] * x[i][j];
                thread_V[j] += m[i] * v[i][j];
            }

        }   //end omp for

        /*** gathering all of the threads contributions to the shared variables: ***/

#pragma omp atomic
        M += thread_M;
#pragma omp atomic
        T += thread_T;
#pragma omp atomic
        U += thread_U;
#pragma omp atomic
        W += thread_W;
        for (descriptor_t j = 0; j < DIM; j++) {
#pragma omp atomic
            X[j] += thread_X[j];
#pragma omp atomic
            V[j] += thread_V[j];
#pragma omp atomic
            L[j] += thread_L[j];
        }

    }   //end omp parallel

    for (descriptor_t j = 0; j < DIM; j++) {
        X[j] /= M;
        V[j] /= M;
    }

    W *= -1;

    E = U + T + W;

#ifdef _USE_T_CUT_
    T_mean /= nterm;
#endif

    fprintf(stderr, "Total mass = %lg\n", M);
    fprintf(stderr, "Center of mass = ( %lg, %lg, %lg )\n", X[0], X[1], X[2]);
    fprintf(stderr, "Mean velocity = ( %lg, %lg, %lg )\n", V[0], V[1], V[2]);
    fprintf(stderr, "Angular momentum = ( %lg, %lg, %lg )\n", L[0], L[1], L[2]);
    fprintf(stderr, "Thermal energy = %lg\n", U);
#ifdef _USE_T_CUT_
    fprintf(stderr, "Temperature = %lg Kelvin\n", T_mean);  //T_unit * (U / M));
#endif
    fprintf(stderr, "Kinetic energy = %lg\n", T);
    fprintf(stderr, "Gravitational energy = %lg\n", W);
    fprintf(stderr, "Total energy = %lg\n", E);
}
