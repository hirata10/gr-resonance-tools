#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "CKerr.h"

/* Computing the J_dot due to self-force from arXiv:1905.00030v2 Eq. (12) */

/* LEGACY CODE (outdated; kept for reference): J_dot_selfforce, the serial self-force.
 * The pipeline does not use it: steps 4 and 7 call J_dot_selfforce_openmp (below). Differences from that routine:
 *   - no mass ratio: it returns the rate for a unit mass ratio, not d(J/mu_body)/dt in physical time;
 *   - the radial and polar grid is fixed at 28 x 14;
 *   - serial, no orbit cache, and no use of the (n,k,m) -> (-n,-k,-m) symmetry of the mode sum.
 * Still called by other legacy code only: J2Jdot and J2Jdot_component (J2J_dot.c), omega_dot (Gamma.c),
 * the Python wrapper j_dot_selfforce, and the serial debug tool J2J_DOT_SF (calling.c, -DIS_J_DOT_SF -DSINGLEVALUE).
 * For tests use the OpenMP tool J2J_DOT_SF_OMP. */
void J_dot_selfforce(int nl, int nmax, int kmax, int mmax, double apo,
                     double rp, double radius_outer, double I, double M,
                     double astar, double *J_dot_sf) {
  int i, i_n, i_k, i_m, il;
  double EQL[3], J[3], Minv[9], Omega[3], info[6], info_outer[6], xuorig[6],
      xuorig_outer[6];
  // int n_res_inner = 1, k_res_inner = 2, m_res_inner = -2;
  double Z_out_square, Z_down_square;
  double omega_nkm, omegagw;
  double rH;
  double J_dot_r_old, J_dot_theta_old, J_dot_phi_old;
  double epsilon, lambda, numer, C2, P, alphankm;
  double *C0, *E0, *E1;

  /* For CIRCULAR orbits */
  if (apo == 0 && rp == 0 && I == 0) {
    /* Orbital frequencies for circular equitorial orbits */
    CKerr_FindEQL_IRCirc(0, radius_outer, EQL, M, astar);
    CKerr_EQL2J(EQL, J, M, astar, NULL);
    CKerr_Minverse(J, Minv, M, astar);
    CKerr_Minv2Omega(Minv, Omega);
    // printf("Minv for circular equitorial case = \n %lg %lg %lg \n", Minv[0],
    // Minv[1], Minv[2]); printf(" %lg %lg %lg \n", Minv[3], Minv[4], Minv[5]);
    // printf(" %lg %lg %lg \n", Minv[6], Minv[7], Minv[8]);
    // printf("Js are: %lg %lg %lg \n", J[0], J[1], J[2]);
    // printf("Omega for circular equatorial is: %lg %lg %lg \n", Omega[0],
    // Omega[1], Omega[2]);

    /* Torus Origin for circular equatorial orbits */
    CKerr_TorusOrigin(J, xuorig, M, astar);
  }

  /* For GENERIC orbits */
  else {
    /* Routine to find location of resonant apocenter */
    // guess1 = rp + 3.;
    // guess2 = radius_outer - 3.;
    // apo = find_resonance_apo(n_res_inner, k_res_inner, m_res_inner,
    // radius_outer, guess1, guess2, rp, I, astar, M);

    /* Get orbital frequencies for inner generic orbit */
    ra_rp_I2EQL(apo, EQL, rp, I, astar, M);
    CKerr_EQL2J(EQL, J, M, astar, NULL);
    CKerr_Minverse(J, Minv, M, astar);
    CKerr_Minv2Omega(Minv, Omega);
    // printf("Minv for generic case = \n %lg %lg %lg \n", Minv[0], Minv[1],
    // Minv[2]); printf(" %lg %lg %lg \n", Minv[3], Minv[4], Minv[5]); printf("
    // %lg %lg %lg \n", Minv[6], Minv[7], Minv[8]); printf("Js are: %lg %lg %lg
    // \n", J[0], J[1], J[2]); printf("Omega for generic is: %lg %lg %lg \n",
    // Omega[0], Omega[1], Omega[2]);

    /* Torus origin for inner body orbit */
    CKerr_TorusOrigin(J, xuorig, M, astar);
    // for(i=0;i<6;i++) printf("xuorig[%d] = %lg\n", i, xuorig[i]);
  }

  /* This takes the M_inv, the origin in phase space, mass and spin of BH,
  {n,k,m}, number of steps in r and theta, number of modes (nl), and the angular
  frequency of the perturber orbit Returns the Z-coeff. and energy eigenvalues
  E[nl]. This is also computed at the specific resonance condition (1, 2, -2).
*/
  E0 = (double *)malloc((size_t)(nl * 10 * sizeof(double)));
  E1 = E0 + nl;
  C0 = E1 + nl;

  // printf("IM HERE 3 \n");

  rH = M + sqrt(M * M - astar * M * astar * M);
  epsilon = sqrt(M * M - astar * M * astar * M) / (4 * M * rH);

  J_dot_sf[0] = 0;
  J_dot_sf[1] = 0;
  J_dot_sf[2] = 0;
  // printf("IM HERE 4 \n");
  /* Compute the J_dots for the inner orbit */
  for (i_n = -nmax; i_n <= nmax; i_n++) {
    for (i_k = -kmax; i_k <= kmax; i_k++) {
      for (i_m = -mmax; i_m <= mmax; i_m++) {
        if (i_n == 0 && i_k == 0 && i_m == 0) {
          continue;
        }
        if (astar == 0 && i_n == 0 && i_k == -i_m) // Schwarzschild case a == 0
          continue;
        CKerr_RadialFunc(Minv, xuorig, M, astar, i_n, i_k, i_m, 28, 14, nl, C0,
                         &omegagw, E0);
        for (il = 0; (il) < nl; il++) {

          omega_nkm = i_n * Omega[0] + i_k * Omega[1] + i_m * Omega[2];
          lambda = E0[il] - 2 * astar * M * i_m * omega_nkm +
                   astar * M * astar * M * omega_nkm * omega_nkm - 2;
          P = omega_nkm - i_m * astar * M / (2 * M * rH);
          numer = 256 * pow(2 * M * rH, 5) * P *
                  (P * P + 4 * epsilon * epsilon) *
                  (P * P + 16 * epsilon * epsilon) * omega_nkm * omega_nkm *
                  omega_nkm;
          C2 = ((lambda + 2) * (lambda + 2) + 4 * i_m * astar * omega_nkm -
                4 * astar * M * astar * M * omega_nkm * omega_nkm) *
                   (lambda * lambda + 36 * i_m * astar * M * omega_nkm -
                    36 * astar * M * astar * M * omega_nkm * omega_nkm) +
               (2 * lambda + 3) *
                   (96 * astar * M * astar * M * omega_nkm * omega_nkm -
                    48 * i_m * astar * M * omega_nkm) +
               144 * omega_nkm * omega_nkm * (M * M - astar * astar * M * M);
          alphankm = numer / C2;

          Z_down_square =
              C0[4 * il] * C0[4 * il] + C0[4 * il + 1] * C0[4 * il + 1];
          Z_out_square =
              C0[4 * il + 2] * C0[4 * il + 2] + C0[4 * il + 3] * C0[4 * il + 3];
#if 0
					printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n, i_k, i_m, il, omegagw, omega_nkm, lambda, alphankm, Z_down_square, Z_out_square);
					fflush(stdout);
					// printf("%i \t %i \t %i \t %i \t %lg \t %lg \n", i_n, i_k, i_m, il, Z_down_square, Z_out_square);
#endif
          // printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t
          // %lg\n", i_n, i_k, i_m, il, C0[4*il], C0[4*il+1], C0[4*il+2],
          // C0[4*il+3], Z_down_square, Z_out_square);

          J_dot_sf[0] += -i_n * (Z_out_square + alphankm * Z_down_square) /
                         (2. * omega_nkm * omega_nkm * omega_nkm);
          J_dot_sf[1] += -i_k * (Z_out_square + alphankm * Z_down_square) /
                         (2. * omega_nkm * omega_nkm * omega_nkm);
          J_dot_sf[2] += -i_m * (Z_out_square + alphankm * Z_down_square) /
                         (2. * omega_nkm * omega_nkm * omega_nkm);
        }
      }
    }
  }
  free((char *)E0);
}

void J_dot_selfforce_openmp(int nl, int nmax, int kmax, int mmax, int Nr, int Nt, double apo,
                            double rp, double I, double mu_body, double M,
                            double astar, double *J_dot_sf) {

  double EQL[3], J[3], Minv[9], Omega[3], info[6], info_outer[6], xuorig[6];
  double radial_time_global = 0.0;
  long radial_calls_global = 0;

  double J_dot_sf_global[3] = {0., 0., 0.};

    // Body apocenter, pericenter, inclination -> EQL and Minv
  ra_rp_I2EQL(apo, EQL, rp, I, astar, M);
  for (int i = 0; i < 3; i++) {
    	if (!isfinite(EQL[i])) {
        printf("J_dot.c: No values of EQL found, exiting...\n");
        exit(1);
		}
	}

  double anc[3];
  CKerr_EQL2J(EQL, J, M, astar, anc);
  for (int i = 0; i < 3; i++) {
    	if (!isfinite(J[i])) {
        printf("J_dot.c: No values of J found, exiting...\n");
        exit(1);
		}
	}

  CKerr_Minverse(J, Minv, M, astar);
  CKerr_TorusOrigin(J, xuorig, M, astar);

  /* Compute outer Kerr horizon and function needed for scattering data
   * (epsilon) */
  double rH = M + sqrt(M * M - astar * M * astar * M);
  double epsilon = sqrt(M * M - astar * M * astar * M) / (4 * M * rH);

#ifndef JDOT_SF_ORIGINAL
  /* Orbit cache (radialfunc_cached.c): torus samples and Y^{-2}_{lm} once per J_dot, l up to max(2,mmax)+nl-1.
   * Build with -DJDOT_SF_ORIGINAL (make JDOT_SF=original) to call CKerr_RadialFunc for every mode instead. */
  CKerr_OrbitCache orbit_cache;
  CKerr_OrbitCacheBuild(Minv, xuorig, M, astar, Nr, Nt, (mmax < 2 ? 2 : mmax) + nl - 1, &orbit_cache);
#endif

  /* This takes the M_inv, the origin in phase space, mass and spin of BH,
  {n,k,m}, number of steps in r and theta, number of modes (nl), and the angular
  frequency of the perturber orbit Returns the Z-coeff. and energy eigenvalues
  E[nl]. */

#pragma omp parallel
  {
    // --------------------------------------------------------
    // Thread-private working buffers
    // Allocated ONCE per thread (not per iteration)
    // --------------------------------------------------------
    // int tid = omp_get_thread_num();     // Get thread ID
    // double thread_t0 = omp_get_wtime(); // Start time at each thread
    double *C0, *E0, *E1;

    E0 = (double *)malloc((size_t)(nl * 10 * sizeof(double)));

    // Safety check that E0 was made correctly
    if (!E0) {
      printf("Memory allocation failed in thread %d\n", omp_get_thread_num());
      exit(1);
    }

    E1 = E0 + nl;
    C0 = E1 + nl;

    // Initialize the values of J_dot_sf to be zero
    double J_dot_sf_local[3] = {0.0, 0.0, 0.0};

    // Time the CKerrRadialFunc
    double radial_time_local = 0.0;
    long radial_calls_local = 0;

// Reduce the nested for loops over the indices and collapse to
// parallelize amongst all indices
#pragma omp for collapse(3) schedule(dynamic, 1)
    for (int i_n = -nmax; i_n <= nmax; i_n++) {
      for (int i_k = -kmax; i_k <= kmax; i_k++) {
        for (int i_m = -mmax; i_m <= mmax; i_m++) {
          if (i_n == 0 && i_k == 0 &&
              i_m == 0) { // Remove 0 mode (time indepedent perturbation)
            continue;
          }
          if (astar == 0 && i_n == 0 && i_k == -i_m) // Schwarzschild case a ==
                                                     // 0
            continue;
          /* (-n,-k,-m) has omega -> -omega and Z -> (-1)^l conj(Z), so its term equals the (n,k,m) term:
           * keep the half whose first nonzero index is positive and count it twice */
          if (i_n < 0 || (i_n == 0 && (i_k < 0 || (i_k == 0 && i_m < 0))))
            continue;
          double omegagw;

          // double radial_t0 = omp_get_wtime();
#ifdef JDOT_SF_ORIGINAL
          CKerr_RadialFunc(Minv, xuorig, M, astar, i_n, i_k, i_m, Nr, Nt, nl,
                           C0, &omegagw, E0);
#else
          CKerr_RadialFunc_cached(Minv, &orbit_cache, M, astar, i_n, i_k, i_m, nl,
                                  C0, &omegagw, E0);
#endif
          // double radial_t1 = omp_get_wtime();
          // radial_time_local +=
          //     (radial_t1 - radial_t0); // Time CKerrRadialFunc per thread
          // radial_calls_local++;

          for (int il = 0; (il) < nl; il++) {

            // omega_nkm = i_n * Omega[0] + i_k * Omega[1] + i_m * Omega[2];
            double omega_nkm = omegagw;

            double lambda = E0[il] - 2 * astar * M * i_m * omega_nkm +
                            astar * M * astar * M * omega_nkm * omega_nkm - 2;

            double P = omega_nkm - i_m * astar * M / (2 * M * rH);

            double numer = 256 * pow(2 * M * rH, 5) * P *
                           (P * P + 4 * epsilon * epsilon) *
                           (P * P + 16 * epsilon * epsilon) * omega_nkm *
                           omega_nkm * omega_nkm;

            double C2 =
                ((lambda + 2) * (lambda + 2) + 4 * i_m * astar * omega_nkm -
                 4 * astar * M * astar * M * omega_nkm * omega_nkm) *
                    (lambda * lambda + 36 * i_m * astar * M * omega_nkm -
                     36 * astar * M * astar * M * omega_nkm * omega_nkm) +
                (2 * lambda + 3) *
                    (96 * astar * M * astar * M * omega_nkm * omega_nkm -
                     48 * i_m * astar * M * omega_nkm) +
                144 * omega_nkm * omega_nkm * (M * M - astar * astar * M * M);

            double alphankm = numer / C2;

            double Z_down_square =
                C0[4 * il] * C0[4 * il] + C0[4 * il + 1] * C0[4 * il + 1];
            double Z_out_square = C0[4 * il + 2] * C0[4 * il + 2] +
                                  C0[4 * il + 3] * C0[4 * il + 3];
#if 0
							printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n, i_k, i_m, il, omegagw, omega_nkm, lambda, alphankm, Z_down_square, Z_out_square);
							fflush(stdout);
							// printf("%i \t %i \t %i \t %i \t %lg \t %lg \n", i_n, i_k, i_m, il, Z_down_square, Z_out_square);
#endif
            // printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg
            // \t %lg\n", i_n, i_k, i_m, il, C0[4*il], C0[4*il+1], C0[4*il+2],
            // C0[4*il+3], Z_down_square, Z_out_square);

            /* Factor 2 counts the (-n,-k,-m) mode skipped above */
            J_dot_sf_local[0] += 2. * mu_body * -i_n *
                                 (Z_out_square + alphankm * Z_down_square) /
                                 (2. * omega_nkm * omega_nkm * omega_nkm);
            J_dot_sf_local[1] += 2. * mu_body * -i_k *
                                 (Z_out_square + alphankm * Z_down_square) /
                                 (2. * omega_nkm * omega_nkm * omega_nkm);
            J_dot_sf_local[2] += 2. * mu_body * -i_m *
                                 (Z_out_square + alphankm * Z_down_square) /
                                 (2. * omega_nkm * omega_nkm * omega_nkm);
          }
        }
      }
    }

    // double thread_t1 = omp_get_wtime(); // End time of thread calculation
    // double loop_time = thread_t1 - thread_t0;

    // --------------------------------------------------------
    // Merge thread-local results
    // --------------------------------------------------------

#pragma omp critical
    {

      #if 0
      printf("Thread %d total time = %.6f s\n", tid,
             thread_t1 - thread_t0); // Thread ID and how long calculation took

      printf("Thread %d calls = %ld\n", tid, radial_calls_local);
      printf("  Thread time   = %.3f s\n", loop_time);
      printf("  Radial time = %.3f s\n", radial_time_local);
      printf("  Fraction    = %.1f %%\n",
             100.0 * radial_time_local / loop_time);

      radial_time_global += radial_time_local;
      radial_calls_global += radial_calls_local;
      #endif

      J_dot_sf_global[0] += J_dot_sf_local[0];
      J_dot_sf_global[1] += J_dot_sf_local[1];
      J_dot_sf_global[2] += J_dot_sf_local[2];
    }
    free((char *)E0);
  }
#ifndef JDOT_SF_ORIGINAL
  CKerr_OrbitCacheFree(&orbit_cache);
#endif
  // ------------------------------------------------------------
  // Final output
  // ------------------------------------------------------------

  #if 0
  double total_t1 = omp_get_wtime();

  double total_runtime = total_t1 - total_t0;

  printf("\n");
  printf("Total runtime              = %.6f s\n", total_runtime);

  printf("CKerr_RadialFunc CPU time  = %.6f s\n", radial_time_global);

  printf("Number of calls            = %ld\n", radial_calls_global);

  printf("Average call time          = %.6e s\n",
         radial_time_global / radial_calls_global);
  #endif

  J_dot_sf[0] = J_dot_sf_global[0];
  J_dot_sf[1] = J_dot_sf_global[1];
  J_dot_sf[2] = J_dot_sf_global[2];
}

/* Computing the J_dot_tidal,i */
void J_dot_tidal(int nl, int N_res, int n_res_inner, int n_res_outer,
                 int k_res_inner, int k_res_outer, int m_res_inner,
                 int m_res_outer, double ra_inner, double rp_inner,
                 double radius_outer, double I_inner, double ra_outer,
                 double rp_outer, double I_outer, double M, double astar,
                 double theta_res_F, double mu_outer, double *J_dot_td) {
  int i, i_n_inner, i_k_inner, i_m_inner, il;
  int i_n_outer, i_k_outer, i_m_outer;
  double EQL_inner[3], J_inner[3], EQL_outer[3], J_outer[3], Minv_inner[9],
      Minv_outer[9], Omega_inner[3], info[6], info_outer[6], xuorig_inner[6],
      xuorig_outer[6], cscat[16], aux[4];
  double term, another_term, last_term;
  double Gamma, sgn_Gamma;
  double omega_nkm, omegagw_inner, omegagw_outer;
  double rH;
  double Rtheta, Itheta;
  double epsilon, lambda, numer, C2, P, alphankm;
  double *C0_inner, *E0_inner, *E1_inner, *C0_outer, *E0_outer, *E1_outer;

  /* Converting apocenter, pericenter, and inclination angle into frequencies
   * for: */
  /* Inner body */
  ra_rp_I2EQL(ra_inner, EQL_inner, rp_inner, I_inner, astar, M);
  CKerr_EQL2J(EQL_inner, J_inner, M, astar, NULL);
  CKerr_Minverse(J_inner, Minv_inner, M, astar);
  CKerr_Minv2Omega(Minv_inner, Omega_inner);
  // printf("EQL_inner: %lg %lg %lg \n", EQL_inner[0], EQL_inner[1],
  // EQL_inner[2]); printf("Minv for inner body = \n %lg %lg %lg \n",
  // Minv_inner[0], Minv_inner[1], Minv_inner[2]); printf(" %lg %lg %lg \n",
  // Minv_inner[3], Minv_inner[4], Minv_inner[5]); printf(" %lg %lg %lg \n",
  // Minv_inner[6], Minv_inner[7], Minv_inner[8]);

  // printf("dtau/dt = %lg \n", EQL_inner[0] - Omega_inner[0]*J_inner[0] -
  // Omega_inner[1]*J_inner[1] - Omega_inner[2]*J_inner[2]);
  /* Outer Body */
  // CKerr_FindEQL_IRCirc(0, radius_outer, EQL_outer, M, astar);
  // CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
  // CKerr_Minverse(J_outer, Minv_outer, M, astar);
  // //printf("Minv for outer body = \n %lg %lg %lg \n", Minv_outer[0],
  // Minv_outer[1], Minv_outer[2]);
  // //printf(" %lg %lg %lg \n", Minv_outer[3], Minv_outer[4], Minv_outer[5]);
  // //printf(" %lg %lg %lg \n", Minv_outer[6], Minv_outer[7], Minv_outer[8]);

  // /* Torus origin for inner body orbit */
  // CKerr_TorusOrigin(J_inner, xuorig_inner, M, astar);

  // /* Torus origin for outer body orbit */
  // CKerr_getData_CircEq(M, astar, radius_outer, info_outer);
  // xuorig_outer[0] = radius_outer; xuorig_outer[1] = M_PI/2.; xuorig_outer[2]
  // = 0.; xuorig_outer[3] = 0.; xuorig_outer[4] = 0.;      xuorig_outer[5] =
  // info_outer[0];

  /* Converting apocenter, pericenter, and inclination angle into frequencies
   * for: */
  /* Outer body, both cases*/
  /* For CIRCULAR orbits */
  if (ra_outer == 0 && rp_outer == 0 && I_outer == 0) {
    /* Orbital frequencies for circular equitorial orbits */
    CKerr_FindEQL_IRCirc(0, radius_outer, EQL_outer, M, astar);
    CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
    CKerr_Minverse(J_outer, Minv_outer, M, astar);
    // CKerr_Minv2Omega(Minv_outer, Omega_outer);
    // printf("Minv for circular equitorial case = \n %lg %lg %lg \n", Minv[0],
    // Minv[1], Minv[2]); printf(" %lg %lg %lg \n", Minv[3], Minv[4], Minv[5]);
    // printf(" %lg %lg %lg \n", Minv[6], Minv[7], Minv[8]);
    // printf("Js are: %lg %lg %lg \n", J[0], J[1], J[2]);
    // printf("Omega for circular equatorial is: %lg %lg %lg \n", Omega[0],
    // Omega[1], Omega[2]);

    /* Torus Origin for circular equatorial orbits */
    // CKerr_TorusOrigin(J, xuorig, M, astar);
  }

  /* For GENERIC orbits */
  else {
    /* Routine to find location of resonant apocenter */
    // guess1 = rp + 3.;
    // guess2 = radius_outer - 3.;
    // apo = find_resonance_apo(n_res_inner, k_res_inner, m_res_inner,
    // radius_outer, guess1, guess2, rp, I, astar, M);

    /* Get orbital frequencies for inner generic orbit */
    ra_rp_I2EQL(ra_outer, EQL_outer, rp_outer, I_outer, astar, M);
    CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
    CKerr_Minverse(J_outer, Minv_outer, M, astar);
    // CKerr_Minv2Omega(Minv_outer, Omega_outer);
    // printf("Minv for generic case = \n %lg %lg %lg \n", Minv[0], Minv[1],
    // Minv[2]); printf(" %lg %lg %lg \n", Minv[3], Minv[4], Minv[5]); printf("
    // %lg %lg %lg \n", Minv[6], Minv[7], Minv[8]); printf("Js are: %lg %lg %lg
    // \n", J[0], J[1], J[2]); printf("Omega for generic is: %lg %lg %lg \n",
    // Omega[0], Omega[1], Omega[2]);

    /* Torus origin for inner body orbit */
    // CKerr_TorusOrigin(J, xuorig, M, astar);
    // for(i=0;i<6;i++) printf("xuorig[%d] = %lg\n", i, xuorig[i]);
  }
#if 0
	/* Converting apocenter, pericenter, and inclination angle into frequencies for: */
	/* Outer body */
	
	ra_rp_I2EQL(ra_outer, EQL_outer, rp_outer, I_outer, astar, M);
	CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
	CKerr_Minverse(J_outer, Minv_outer, M, astar);
	CKerr_Minv2Omega(Minv_outer, Omega_outer);
#endif

#if 0
	/* Outer Body */
	CKerr_FindEQL_IRCirc(0, radius_outer, EQL_outer, M, astar);
	CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
	CKerr_Minverse(J_outer, Minv_outer, M, astar);
	//printf("Minv for outer body = \n %lg %lg %lg \n", Minv_outer[0], Minv_outer[1], Minv_outer[2]);
	//printf(" %lg %lg %lg \n", Minv_outer[3], Minv_outer[4], Minv_outer[5]);
	//printf(" %lg %lg %lg \n", Minv_outer[6], Minv_outer[7], Minv_outer[8]);
	/* Torus origin for outer body orbit */
	CKerr_getData_CircEq(M, astar, radius_outer, info_outer);
	xuorig_outer[0] = radius_outer; xuorig_outer[1] = M_PI/2.; xuorig_outer[2] = 0.;
  	xuorig_outer[3] = 0.; xuorig_outer[4] = 0.;      xuorig_outer[5] = info_outer[0];
#endif

  /* Torus origin for inner/outer body orbit */
  CKerr_TorusOrigin(J_inner, xuorig_inner, M, astar);
  CKerr_TorusOrigin(J_outer, xuorig_outer, M, astar);

  E0_inner = (double *)malloc((size_t)(nl * 10 * sizeof(double)));
  E1_inner = E0_inner + nl;
  C0_inner = E1_inner + nl;

  E0_outer = (double *)malloc((size_t)(nl * 10 * sizeof(double)));
  E1_outer = E0_outer + nl;
  C0_outer = E1_outer + nl;

  rH = M + sqrt(M * M - astar * astar * M * M);
  epsilon = sqrt(M * M - astar * astar * M * M) / (4 * M * rH);
  // Gamma = omega_dot(nl, n_res_inner, k_res_inner, m_res_inner, apo, rp, I,
  // astar, M, radius_outer, 1e-4); sgn_Gamma = Gamma/fabs(Gamma); printf("Gamma
  // and its sign are: %lg %lg \n", Gamma, sgn_Gamma);

#if 0
	Delta_J_r_tidal[0] = 0;
	Delta_J_r_tidal[1] = 0; 
	Delta_J_theta_tidal[0] = 0;
	Delta_J_theta_tidal[1] = 0;	
	Delta_J_phi_tidal[0] = 0;
	Delta_J_phi_tidal[1] = 0;
#endif
  J_dot_td[0] = 0;
  J_dot_td[1] = 0;
  J_dot_td[2] = 0;

  // printf("About to start for-loops \n");
  // printf("il \t i_n \t i_k \t i_m \t C01_in \t C02_in \t C03_in \t C04_in \t
  // C01_out \t C02_out \t term \t another_term \t last_term \n");
  /* Resonant sum, so it is over integer mulitples of the initial resonance
   * modes for the inner body */
  for (i = -N_res; i <= N_res; i++) {
    i_n_inner = i * n_res_inner;
    i_k_inner = i * k_res_inner;
    i_m_inner = i * m_res_inner;
    i_n_outer = i * n_res_outer;
    i_k_outer = i * k_res_outer;
    i_m_outer = i * m_res_outer;
    Rtheta = cos(i * theta_res_F);
    Itheta = sin(i * theta_res_F);
    if ((i_n_outer == 0 && i_k_outer == 0 && i_m_outer == 0) ||
        (i_n_inner == 0 && i_k_inner == 0 && i_m_inner == 0) ||
        abs(i_n_outer) > 20 || abs(i_k_outer) > 20 || abs(i_n_inner) > 20 ||
        abs(i_k_inner) > 20)
      continue;
    /* Create the amplitude data for inner body*/
    /* Create the amplitude data for outer body*/
    CKerr_RadialFunc(Minv_outer, xuorig_outer, M, astar, i_n_outer, i_k_outer,
                     i_m_outer, 42, 42, nl, C0_outer, &omegagw_outer, E0_outer);
    CKerr_RadialFunc(Minv_inner, xuorig_inner, M, astar, i_n_inner, i_k_inner,
                     i_m_inner, 42, 42, nl, C0_inner, &omegagw_inner, E0_inner);
    for (il = 0; (il) < nl; il++) {
      omega_nkm = i_n_inner * Omega_inner[0] + i_k_inner * Omega_inner[1] +
                  i_m_inner * Omega_inner[2];
      /* Get the scattering coefficients of the outer body radiation */
      CKerr_GWScatMatrix(M, astar, omega_nkm, i_m_inner, E0_outer[il], cscat,
                         aux);
      /* Define frequency of inner body (omega_nkm) and other functions that
       * depends on n,k,m,l */
      lambda = E0_inner[il] - 2 * astar * M * i_m_inner * omega_nkm +
               astar * M * astar * M * omega_nkm * omega_nkm - 2;
      P = omega_nkm - i_m_inner * astar * M / (2 * M * rH);
      numer = 256 * pow(2 * M * rH, 5) * P * (P * P + 4 * epsilon * epsilon) *
              (P * P + 16 * epsilon * epsilon) * omega_nkm * omega_nkm *
              omega_nkm;
      C2 =
          ((lambda + 2) * (lambda + 2) + 4 * i_m_inner * astar * M * omega_nkm -
           4 * astar * astar * M * M * omega_nkm * omega_nkm) *
              (lambda * lambda + 36 * i_m_inner * astar * M * omega_nkm -
               36 * astar * astar * M * M * omega_nkm * omega_nkm) +
          (2 * lambda + 3) *
              (96 * astar * astar * M * M * omega_nkm * omega_nkm -
               48 * i_m_inner * astar * M * omega_nkm) +
          144 * omega_nkm * omega_nkm * (M * M - astar * astar * M * M);
      alphankm = numer / C2;
#if 0
			printf("Frequencies from resonance and RadialFunc: %lg \t %lg \t %lg\n", omega_nkm, omegagw_inner, omegagw_outer);
			printf("%i \t %lg \t %lg \t %lg\n", il, omega_nkm, lambda, alphankm);

			printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n_inner, i_k_inner, i_m_inner, il,  C0_inner[4*il], C0_inner[4*il+1], C0_inner[4*il+2], C0_inner[4*il+3], C0_outer[4*il], C0_outer[4*il+1], C0_outer[4*il+2], C0_outer[4*il+3]);
			//printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n, i_k, i_m, il, C0[4*il], C0[4*il+1], C0[4*il+2], C0[4*il+3], Z_down_square, Z_out_square);
#endif
      /* For INNER/OUTER Body: */
      /* C0_{inner/outer}[4*il    ] = Re[Z_down, inner] */
      /* C0_{inner/outer}[4*il + 1] = Im[Z_down, inner] */
      /* C0_{inner/outer}[4*il + 2] = Re[Z_out, inner] */
      /* C0_{inner/outer}[4*il + 3] = Im[Z_out, inner] */

      /* Expanded out the terms in our expression for \dot{J}_{td} and kept only
       * the real parts*/
      term = C0_inner[4 * il + 2] * cscat[0] * Rtheta * C0_outer[4 * il] +
             C0_outer[4 * il + 1] * cscat[0] * C0_inner[4 * il + 2] * Itheta +
             C0_outer[4 * il] * cscat[1] * C0_inner[4 * il + 3] * Rtheta +
             C0_outer[4 * il + 1] * cscat[1] * C0_inner[4 * il + 3] * Itheta;
      another_term =
          -C0_outer[4 * il] * C0_inner[4 * il + 3] * cscat[0] * Itheta +
          C0_outer[4 * il + 1] * C0_inner[4 * il + 3] * cscat[0] * Rtheta +
          C0_outer[4 * il] * cscat[1] * C0_inner[4 * il + 2] * Itheta -
          C0_outer[4 * il + 1] * cscat[1] * C0_inner[4 * il + 2] * Rtheta;
      last_term = C0_outer[4 * il] * C0_inner[4 * il] * Rtheta +
                  C0_outer[4 * il + 1] * C0_inner[4 * il] * Itheta -
                  C0_outer[4 * il] * C0_inner[4 * il + 1] * Itheta +
                  C0_outer[4 * il + 1] * C0_inner[4 * il + 1] * Rtheta;
      // printf("%lg \t %lg \t %lg \n", term, another_term, alphankm *
      // last_term);

      /* J_dot of inner body due to tidal field of outer body */
      J_dot_td[0] += -i_n_inner * mu_outer *
                     (term + another_term + alphankm * last_term) /
                     (omega_nkm * omega_nkm * omega_nkm);
      J_dot_td[1] += -i_k_inner * mu_outer *
                     (term + another_term + alphankm * last_term) /
                     (omega_nkm * omega_nkm * omega_nkm);
      J_dot_td[2] += -i_m_inner * mu_outer *
                     (term + another_term + alphankm * last_term) /
                     (omega_nkm * omega_nkm * omega_nkm);
// printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg
// \t %lg \t %lg \n", il, i_n_inner, i_k_inner, i_m_inner, C0_inner[4*il],
// C0_inner[4*il+1], C0_inner[4*il+2], C0_inner[4*il+3], C0_outer[4*il],
// C0_outer[4*il+1], term, another_term, alphankm * last_term); printf("%lg \t
// %lg \t %lg \n", *J_dot_r_tidal, *J_dot_theta_tidal, *J_dot_phi_tidal);
#if 0
			/* \Delta J_td of inner body due to tidal field of outer body, for the real [0] and imaginary [1] parts */
			Delta_J_r_tidal[0] += -i_n_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm) * cos(sgn_Gamma * M_PI/4.) * sqrt(2. * M_PI / fabs(Gamma));
			Delta_J_r_tidal[1] += -i_n_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm) * sin(sgn_Gamma * M_PI/4.) * sqrt(2. * M_PI / fabs(Gamma));
			Delta_J_theta_tidal[0] += -i_k_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm) * cos(sgn_Gamma * M_PI/4.) * sqrt(2. * M_PI / fabs(Gamma));
			Delta_J_theta_tidal[1] += -i_k_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm) * sin(sgn_Gamma * M_PI/4.) * sqrt(2. * M_PI / fabs(Gamma));
			Delta_J_phi_tidal[0] += -i_m_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm) * cos(sgn_Gamma * M_PI/4.) * sqrt(2. * M_PI / fabs(Gamma));
			Delta_J_phi_tidal[1] += -i_m_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm) * sin(sgn_Gamma * M_PI/4.) * sqrt(2. * M_PI / fabs(Gamma));
			printf("%lg \t %lg \t %lg \t %lg \t %lg \t %lg \n", Delta_J_r_tidal[0], Delta_J_r_tidal[1], Delta_J_theta_tidal[0], Delta_J_theta_tidal[1], Delta_J_phi_tidal[0], Delta_J_phi_tidal[1]);
#endif
    }
  }
  free((char *)E0_inner);
  free((char *)E0_outer);
}

/*
------------------------------------------------------------
Compute J_dot_tidal using OpenMP parallelization
over independent resonance modes (i-loop)
------------------------------------------------------------
This is the primary function to use
J_dot_inner -> EMRI is interior to perturber, scaled(J_dot_tidal_inner) ~ mu_outer
J_dot_outer -> EMRI is exterior to perturber, scaled(J_dot_tidal_outer) ~ mu_inner
*/

/* Modes skipped in J_dot_tidal_openmp: zero modes, or |n|, |k| > 20 for either body */
static int tidal_mode_skipped(int i, int n_res_inner, int k_res_inner, int m_res_inner,
                              int n_res_outer, int k_res_outer, int m_res_outer) {
  int i_n_inner = i * n_res_inner, i_k_inner = i * k_res_inner, i_m_inner = i * m_res_inner;
  int i_n_outer = i * n_res_outer, i_k_outer = i * k_res_outer, i_m_outer = i * m_res_outer;
  return (i_n_outer == 0 && i_k_outer == 0 && i_m_outer == 0) ||
         (i_n_inner == 0 && i_k_inner == 0 && i_m_inner == 0) ||
         abs(i_n_outer) > 20 || abs(i_k_outer) > 20 || abs(i_n_inner) > 20 ||
         abs(i_k_inner) > 20;
}

void J_dot_tidal_openmp(int nl, int N_res, int n_res_inner, int n_res_outer,
                        int k_res_inner, int k_res_outer, int m_res_inner,
                        int m_res_outer, double ra_inner, double rp_inner,
                        double radius_outer, double I_inner, double ra_outer,
                        double rp_outer, double I_outer, double M, double astar,
                        double theta_res_F, double mu_inner, double mu_outer, 
                        double *J_dot_td_inner, double *J_dot_td_outer) {
  // --------------------------------------------------------
  // Global accumulation (final result)
  // --------------------------------------------------------
  double J_global_inner[3] = {0.0, 0.0, 0.0};
  double J_global_outer[3] = {0.0, 0.0, 0.0};

  // --------------------------------------------------------
  // Precompute orbit quantities (NOT parallel region)
  // --------------------------------------------------------
  double EQL_inner[3], J_inner[3], Minv_inner[9], Omega_inner[3];
  double EQL_outer[3], J_outer[3], Minv_outer[9];
  double JRe_global_inner[3] = {0.0, 0.0, 0.0}, JIm_global_inner[3] = {0.0, 0.0, 0.0};
  double JRe_global_outer[3] = {0.0, 0.0, 0.0}, JIm_global_outer[3] = {0.0, 0.0, 0.0};

  // Inner body apocenter, pericenter, inclination -> EQL and Minv
  ra_rp_I2EQL(ra_inner, EQL_inner, rp_inner, I_inner, astar, M);
  for (int i = 0; i < 3; i++) {
    	if (!isfinite(EQL_inner[i])) {
        printf("J_dot.c: No values of EQL_inner found, exiting...\n");
        exit(1);
		}
	}

  CKerr_EQL2J(EQL_inner, J_inner, M, astar, NULL);
  for (int i = 0; i < 3; i++) {
    	if (!isfinite(J_inner[i])) {
        printf("J_dot.c: No values of J_inner found, exiting...\n");
        exit(1);
		}
	}

  CKerr_Minverse(J_inner, Minv_inner, M, astar);
  // CKerr_Minv2Omega(Minv_inner, Omega_inner);

  ra_rp_I2EQL(ra_outer, EQL_outer, rp_outer, I_outer, astar, M);
  for (int i = 0; i < 3; i++) {
    	  if (!isfinite(EQL_outer[i])) {
          printf("J_dot.c: No values of EQL_outer found, exiting...\n");
          exit(1);
	    }
	  }

  CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
  for (int i = 0; i < 3; i++) {
    	if (!isfinite(J_outer[i])) {
        printf("J_dot.c: No values of J_outer found, exiting...\n");
        exit(1);
		}
	}
  
  CKerr_Minverse(J_outer, Minv_outer, M, astar);

  // Set up initial conditions to integrate on tori
  double xuorig_inner[6], xuorig_outer[6];
  CKerr_TorusOrigin(J_inner, xuorig_inner, M, astar);
  CKerr_TorusOrigin(J_outer, xuorig_outer, M, astar);

  double rH = M + sqrt(M * M - astar * astar * M * M);
  double epsilon = sqrt(M * M - astar * astar * M * M) / (4 * M * rH);

// --------------------------------------------------------
// Main OpenMP region
// --------------------------------------------------------
#if 0
  /* Original version: one task per i = -N_res..N_res, each calling CKerr_RadialFunc for both bodies in sequence.
   * Replaced by the two-step version below (modes i > 0 counted twice, one task per RadialFunc call). */
#pragma omp parallel
  {
    // --------------------------------------------------------
    // Thread-private working buffers
    // Allocated ONCE per thread (not per iteration)
    // --------------------------------------------------------
    // int tid = omp_get_thread_num();     // Get thread ID
    // double thread_t0 = omp_get_wtime(); // Start time at each thread

    double *C0_inner, *E0_inner, *E1_inner, *C0_outer, *E0_outer, *E1_outer;
    E0_inner = (double *)malloc((size_t)(nl * 10 * sizeof(double)));
    E1_inner = E0_inner + nl;
    C0_inner = E1_inner + nl;

    E0_outer = (double *)malloc((size_t)(nl * 10 * sizeof(double)));
    E1_outer = E0_outer + nl;
    C0_outer = E1_outer + nl;

    // Safety check (good HPC habit)
    if (!E0_inner || !E0_outer) {
      printf("Memory allocation failed in thread %d\n", omp_get_thread_num());
      exit(1);
    }
    double J_local_inner[3] = {0.0, 0.0, 0.0};
    double J_local_outer[3] = {0.0, 0.0, 0.0};
    // double JRe_local_inner[3] = {0.0, 0.0, 0.0};
    // double JIm_local_inner[3] = {0.0, 0.0, 0.0};
    double JRe_local_outer[3] = {0.0, 0.0, 0.0};
    double JIm_local_outer[3] = {0.0, 0.0, 0.0};

#pragma omp for schedule(dynamic, 1)
    for (int i = -N_res; i <= N_res; i++) {

      int i_n_inner = i * n_res_inner;
      int i_k_inner = i * k_res_inner;
      int i_m_inner = i * m_res_inner;

      int i_n_outer = i * n_res_outer;
      int i_k_outer = i * k_res_outer;
      int i_m_outer = i * m_res_outer;

      // Skip invalid / out-of-range modes
      if ((i_n_outer == 0 && i_k_outer == 0 && i_m_outer == 0) ||
          (i_n_inner == 0 && i_k_inner == 0 && i_m_inner == 0) ||
          abs(i_n_outer) > 20 || abs(i_k_outer) > 20 || abs(i_n_inner) > 20 ||
          abs(i_k_inner) > 20) {
        continue;
      }

      double Rtheta = cos(i * theta_res_F);
      double Itheta = sin(i * theta_res_F);

      double omegagw_inner, omegagw_outer;
      double cscat[16], aux[4];

      // Radial functions (mode-dependent setup)
      CKerr_RadialFunc(Minv_outer, xuorig_outer, M, astar, i_n_outer, i_k_outer,
                       i_m_outer, 42, 42, nl, C0_outer, &omegagw_outer,
                       E0_outer);

      CKerr_RadialFunc(Minv_inner, xuorig_inner, M, astar, i_n_inner, i_k_inner,
                       i_m_inner, 42, 42, nl, C0_inner, &omegagw_inner,
                       E0_inner);

      // double omega_nkm =
      //         i_n_inner * Omega_inner[0] +
      //         i_k_inner * Omega_inner[1] +
      //         i_m_inner * Omega_inner[2];

      // Set the frequency to the inner body value,
      // but either body is fine since the system is in resonance
      double omega_nkm = omegagw_inner;

      double denom = fabs(omegagw_inner) +
                     fabs(omegagw_outer); // absolute sum of both outputs from
                                          // inner and outer bodies

      double rel_diff;

      if (denom > 0.0) {
        rel_diff = 2.0 * fabs(omegagw_inner - omegagw_outer) / denom;
      } else {
        rel_diff = 0.0; // both are zero
      }

      double delta_omega_tol = 1.0e-4; // 0.1%

      if (rel_diff > delta_omega_tol) {
        printf("WARNING: GW frequencies for inner and outer bodies differ by "
               "%.6e (tol %.6e)\n",
               rel_diff, delta_omega_tol);
      }

      // ----------------------------------------------------
      // Inner sum over l
      // ----------------------------------------------------
      for (int il = 0; il < nl; il++) {

        CKerr_GWScatMatrix(M, astar, omega_nkm, i_m_inner, E0_outer[il], cscat,
                           aux);

        double lambda = E0_inner[il] - 2 * astar * M * i_m_inner * omega_nkm +
                        astar * astar * M * M * omega_nkm * omega_nkm - 2;

        double P =
            omega_nkm - i_m_inner * astar * M /
                            (2 * M * (M + sqrt(M * M - astar * astar * M * M)));

        double numer = 256 * pow(2 * M * rH, 5) * P *
                       (P * P + 4 * epsilon * epsilon) *
                       (P * P + 16 * epsilon * epsilon) * omega_nkm *
                       omega_nkm * omega_nkm;

        double C2 =
            ((lambda + 2) * (lambda + 2) +
             4 * i_m_inner * astar * M * omega_nkm -
             4 * astar * astar * M * M * omega_nkm * omega_nkm) *
                (lambda * lambda + 36 * i_m_inner * astar * M * omega_nkm -
                 36 * astar * astar * M * M * omega_nkm * omega_nkm) +
            (2 * lambda + 3) *
                (96 * astar * astar * M * M * omega_nkm * omega_nkm -
                 48 * i_m_inner * astar * M * omega_nkm) +
            144 * omega_nkm * omega_nkm * (M * M - astar * astar * M * M);

        double alphankm = numer / C2;

#if 0
				printf("Frequencies from resonance and RadialFunc: %lg \t %lg \t %lg\n", omega_nkm, omegagw_inner, omegagw_outer);
				printf("%i \t %lg \t %lg \t %lg\n", il, omega_nkm, lambda, alphankm);

				printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n_inner, i_k_inner, i_m_inner, il,  C0_inner[4*il], C0_inner[4*il+1], C0_inner[4*il+2], C0_inner[4*il+3], C0_outer[4*il], C0_outer[4*il+1], C0_outer[4*il+2], C0_outer[4*il+3]);
				//printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n, i_k, i_m, il, C0[4*il], C0[4*il+1], C0[4*il+2], C0[4*il+3], Z_down_square, Z_out_square);
#endif
        /* REAL PART of J_dot_tidal_i */
        /* For INNER/OUTER Body: */
        /* C0_{inner/outer}[4*il    ] = Re[Z_down, inner] */
        /* C0_{inner/outer}[4*il + 1] = Im[Z_down, inner] */
        /* C0_{inner/outer}[4*il + 2] = Re[Z_out, inner] */
        /* C0_{inner/outer}[4*il + 3] = Im[Z_out, inner] */
        
        #if 0
        double term =
            C0_inner[4 * il + 2] * cscat[0] * Rtheta * C0_outer[4 * il] +
            C0_outer[4 * il + 1] * cscat[0] * C0_inner[4 * il + 2] * Itheta +
            C0_outer[4 * il] * cscat[1] * C0_inner[4 * il + 3] * Rtheta +
            C0_outer[4 * il + 1] * cscat[1] * C0_inner[4 * il + 3] * Itheta;

        double another_term =
            -C0_outer[4 * il] * C0_inner[4 * il + 3] * cscat[0] * Itheta +
            C0_outer[4 * il + 1] * C0_inner[4 * il + 3] * cscat[0] * Rtheta +
            C0_outer[4 * il] * cscat[1] * C0_inner[4 * il + 2] * Itheta +
            -C0_outer[4 * il + 1] * cscat[1] * C0_inner[4 * il + 2] * Rtheta;

        double last_term = C0_outer[4 * il] * C0_inner[4 * il] * Rtheta +
                           C0_outer[4 * il + 1] * C0_inner[4 * il] * Itheta -
                           C0_outer[4 * il] * C0_inner[4 * il + 1] * Itheta +
                           C0_outer[4 * il + 1] * C0_inner[4 * il + 1] * Rtheta;

        double contrib = (term + another_term + alphankm * last_term) /
                         (omega_nkm * omega_nkm * omega_nkm);
        #endif
        
        /* J_dot_tidal_inner */
        double termR1_inner = 
              C0_inner[4 * il + 2] * cscat[0] * C0_outer[4 * il] +
              C0_outer[4 * il] * cscat[1] * C0_inner[4 * il + 3] +
              C0_outer[4 * il + 1] * C0_inner[4 * il + 3] * cscat[0] -
              C0_outer[4 * il + 1] * cscat[1] * C0_inner[4 * il + 2] +
              alphankm * (C0_outer[4 * il] * C0_inner[4 * il] + 
                C0_outer[4 * il + 1] * C0_inner[4 * il + 1]);
        
        double termR2_inner = 
              C0_outer[4 * il + 1] * cscat[0] * C0_inner[4 * il + 2] +
              C0_outer[4 * il + 1] * cscat[1] * C0_inner[4 * il + 3] -
              C0_outer[4 * il] * C0_inner[4 * il + 3] * cscat[0] +
              C0_outer[4 * il] * cscat[1] * C0_inner[4 * il + 2] +
              alphankm * (C0_outer[4 * il + 1] * C0_inner[4 * il] -
              C0_outer[4 * il] * C0_inner[4 * il + 1]);
        
        double contrib_Re_inner = (termR1_inner * Rtheta + termR2_inner * Itheta) / (omega_nkm * omega_nkm * omega_nkm);
        printf("J_dot.c: contrib_Re_inner = %.15e \n", contrib_Re_inner);
        
        #if 0
        double termIm1_inner = termR1_inner * Itheta;
        double termIm2_inner = -termR2_inner * Rtheta;
        double contrib_Im_inner = (termIm1_inner + termIm2_inner) / (omega_nkm * omega_nkm * omega_nkm);
        #endif

        /* TODO: J_dot_tidal_outer */
        double termR1_outer = 
        C0_inner[4 * il + 2] * C0_outer[4 * il + 2] +
        C0_inner[4 * il + 3] * C0_outer[4 * il + 3] -
        C0_inner[4 * il + 2] * cscat[0] * C0_outer[4 * il] +
        C0_inner[4 * il + 2] * cscat[1] * C0_outer[4 * il + 1] -
        C0_inner[4 * il + 3] * cscat[1] * C0_outer[4 * il] -
        C0_inner[4 * il + 3] * cscat[0] * C0_outer[4 * il + 1];

        double termR2_outer = 
        C0_inner[4 * il + 2] * C0_outer[4 * il + 3] -
        C0_inner[4 * il + 3] * C0_outer[4 * il + 2] +
        C0_inner[4 * il + 3] * cscat[0] * C0_outer[4 * il] -
        C0_inner[4 * il + 3] * cscat[1] * C0_outer[4 * il + 1] -
        C0_inner[4 * il + 2] * cscat[1] * C0_outer[4 * il] -
        C0_inner[4 * il + 2] * cscat[0] * C0_outer[4 * il + 1];

        double contrib_Re_outer = (termR1_outer * Rtheta + termR2_outer * Itheta) / (omega_nkm * omega_nkm * omega_nkm);
        printf("J_dot.c: contrib_Re_outer = %.15e \n", contrib_Re_outer);
        
        #if 1
        double termIm1_outer = termR1_outer * Itheta;
        double termIm2_outer = -termR2_outer * Rtheta;
        double contrib_Im_outer = (termIm1_outer + termIm2_outer) / (omega_nkm * omega_nkm * omega_nkm);
        #endif

        // ------------------------------------------------
        // LOCAL accumulation (thread-safe)
        // ------------------------------------------------
        J_local_inner[0] += -i_n_inner * mu_outer * contrib_Re_inner;
        J_local_inner[1] += -i_k_inner * mu_outer * contrib_Re_inner;
        J_local_inner[2] += -i_m_inner * mu_outer * contrib_Re_inner;

        J_local_outer[0] += -i_n_outer * mu_inner * contrib_Re_outer;
        J_local_outer[1] += -i_k_outer * mu_inner * contrib_Re_outer;
        J_local_outer[2] += -i_m_outer * mu_inner * contrib_Re_outer;
        
        #if 0
        /* Check that this real part agrees with previous result */
        JRe_local_inner[0] += -i_n_inner * mu_outer * contrib_Re_inner;
        JRe_local_inner[1] += -i_k_inner * mu_outer * contrib_Re_inner;
        JRe_local_inner[2] += -i_m_inner * mu_outer * contrib_Re_inner;

        /* Check that the imaginary part is zero due to selection rules */
        JIm_local_inner[0] += -i_n_inner * mu_outer * contrib_Im_inner;
        JIm_local_inner[1] += -i_k_inner * mu_outer * contrib_Im_inner;
        JIm_local_inner[2] += -i_m_inner * mu_outer * contrib_Im_inner;
        #endif
        #if 1
        /* Check that this real part agrees with previous result */
        JRe_local_outer[0] += -i_n_outer * mu_inner * contrib_Re_outer;
        JRe_local_outer[1] += -i_k_outer * mu_inner * contrib_Re_outer;
        JRe_local_outer[2] += -i_m_outer * mu_inner * contrib_Re_outer;

        /* Check that the imaginary part is zero due to selection rules */
        JIm_local_outer[0] += -i_n_outer * mu_inner * contrib_Im_outer;
        JIm_local_outer[1] += -i_k_outer * mu_inner * contrib_Im_outer;
        JIm_local_outer[2] += -i_m_outer * mu_inner * contrib_Im_outer;
        #endif
      }
    }

    double thread_t1 = omp_get_wtime(); // End time of thread calculation

// --------------------------------------------------------
// Merge thread-local results
// --------------------------------------------------------
#pragma omp critical
    {
      // printf("Thread %d total time = %.6f s\n", tid,
      //        thread_t1 - thread_t0); // Thread ID and how long calculation took

      J_global_inner[0] += J_local_inner[0];
      J_global_inner[1] += J_local_inner[1];
      J_global_inner[2] += J_local_inner[2];

      J_global_outer[0] += J_local_outer[0];
      J_global_outer[1] += J_local_outer[1];
      J_global_outer[2] += J_local_outer[2];

      #if 0
      JRe_global_inner[0] += JRe_local_inner[0];
      JRe_global_inner[1] += JRe_local_inner[1];
      JRe_global_inner[2] += JRe_local_inner[2];

      JIm_global_inner[0] += JIm_local_inner[0];
      JIm_global_inner[1] += JIm_local_inner[1];
      JIm_global_inner[2] += JIm_local_inner[2];
      #endif
      #if 1
      JRe_global_outer[0] += JRe_local_outer[0];
      JRe_global_outer[1] += JRe_local_outer[1];
      JRe_global_outer[2] += JRe_local_outer[2];

      JIm_global_outer[0] += JIm_local_outer[0];
      JIm_global_outer[1] += JIm_local_outer[1];
      JIm_global_outer[2] += JIm_local_outer[2];
      #endif
    }

    free((char *)E0_inner);
    free((char *)E0_outer);
  }
#endif

  /* Modes i and -i give the same real contribution (omega -> -omega, Z -> (-1)^l conj(Z), cos(i theta) even,
   * sin(i theta) odd), so only i = 1..N_res is computed and counted twice. Their imaginary parts cancel
   * exactly, so the imaginary-part check is not accumulated here. */
  /* Step 1: every CKerr_RadialFunc call is its own task, t = 2*(i-1) + body (body 0 = inner, 1 = outer).
   * Each task keeps the E0/E1/C0 layout of the original per-thread buffer (10*nl doubles). */
  int n_tasks = 2 * N_res;
  double *mode_buf = (double *)malloc((size_t)(n_tasks * nl * 10) * sizeof(double));
  double *omegagw_buf = (double *)malloc((size_t)n_tasks * sizeof(double));
  if (!mode_buf || !omegagw_buf) {
    printf("Memory allocation failed in J_dot_tidal_openmp\n");
    exit(1);
  }

#pragma omp parallel for schedule(dynamic, 1)
  for (int t = 0; t < n_tasks; t++) {
    int i = t / 2 + 1, body = t % 2;
    if (tidal_mode_skipped(i, n_res_inner, k_res_inner, m_res_inner, n_res_outer, k_res_outer, m_res_outer)) continue;

    double *E0 = mode_buf + (size_t)t * nl * 10;
    double *C0 = E0 + 2 * nl;
    if (body == 0)
      CKerr_RadialFunc(Minv_inner, xuorig_inner, M, astar, i * n_res_inner, i * k_res_inner,
                       i * m_res_inner, 42, 42, nl, C0, &omegagw_buf[t], E0);
    else
      CKerr_RadialFunc(Minv_outer, xuorig_outer, M, astar, i * n_res_outer, i * k_res_outer,
                       i * m_res_outer, 42, 42, nl, C0, &omegagw_buf[t], E0);
  }

  /* Step 2: sum over (i, l) pairs */
#pragma omp parallel
  {
    double J_local_inner[3] = {0.0, 0.0, 0.0};
    double J_local_outer[3] = {0.0, 0.0, 0.0};
    double JRe_local_outer[3] = {0.0, 0.0, 0.0};

#pragma omp for schedule(dynamic, 1)
    for (int p = 0; p < N_res * nl; p++) {
      int i = p / nl + 1, il = p % nl;
      if (tidal_mode_skipped(i, n_res_inner, k_res_inner, m_res_inner, n_res_outer, k_res_outer, m_res_outer)) continue;

      int i_n_inner = i * n_res_inner;
      int i_k_inner = i * k_res_inner;
      int i_m_inner = i * m_res_inner;

      int i_n_outer = i * n_res_outer;
      int i_k_outer = i * k_res_outer;
      int i_m_outer = i * m_res_outer;

      double Rtheta = cos(i * theta_res_F);
      double Itheta = sin(i * theta_res_F);

      double *E0_inner = mode_buf + (size_t)(2 * (i - 1)) * nl * 10, *C0_inner = E0_inner + 2 * nl;
      double *E0_outer = mode_buf + (size_t)(2 * (i - 1) + 1) * nl * 10, *C0_outer = E0_outer + 2 * nl;
      double omegagw_inner = omegagw_buf[2 * (i - 1)], omegagw_outer = omegagw_buf[2 * (i - 1) + 1];
      double cscat[16], aux[4];

      // Set the frequency to the inner body value,
      // but either body is fine since the system is in resonance
      double omega_nkm = omegagw_inner;

      /* Frequency check, once per mode */
      if (il == 0) {
        double denom = fabs(omegagw_inner) + fabs(omegagw_outer);
        double rel_diff = denom > 0.0 ? 2.0 * fabs(omegagw_inner - omegagw_outer) / denom : 0.0;
        double delta_omega_tol = 1.0e-4; // 0.01%
        if (rel_diff > delta_omega_tol) {
          printf("WARNING: GW frequencies for inner and outer bodies differ by "
                 "%.6e (tol %.6e)\n",
                 rel_diff, delta_omega_tol);
        }
      }

      CKerr_GWScatMatrix(M, astar, omega_nkm, i_m_inner, E0_outer[il], cscat,
                         aux);

      double lambda = E0_inner[il] - 2 * astar * M * i_m_inner * omega_nkm +
                      astar * astar * M * M * omega_nkm * omega_nkm - 2;

      double P =
          omega_nkm - i_m_inner * astar * M /
                          (2 * M * (M + sqrt(M * M - astar * astar * M * M)));

      double numer = 256 * pow(2 * M * rH, 5) * P *
                     (P * P + 4 * epsilon * epsilon) *
                     (P * P + 16 * epsilon * epsilon) * omega_nkm *
                     omega_nkm * omega_nkm;

      double C2 =
          ((lambda + 2) * (lambda + 2) +
           4 * i_m_inner * astar * M * omega_nkm -
           4 * astar * astar * M * M * omega_nkm * omega_nkm) *
              (lambda * lambda + 36 * i_m_inner * astar * M * omega_nkm -
               36 * astar * astar * M * M * omega_nkm * omega_nkm) +
          (2 * lambda + 3) *
              (96 * astar * astar * M * M * omega_nkm * omega_nkm -
               48 * i_m_inner * astar * M * omega_nkm) +
          144 * omega_nkm * omega_nkm * (M * M - astar * astar * M * M);

      double alphankm = numer / C2;

      /* J_dot_tidal_inner */
      double termR1_inner =
            C0_inner[4 * il + 2] * cscat[0] * C0_outer[4 * il] +
            C0_outer[4 * il] * cscat[1] * C0_inner[4 * il + 3] +
            C0_outer[4 * il + 1] * C0_inner[4 * il + 3] * cscat[0] -
            C0_outer[4 * il + 1] * cscat[1] * C0_inner[4 * il + 2] +
            alphankm * (C0_outer[4 * il] * C0_inner[4 * il] +
              C0_outer[4 * il + 1] * C0_inner[4 * il + 1]);

      double termR2_inner =
            C0_outer[4 * il + 1] * cscat[0] * C0_inner[4 * il + 2] +
            C0_outer[4 * il + 1] * cscat[1] * C0_inner[4 * il + 3] -
            C0_outer[4 * il] * C0_inner[4 * il + 3] * cscat[0] +
            C0_outer[4 * il] * cscat[1] * C0_inner[4 * il + 2] +
            alphankm * (C0_outer[4 * il + 1] * C0_inner[4 * il] -
            C0_outer[4 * il] * C0_inner[4 * il + 1]);

      double contrib_Re_inner = (termR1_inner * Rtheta + termR2_inner * Itheta) / (omega_nkm * omega_nkm * omega_nkm);

      /* J_dot_tidal_outer */
      double termR1_outer =
      C0_inner[4 * il + 2] * C0_outer[4 * il + 2] +
      C0_inner[4 * il + 3] * C0_outer[4 * il + 3] -
      C0_inner[4 * il + 2] * cscat[0] * C0_outer[4 * il] +
      C0_inner[4 * il + 2] * cscat[1] * C0_outer[4 * il + 1] -
      C0_inner[4 * il + 3] * cscat[1] * C0_outer[4 * il] -
      C0_inner[4 * il + 3] * cscat[0] * C0_outer[4 * il + 1];

      double termR2_outer =
      C0_inner[4 * il + 2] * C0_outer[4 * il + 3] -
      C0_inner[4 * il + 3] * C0_outer[4 * il + 2] +
      C0_inner[4 * il + 3] * cscat[0] * C0_outer[4 * il] -
      C0_inner[4 * il + 3] * cscat[1] * C0_outer[4 * il + 1] -
      C0_inner[4 * il + 2] * cscat[1] * C0_outer[4 * il] -
      C0_inner[4 * il + 2] * cscat[0] * C0_outer[4 * il + 1];

      double contrib_Re_outer = (termR1_outer * Rtheta + termR2_outer * Itheta) / (omega_nkm * omega_nkm * omega_nkm);

      // ------------------------------------------------
      // LOCAL accumulation (thread-safe); factor 2 counts mode -i
      // ------------------------------------------------
      J_local_inner[0] += 2. * -i_n_inner * mu_outer * contrib_Re_inner;
      J_local_inner[1] += 2. * -i_k_inner * mu_outer * contrib_Re_inner;
      J_local_inner[2] += 2. * -i_m_inner * mu_outer * contrib_Re_inner;

      J_local_outer[0] += 2. * -i_n_outer * mu_inner * contrib_Re_outer;
      J_local_outer[1] += 2. * -i_k_outer * mu_inner * contrib_Re_outer;
      J_local_outer[2] += 2. * -i_m_outer * mu_inner * contrib_Re_outer;

      /* Check that this real part agrees with previous result */
      JRe_local_outer[0] += 2. * -i_n_outer * mu_inner * contrib_Re_outer;
      JRe_local_outer[1] += 2. * -i_k_outer * mu_inner * contrib_Re_outer;
      JRe_local_outer[2] += 2. * -i_m_outer * mu_inner * contrib_Re_outer;
    }

// --------------------------------------------------------
// Merge thread-local results
// --------------------------------------------------------
#pragma omp critical
    {
      J_global_inner[0] += J_local_inner[0];
      J_global_inner[1] += J_local_inner[1];
      J_global_inner[2] += J_local_inner[2];

      J_global_outer[0] += J_local_outer[0];
      J_global_outer[1] += J_local_outer[1];
      J_global_outer[2] += J_local_outer[2];

      JRe_global_outer[0] += JRe_local_outer[0];
      JRe_global_outer[1] += JRe_local_outer[1];
      JRe_global_outer[2] += JRe_local_outer[2];
    }
  }

  free(mode_buf);
  free(omegagw_buf);

  // ------------------------------------------------------------
  // Final output
  // ------------------------------------------------------------
  J_dot_td_inner[0] = J_global_inner[0];
  J_dot_td_inner[1] = J_global_inner[1];
  J_dot_td_inner[2] = J_global_inner[2];

  J_dot_td_outer[0] = J_global_outer[0];
  J_dot_td_outer[1] = J_global_outer[1];
  J_dot_td_outer[2] = J_global_outer[2];

  #if 0
  printf("J_dot.c: Real part of tidal resonance (inner) \n Jdotr Jdottheta Jdotphi: %.15e %.15e %.15e \n", JRe_global_inner[0], JRe_global_inner[1], JRe_global_inner[2]);
  printf("J_dot.c: Imag. part of tidal resonance (inner) \n Jdotr Jdottheta Jdotphi: %.15e %.15e %.15e \n", JIm_global_inner[0], JIm_global_inner[1], JIm_global_inner[2]);
  printf("J_dot.c: (calc(J_dot_td_inner) - J_dot_td_Re_inner) / calc(J_dot_td_inner): %.15e %.15e %.15e \n", 
    (JRe_global_inner[0] - J_global_inner[0]) / J_global_inner[0], (JRe_global_inner[1] - J_global_inner[1]) / J_global_inner[1], (JRe_global_inner[2] - J_global_inner[2]) / J_global_inner[2]);
  #endif
  #if 0
  printf("J_dot.c: Real part of tidal resonance (outer) \n Jdotr Jdottheta Jdotphi: %.15e %.15e %.15e \n", JRe_global_outer[0], JRe_global_outer[1], JRe_global_outer[2]);
  printf("J_dot.c: Imag. part of tidal resonance (outer): not computed, it cancels exactly between modes i and -i \n");
  printf("J_dot.c: (calc(J_dot_td_outer) - J_dot_td_Re_outer) / calc(J_dot_td_outer): %.15e %.15e %.15e \n", 
    (JRe_global_outer[0] - J_global_outer[0]) / J_global_outer[0], (JRe_global_outer[1] - J_global_outer[1]) / J_global_outer[1], (JRe_global_outer[2] - J_global_outer[2]) / J_global_outer[2]);
  #endif
}

/*
------------------------------------------------------------
Tidal resonance amplitudes: the Theta- and Gamma-independent part of the tidal J_dot and of Delta J
------------------------------------------------------------
The tidal J_dot of Eq. (36) of arXiv:2207.07733 is a Fourier series in the resonant angle Theta, one harmonic s per
multiple s*(n,k,m) of the resonance. J_dot_tidal_openmp (above) evaluates that series at one Theta. This routine does
the same Teukolsky work (same CKerr_RadialFunc calls, same expressions for alpha, lambda, termR1, termR2) but stops
before the angle is applied and returns, for each harmonic s = 1..N_res, the four real numbers

  amp[4*(s-1) + 0] = R1_inner = sum_l termR1_inner / omega_s^3      amp[4*(s-1) + 2] = R1_outer
  amp[4*(s-1) + 1] = R2_inner = sum_l termR2_inner / omega_s^3      amp[4*(s-1) + 3] = R2_outer

so that (J_dot_tidal_from_amplitudes, below)
  J_dot_td_inner,j(Theta) = 2 sum_s (-s N_inner,j) mu_outer (R1_inner cos(s Theta) + R2_inner sin(s Theta))
  J_dot_td_outer,j(Theta) = 2 sum_s (-s N_outer,j) mu_inner (R1_outer cos(s Theta) + R2_outer sin(s Theta))
and the jump across the resonance for any Gamma follows without new Teukolsky work (Delta_J_from_amplitudes, Delta_J.c).

The factor 2 counts harmonic -s. Under (n,k,m) -> -(n,k,m): omega -> -omega, alpha and Re c13 even, Im c13 odd, and
the amplitudes Z -> (+/-)(-1)^k conj(Z), so the product of an inner and an outer amplitude picks up (-1)^{s (k_inner + k_outer)}.
The -s term therefore EQUALS the +s term when s (k_inner - k_outer) is even, and CANCELS it when odd: the selection rule
k_inner - k_outer = even, Eq. (34) of arXiv:2207.07733. Harmonics with odd s (k_inner - k_outer) are set to zero here
(the full sum over +/-s gives zero for them to ~1e-13 of a single term; checked 2026-11-17 against Delta_J_tidal2).
Resonances written by the resonance finder always have even k_inner - k_outer.

Differences from J_dot_tidal_openmp, besides the output:
  - Nr, Nt (orbit grid of CKerr_RadialFunc) are arguments; J_dot_tidal_openmp has 42 x 42 written in;
  - harmonics are skipped when 2|s n| >= Nr or 2|s k| >= Nt for either body (the grid cannot resolve them). For
    Nr = Nt = 42 this is the same rule as J_dot_tidal_openmp's |s n|, |s k| > 20;
  - nothing is printed to stdout and exit() is never called: warnings go to stderr and failures are a return value,
    because pipeline executables print their result as one line on stdout;
  - the (s, l) terms are stored and added in a fixed order, so the result does not depend on the number of threads;
  - the mode functions come from the orbit cache (one CKerr_OrbitCacheBuild per body, CKerr_RadialFunc_cached per
    harmonic; radialfunc_cached.c). Build with -DTIDAL_AMP_ORIGINAL (make TIDAL_AMP=original) to call
    CKerr_RadialFunc for every harmonic instead: same numbers to ~1e-13, about ten times the CPU;
  - the l range of a harmonic can be extended (l_extend, next paragraph).

The l range. CKerr_RadialFunc returns l = l0, l0+1, ... with l0 = max(2, s|m|). J_dot_tidal_openmp and Delta_J_tidal2 take
the first nl of them for every harmonic. But harmonic s of a body with polar harmonic k only couples to multipoles
l >= s|m + k| (in the Newtonian limit the orbit enters (l, m) through the orbital-plane harmonic m' = m + k, and |m'| <= l;
checked numerically 2026-11-18: terms below that l are <= 7e-6 of the leading term at large radii, <= 8e-3 for
pericentres of a few M). So the first l that matters for the pair of bodies is
    l_th(s) = max(2, s|m|, s|m + k_inner|, s|m + k_outer|) [For a prograde orbit (L>0) ONLY! Must use s|m-k| for retrograde orbits (L<0)],
and when l_th(s) > l0 + nl - 1 the fixed range misses the harmonic altogether.
  l_extend = 0: l = l0 .. l0 + nl - 1 for every harmonic (the range of J_dot_tidal_openmp and Delta_J_tidal2);
  l_extend != 0: l = l0 .. l_th(s) + nl - 1, i.e. nl values of l counted from l_th(s), plus everything below. Terms are
                 only ever added, so a harmonic with l_th(s) = l0 is unchanged.
J_dot_tidal_nl_harmonic gives the number of l values of harmonic s, J_dot_tidal_nl_max the largest over s = 1..N_res.

Inputs: nl = number of l per harmonic (see l_extend), N_res = number of harmonics, mode vectors (n,k,m) of both bodies
        (m must be equal), (ra, rp, I) of both bodies, M, astar, Nr, Nt.
Outputs: amp[4*N_res] (zeros for skipped harmonics);
         optional (NULL to skip) amp_l[4*N_res*nl_max], nl_max = J_dot_tidal_nl_max(...): the four numbers for each
         (s, l), index 4*((s-1)*nl_max + il) + 0..3, l = max(2, s|m|) + il (zero for il beyond the harmonic's range);
         optional omegagw[2*N_res]: s N.Omega of the inner body, then of the outer body, as CKerr_RadialFunc returns them;
         optional *n_used: number of harmonics that were summed.
Returns 0, or 1 (bad arguments), 2 (no orbit found for the given ra, rp, I), 3 (out of memory).
*/

/* Y table of one orbit cache (radialfunc_cached.c) is 16 * Nr * Nt * ((lc+1)^2 - 4) bytes. If the lc that the summed
 * harmonics need would exceed this many bytes per body, lc is lowered to fit and the modes beyond it are evaluated
 * without the table (CKerr_SpheroidalYSLM inside the cached routine: same numbers, slower). 42 x 42 never gets there
 * (lc 133); 168 x 168 does at lc 32. */
#ifndef TIDAL_CACHE_MAX_BYTES
#define TIDAL_CACHE_MAX_BYTES 5.0e8
#endif

/* Number of l values summed for harmonic s: nl, plus (when l_extend != 0) the distance from l0 = max(2, s|m|) to
 * l_th(s) = max(2, s|m|, s|m + k_inner|, s|m + k_outer|). See the comment above. */
int J_dot_tidal_nl_harmonic(int nl, int l_extend, int s, int k_res_inner, int k_res_outer, int m_res) {
  int l0 = abs(s * m_res) < 2 ? 2 : abs(s * m_res), l_th = l0;
  if (!l_extend) return nl;
  if (abs(s * (m_res + k_res_inner)) > l_th) l_th = abs(s * (m_res + k_res_inner));
  if (abs(s * (m_res + k_res_outer)) > l_th) l_th = abs(s * (m_res + k_res_outer));
  return nl + (l_th - l0);
}

/* Largest J_dot_tidal_nl_harmonic over s = 1..N_res: the l stride of amp_l */
int J_dot_tidal_nl_max(int nl, int l_extend, int N_res, int k_res_inner, int k_res_outer, int m_res) {
  int s, nl_max = nl;
  for (s = 1; s <= N_res; s++) {
    int nl_s = J_dot_tidal_nl_harmonic(nl, l_extend, s, k_res_inner, k_res_outer, m_res);
    if (nl_s > nl_max) nl_max = nl_s;
  }
  return nl_max;
}

/* 0 = harmonic s is summed; 1 = skipped (zero mode, or not resolved by the Nr x Nt grid);
 * 2 = zero by the selection rule, s (k_inner - k_outer) odd */
static int tidal_harmonic_status(int s, int n_res_inner, int k_res_inner, int m_res_inner,
                                 int n_res_outer, int k_res_outer, int m_res_outer, int Nr, int Nt) {
  int s_n_inner = s * n_res_inner, s_k_inner = s * k_res_inner, s_m_inner = s * m_res_inner;
  int s_n_outer = s * n_res_outer, s_k_outer = s * k_res_outer, s_m_outer = s * m_res_outer;
  if ((s_n_outer == 0 && s_k_outer == 0 && s_m_outer == 0) ||
      (s_n_inner == 0 && s_k_inner == 0 && s_m_inner == 0) ||
      2 * abs(s_n_outer) >= Nr || 2 * abs(s_k_outer) >= Nt || 2 * abs(s_n_inner) >= Nr ||
      2 * abs(s_k_inner) >= Nt)
    return 1;
  if (abs(s_k_inner - s_k_outer) % 2 == 1) return 2;
  return 0;
}

int J_dot_tidal_amplitudes_openmp(int nl, int l_extend, int N_res, int n_res_inner, int n_res_outer,
                                  int k_res_inner, int k_res_outer, int m_res_inner,
                                  int m_res_outer, double ra_inner, double rp_inner,
                                  double I_inner, double ra_outer, double rp_outer,
                                  double I_outer, double M, double astar, int Nr, int Nt,
                                  double *amp, double *amp_l, double *omegagw, int *n_used) {
  if (n_used) *n_used = 0;
  if (nl < 1 || N_res < 1 || Nr < 1 || Nt < 1 || !amp) {
    fprintf(stderr, "J_dot.c: J_dot_tidal_amplitudes_openmp needs nl, N_res, Nr, Nt >= 1\n");
    return 1;
  }
  /* l stride of the buffers: the largest number of l values any harmonic sums (= nl when l_extend is 0) */
  int nl_max = J_dot_tidal_nl_max(nl, l_extend, N_res, k_res_inner, k_res_outer, m_res_inner);
  for (int j = 0; j < 4 * N_res; j++) amp[j] = 0.0;
  if (amp_l) for (int j = 0; j < 4 * N_res * nl_max; j++) amp_l[j] = 0.0;
  if (omegagw) for (int j = 0; j < 2 * N_res; j++) omegagw[j] = 0.0;
  if (m_res_inner != m_res_outer) {
    /* CKerr_RadialFunc starts its l list at max(2,|m|), so with different m the same il would be a different l */
    fprintf(stderr, "J_dot.c: J_dot_tidal_amplitudes_openmp needs m_res_inner = m_res_outer (got %d, %d)\n", m_res_inner, m_res_outer);
    return 1;
  }

  // --------------------------------------------------------
  // Orbit quantities (as in J_dot_tidal_openmp)
  // --------------------------------------------------------
  double EQL_inner[3], J_inner[3], Minv_inner[9];
  double EQL_outer[3], J_outer[3], Minv_outer[9];
  double xuorig_inner[6], xuorig_outer[6];

  /* ra_rp_I2EQL can return without setting EQL (no real L/E solution), so start from NaN */
  for (int i = 0; i < 3; i++) EQL_inner[i] = EQL_outer[i] = NAN;

  ra_rp_I2EQL(ra_inner, EQL_inner, rp_inner, I_inner, astar, M);
  ra_rp_I2EQL(ra_outer, EQL_outer, rp_outer, I_outer, astar, M);
  for (int i = 0; i < 3; i++) {
    if (!isfinite(EQL_inner[i]) || !isfinite(EQL_outer[i])) {
      fprintf(stderr, "J_dot.c: J_dot_tidal_amplitudes_openmp found no EQL for the %s orbit\n", isfinite(EQL_inner[i]) ? "outer" : "inner");
      return 2;
    }
  }

  CKerr_EQL2J(EQL_inner, J_inner, M, astar, NULL);
  CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
  for (int i = 0; i < 3; i++) {
    if (!isfinite(J_inner[i]) || !isfinite(J_outer[i])) {
      fprintf(stderr, "J_dot.c: J_dot_tidal_amplitudes_openmp found no J for the %s orbit\n", isfinite(J_inner[i]) ? "outer" : "inner");
      return 2;
    }
  }

  CKerr_Minverse(J_inner, Minv_inner, M, astar);
  CKerr_Minverse(J_outer, Minv_outer, M, astar);

  // Set up initial conditions to integrate on tori
  CKerr_TorusOrigin(J_inner, xuorig_inner, M, astar);
  CKerr_TorusOrigin(J_outer, xuorig_outer, M, astar);

  double rH = M + sqrt(M * M - astar * astar * M * M);
  double epsilon = sqrt(M * M - astar * astar * M * M) / (4 * M * rH);

  /* Step 1: every mode-function call is its own task, t = 2*(s-1) + body (body 0 = inner, 1 = outer).
   * Each task has the E0/E1/C0 layout of J_dot_tidal_openmp with nl_max in place of nl (10*nl_max doubles), of which
   * harmonic s fills its first nl_s = J_dot_tidal_nl_harmonic(...) values of l. calloc: the rest stays zero. */
  int n_tasks = 2 * N_res;
  double *mode_buf = (double *)calloc((size_t)n_tasks * nl_max * 10, sizeof(double));
  double *omegagw_buf = (double *)calloc((size_t)n_tasks, sizeof(double));
  double *term_buf = amp_l ? amp_l : (double *)calloc((size_t)4 * N_res * nl_max, sizeof(double));
  if (!mode_buf || !omegagw_buf || !term_buf) {
    fprintf(stderr, "J_dot.c: memory allocation failed in J_dot_tidal_amplitudes_openmp\n");
    free(mode_buf);
    free(omegagw_buf);
    if (!amp_l) free(term_buf);
    return 3;
  }

#ifndef TIDAL_AMP_ORIGINAL
  /* Orbit caches: torus samples and Y^{-2}_{lm} of each body, computed once and shared by its harmonics. lc is the
   * largest l any summed harmonic reads from the Y table (CKerr_RadialFunc_cached_lmax), lowered if the table would
   * exceed TIDAL_CACHE_MAX_BYTES. Nothing is built when no harmonic is summed. */
  CKerr_OrbitCache cache_inner, cache_outer;
  long lc_inner = 0, lc_outer = 0;
  for (int s = 1; s <= N_res; s++) {
    if (tidal_harmonic_status(s, n_res_inner, k_res_inner, m_res_inner, n_res_outer, k_res_outer, m_res_outer, Nr, Nt)) continue;
    int nl_s = J_dot_tidal_nl_harmonic(nl, l_extend, s, k_res_inner, k_res_outer, m_res_inner);
    long l_in = CKerr_RadialFunc_cached_lmax(Minv_inner, M, astar, s * n_res_inner, s * k_res_inner, s * m_res_inner, nl_s);
    long l_out = CKerr_RadialFunc_cached_lmax(Minv_outer, M, astar, s * n_res_outer, s * k_res_outer, s * m_res_outer, nl_s);
    if (l_in > lc_inner) lc_inner = l_in;
    if (l_out > lc_outer) lc_outer = l_out;
  }
  int have_cache = (lc_inner > 0);
  if (have_cache) {
    long lc_cap = (long)sqrt(TIDAL_CACHE_MAX_BYTES / (16.0 * (double)Nr * (double)Nt) + 4.0) - 1;
    if (lc_cap < 2) lc_cap = 2;
    if (lc_inner > lc_cap || lc_outer > lc_cap) {
      fprintf(stderr, "J_dot.c: note: Y table limited to l <= %ld (needed %ld, %ld) to stay under %.0f MB per body on the %d x %d grid; "
              "higher l are computed without the table (slower, same numbers)\n",
              lc_cap, lc_inner, lc_outer, TIDAL_CACHE_MAX_BYTES / 1.0e6, Nr, Nt);
      if (lc_inner > lc_cap) lc_inner = lc_cap;
      if (lc_outer > lc_cap) lc_outer = lc_cap;
    }
    CKerr_OrbitCacheBuild(Minv_inner, xuorig_inner, M, astar, Nr, Nt, lc_inner, &cache_inner);
    CKerr_OrbitCacheBuild(Minv_outer, xuorig_outer, M, astar, Nr, Nt, lc_outer, &cache_outer);
  }
#endif

#pragma omp parallel for schedule(dynamic, 1)
  for (int t = 0; t < n_tasks; t++) {
    int s = t / 2 + 1, body = t % 2;
    if (tidal_harmonic_status(s, n_res_inner, k_res_inner, m_res_inner, n_res_outer, k_res_outer, m_res_outer, Nr, Nt)) continue;

    int nl_s = J_dot_tidal_nl_harmonic(nl, l_extend, s, k_res_inner, k_res_outer, m_res_inner);
    double *E0 = mode_buf + (size_t)t * nl_max * 10;
    double *C0 = E0 + 2 * nl_max;
#ifdef TIDAL_AMP_ORIGINAL
    if (body == 0)
      CKerr_RadialFunc(Minv_inner, xuorig_inner, M, astar, s * n_res_inner, s * k_res_inner,
                       s * m_res_inner, Nr, Nt, nl_s, C0, &omegagw_buf[t], E0);
    else
      CKerr_RadialFunc(Minv_outer, xuorig_outer, M, astar, s * n_res_outer, s * k_res_outer,
                       s * m_res_outer, Nr, Nt, nl_s, C0, &omegagw_buf[t], E0);
#else
    if (body == 0)
      CKerr_RadialFunc_cached(Minv_inner, &cache_inner, M, astar, s * n_res_inner, s * k_res_inner,
                              s * m_res_inner, nl_s, C0, &omegagw_buf[t], E0);
    else
      CKerr_RadialFunc_cached(Minv_outer, &cache_outer, M, astar, s * n_res_outer, s * k_res_outer,
                              s * m_res_outer, nl_s, C0, &omegagw_buf[t], E0);
#endif
  }
#ifndef TIDAL_AMP_ORIGINAL
  if (have_cache) {
    CKerr_OrbitCacheFree(&cache_inner);
    CKerr_OrbitCacheFree(&cache_outer);
  }
#endif

  /* Step 2: one (s, l) term per iteration, each written to its own slot of term_buf.
   * The expressions are those of step 2 of J_dot_tidal_openmp, without the cos(s Theta), sin(s Theta) factors. */
#pragma omp parallel for schedule(dynamic, 1)
  for (int p = 0; p < N_res * nl_max; p++) {
    int s = p / nl_max + 1, il = p % nl_max;
    if (tidal_harmonic_status(s, n_res_inner, k_res_inner, m_res_inner, n_res_outer, k_res_outer, m_res_outer, Nr, Nt)) continue;
    if (il >= J_dot_tidal_nl_harmonic(nl, l_extend, s, k_res_inner, k_res_outer, m_res_inner)) continue;

    int i_m_inner = s * m_res_inner;

    double *E0_inner = mode_buf + (size_t)(2 * (s - 1)) * nl_max * 10, *C0_inner = E0_inner + 2 * nl_max;
    double *E0_outer = mode_buf + (size_t)(2 * (s - 1) + 1) * nl_max * 10, *C0_outer = E0_outer + 2 * nl_max;
    double cscat[16], aux[4];

    // Set the frequency to the inner body value,
    // but either body is fine since the system is in resonance
    double omega_nkm = omegagw_buf[2 * (s - 1)];
    
    // One scattering matrix c_ij(l, m, omega) for both bodies (Hirata II Eqs. (94) and (98) use the same c13);
    // eigenvalue from the outer body's mode (equal to the inner one's at resonance)
    CKerr_GWScatMatrix(M, astar, omega_nkm, i_m_inner, E0_outer[il], cscat,
                       aux);

    double lambda = E0_inner[il] - 2 * astar * M * i_m_inner * omega_nkm +
                    astar * astar * M * M * omega_nkm * omega_nkm - 2;

    double P =
        omega_nkm - i_m_inner * astar * M /
                        (2 * M * (M + sqrt(M * M - astar * astar * M * M)));

    double numer = 256 * pow(2 * M * rH, 5) * P *
                   (P * P + 4 * epsilon * epsilon) *
                   (P * P + 16 * epsilon * epsilon) * omega_nkm *
                   omega_nkm * omega_nkm;

    double C2 =
        ((lambda + 2) * (lambda + 2) +
         4 * i_m_inner * astar * M * omega_nkm -
         4 * astar * astar * M * M * omega_nkm * omega_nkm) *
            (lambda * lambda + 36 * i_m_inner * astar * M * omega_nkm -
             36 * astar * astar * M * M * omega_nkm * omega_nkm) +
        (2 * lambda + 3) *
            (96 * astar * astar * M * M * omega_nkm * omega_nkm -
             48 * i_m_inner * astar * M * omega_nkm) +
        144 * omega_nkm * omega_nkm * (M * M - astar * astar * M * M);

    double alphankm = numer / C2;

    /* For INNER/OUTER Body: */
    /* C0_{inner/outer}[4*il    ] = Re[Z_down] */
    /* C0_{inner/outer}[4*il + 1] = Im[Z_down] */
    /* C0_{inner/outer}[4*il + 2] = Re[Z_out] */
    /* C0_{inner/outer}[4*il + 3] = Im[Z_out] */

    /* J_dot_tidal_inner (2207.07733 Eq. (35), Hirata II Eq. (94)), with cscat above */
    double termR1_inner =
          C0_inner[4 * il + 2] * cscat[0] * C0_outer[4 * il] +
          C0_outer[4 * il] * cscat[1] * C0_inner[4 * il + 3] +
          C0_outer[4 * il + 1] * C0_inner[4 * il + 3] * cscat[0] -
          C0_outer[4 * il + 1] * cscat[1] * C0_inner[4 * il + 2] +
          alphankm * (C0_outer[4 * il] * C0_inner[4 * il] +
            C0_outer[4 * il + 1] * C0_inner[4 * il + 1]);

    double termR2_inner =
          C0_outer[4 * il + 1] * cscat[0] * C0_inner[4 * il + 2] +
          C0_outer[4 * il + 1] * cscat[1] * C0_inner[4 * il + 3] -
          C0_outer[4 * il] * C0_inner[4 * il + 3] * cscat[0] +
          C0_outer[4 * il] * cscat[1] * C0_inner[4 * il + 2] +
          alphankm * (C0_outer[4 * il + 1] * C0_inner[4 * il] -
          C0_outer[4 * il] * C0_inner[4 * il + 1]);

    /* J_dot_tidal_outer (Hirata II Eq. (98)), with the same cscat */
    double termR1_outer =
    C0_inner[4 * il + 2] * C0_outer[4 * il + 2] +
    C0_inner[4 * il + 3] * C0_outer[4 * il + 3] -
    C0_inner[4 * il + 2] * cscat[0] * C0_outer[4 * il] +
    C0_inner[4 * il + 2] * cscat[1] * C0_outer[4 * il + 1] -
    C0_inner[4 * il + 3] * cscat[1] * C0_outer[4 * il] -
    C0_inner[4 * il + 3] * cscat[0] * C0_outer[4 * il + 1];

    double termR2_outer =
    C0_inner[4 * il + 2] * C0_outer[4 * il + 3] -
    C0_inner[4 * il + 3] * C0_outer[4 * il + 2] +
    C0_inner[4 * il + 3] * cscat[0] * C0_outer[4 * il] -
    C0_inner[4 * il + 3] * cscat[1] * C0_outer[4 * il + 1] -
    C0_inner[4 * il + 2] * cscat[1] * C0_outer[4 * il] -
    C0_inner[4 * il + 2] * cscat[0] * C0_outer[4 * il + 1];

    double omega_cubed = omega_nkm * omega_nkm * omega_nkm;
    double *term = term_buf + (size_t)4 * ((s - 1) * nl_max + il);
    term[0] = termR1_inner / omega_cubed;
    term[1] = termR2_inner / omega_cubed;
    term[2] = termR1_outer / omega_cubed;
    term[3] = termR2_outer / omega_cubed;
  }

  /* Serial part: add the l terms of each harmonic in a fixed order, and check the two frequencies */
  int used = 0;
  for (int s = 1; s <= N_res; s++) {
    int status = tidal_harmonic_status(s, n_res_inner, k_res_inner, m_res_inner, n_res_outer, k_res_outer, m_res_outer, Nr, Nt);
    if (status == 2)
      fprintf(stderr, "J_dot.c: note: harmonic s = %d has odd s (k_inner - k_outer); it vanishes by the selection rule and is set to zero\n", s);
    if (status) continue;
    used++;

    /* il beyond the harmonic's own range holds zeros, so the whole stride can be added */
    for (int il = 0; il < nl_max; il++)
      for (int c = 0; c < 4; c++) amp[4 * (s - 1) + c] += term_buf[(size_t)4 * ((s - 1) * nl_max + il) + c];

    double omegagw_inner = omegagw_buf[2 * (s - 1)], omegagw_outer = omegagw_buf[2 * (s - 1) + 1];
    if (omegagw) {
      omegagw[2 * (s - 1)] = omegagw_inner;
      omegagw[2 * (s - 1) + 1] = omegagw_outer;
    }

    /* Frequency check (same test and tolerance as J_dot_tidal_openmp, but on stderr). Both frequencies scale with s,
     * so the fractional difference is the same for every harmonic: it is tested, and reported, once per call. */
    if (used == 1) {
      double denom = fabs(omegagw_inner) + fabs(omegagw_outer);
      double rel_diff = denom > 0.0 ? 2.0 * fabs(omegagw_inner - omegagw_outer) / denom : 0.0;
      double delta_omega_tol = 1.0e-4; // 0.01%
      if (rel_diff > delta_omega_tol) {
        fprintf(stderr, "WARNING: GW frequencies for inner and outer bodies differ by "
                "%.6e (tol %.6e)\n",
                rel_diff, delta_omega_tol);
      }
    }
  }
  if (n_used) *n_used = used;

  free(mode_buf);
  free(omegagw_buf);
  if (!amp_l) free(term_buf);
  return 0;
}

/* Tidal J_dot of both bodies at resonant angle theta_res_F, from the amplitudes of J_dot_tidal_amplitudes_openmp.
 * Same quantity as J_dot_tidal_openmp returns: d(J/mu_body)/dt in physical time, inner ~ mu_outer, outer ~ mu_inner.
 * Cheap (no Teukolsky work), so it can be called for many angles from one set of amplitudes. */
void J_dot_tidal_from_amplitudes(int N_res, int n_res_inner, int n_res_outer, int k_res_inner,
                                 int k_res_outer, int m_res_inner, int m_res_outer,
                                 const double *amp, double theta_res_F, double mu_inner,
                                 double mu_outer, double *J_dot_td_inner, double *J_dot_td_outer) {
  for (int j = 0; j < 3; j++) J_dot_td_inner[j] = J_dot_td_outer[j] = 0.0;

  for (int s = 1; s <= N_res; s++) {
    double Rtheta = cos(s * theta_res_F);
    double Itheta = sin(s * theta_res_F);

    double contrib_Re_inner = amp[4 * (s - 1)] * Rtheta + amp[4 * (s - 1) + 1] * Itheta;
    double contrib_Re_outer = amp[4 * (s - 1) + 2] * Rtheta + amp[4 * (s - 1) + 3] * Itheta;

    /* factor 2 counts harmonic -s */
    J_dot_td_inner[0] += 2. * -(s * n_res_inner) * mu_outer * contrib_Re_inner;
    J_dot_td_inner[1] += 2. * -(s * k_res_inner) * mu_outer * contrib_Re_inner;
    J_dot_td_inner[2] += 2. * -(s * m_res_inner) * mu_outer * contrib_Re_inner;

    J_dot_td_outer[0] += 2. * -(s * n_res_outer) * mu_inner * contrib_Re_outer;
    J_dot_td_outer[1] += 2. * -(s * k_res_outer) * mu_inner * contrib_Re_outer;
    J_dot_td_outer[2] += 2. * -(s * m_res_outer) * mu_inner * contrib_Re_outer;
  }
}

/* Keplerian expression for evolution of L_z, while only keeping the resonant
 * prefactor terms. We will compare this to our general machine for computing
 * J_phi_dot since it equals L_z_dot */
double J_dot_phi_Kepler(double mu_outer, double r_outer, double apo,
                        double peri, double incline, double Theta_res) {
  // double resonant_angle;
  double x_orbitplane_square, y_orbitplane_square;
  double a, eccen;
  double torque_res;

  /* Semi-major axis and eccentricity as a function of apocenter and pericenter
   */
  a = (apo + peri) / 2;
  eccen = (apo - peri) / (apo + peri);
  // resonant_angle = 2 * (Omega + omega);

  /* x^2 and y^2 in the orbit plane after being averaged over mean anomly, M */
  x_orbitplane_square = a * a * (1 + 4 * eccen * eccen) / 2;
  y_orbitplane_square = a * a * (1 - eccen * eccen) / 2;

  // printf("Semi-major axis: %lg \n", a);
  // printf("Eccentricity: %lg \n", eccen);
  // printf("<x^2>_orbitplane = %lg \n", x_orbitplane_square);
  // printf("<y^2>_orbitplane = %lg \n", y_orbitplane_square);

  /* Expression for L_z_dot in the observer frame and only keeping the prefactor
   * of the term sin(-resonant_angle) */
  torque_res = -3 * mu_outer / (8 * r_outer * r_outer * r_outer) *
               (x_orbitplane_square - y_orbitplane_square) *
               (1 + cos(incline)) * (1 + cos(incline)) * sin(-2 * Theta_res);
  return (torque_res);
}
