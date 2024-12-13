#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CKerr.h"

/* Computing the J_dot due to self-force from arXiv:1905.00030v2 Eq. (12) */

void J_dot_selfforce(int nl, int nmax, int kmax, int mmax, double apo, double rp, double radius_outer, double I, double M, double astar, double *J_dot_sf){
	int i, i_n, i_k, i_m, il;
	double EQL[3], J[3], Minv[9], Omega[3], info[6], info_outer[6], xuorig[6], xuorig_outer[6];
	//int n_res_inner = 1, k_res_inner = 2, m_res_inner = -2;
	double Z_out_square, Z_down_square;
	double omega_nkm, omegagw;
	double rH;
	double J_dot_r_old, J_dot_theta_old, J_dot_phi_old;
	double epsilon, lambda, numer, C2, P, alphankm;
	double *C0, *E0, *E1;

	/* For CIRCULAR orbits */
	if(apo == 0 && rp == 0 && I == 0){
		/* Orbital frequencies for circular equitorial orbits */
		CKerr_FindEQL_IRCirc(0, radius_outer, EQL, M, astar);
		CKerr_EQL2J(EQL, J, M, astar, NULL);
		CKerr_Minverse(J, Minv, M, astar);
		CKerr_Minv2Omega(Minv, Omega);
		//printf("Minv for circular equitorial case = \n %lg %lg %lg \n", Minv[0], Minv[1], Minv[2]);
		//printf(" %lg %lg %lg \n", Minv[3], Minv[4], Minv[5]);
		//printf(" %lg %lg %lg \n", Minv[6], Minv[7], Minv[8]);
		//printf("Js are: %lg %lg %lg \n", J[0], J[1], J[2]);
		//printf("Omega for circular equatorial is: %lg %lg %lg \n", Omega[0], Omega[1], Omega[2]);

		/* Torus Origin for circular equatorial orbits */
		CKerr_TorusOrigin(J, xuorig, M, astar);
	}

	/* For GENERIC orbits */
	else{
		/* Routine to find location of resonant apocenter */
		//guess1 = rp + 3.;
		//guess2 = radius_outer - 3.;
		//apo = find_resonance_apo(n_res_inner, k_res_inner, m_res_inner, radius_outer, guess1, guess2, rp, I, astar, M);
		
		/* Get orbital frequencies for inner generic orbit */
		ra_rp_I2EQL(apo, EQL, rp, I, astar, M);
		CKerr_EQL2J(EQL, J, M, astar, NULL);
		CKerr_Minverse(J, Minv, M, astar);
		CKerr_Minv2Omega(Minv, Omega);
		//printf("Minv for generic case = \n %lg %lg %lg \n", Minv[0], Minv[1], Minv[2]);
		//printf(" %lg %lg %lg \n", Minv[3], Minv[4], Minv[5]);
		//printf(" %lg %lg %lg \n", Minv[6], Minv[7], Minv[8]);
		//printf("Js are: %lg %lg %lg \n", J[0], J[1], J[2]);
		//printf("Omega for generic is: %lg %lg %lg \n", Omega[0], Omega[1], Omega[2]);

		/* Torus origin for inner body orbit */
		CKerr_TorusOrigin(J, xuorig, M, astar);
		//for(i=0;i<6;i++) printf("xuorig[%d] = %lg\n", i, xuorig[i]);
	}
	
	

  	/* This takes the M_inv, the origin in phase space, mass and spin of BH, {n,k,m}, 
	number of steps in r and theta, number of modes (nl), and the angular frequency of the perturber orbit
	Returns the Z-coeff. and energy eigenvalues E[nl]. This is also
	computed at the specific resonance condition (1, 2, -2). */
	E0 = (double*)malloc((size_t)(nl*10*sizeof(double)));
  	E1 = E0 + nl;
  	C0 = E1 + nl;
	
	//printf("IM HERE 3 \n");

	rH = M + sqrt(M*M - astar*M*astar*M);
	epsilon = sqrt(M*M - astar*M*astar*M) / (4 * M * rH);

	J_dot_sf[0] = 0;
	J_dot_sf[1] = 0; 
	J_dot_sf[2] = 0;
	//printf("IM HERE 4 \n");
	/* Compute the J_dots for the inner orbit */
	for (i_n = -nmax; i_n <= nmax; i_n++){
		for (i_k = -kmax; i_k <= kmax; i_k++){
			for (i_m = -mmax; i_m <= mmax; i_m++){
				if(i_n == 0 && i_k == 0 && i_m == 0)
					continue;
				if(astar == 0 && i_n == 0 && i_k == -i_m) //Schwarzschild case a == 0
					continue;
				CKerr_RadialFunc(Minv, xuorig, M, astar, i_n, i_k, i_m, 14, 14, nl, C0, &omegagw, E0);
				for (il = 0; (il) < nl; il++){
					
					omega_nkm = i_n * Omega[0] + i_k * Omega[1] + i_m * Omega[2];
					lambda = E0[il] - 2 * astar * M * i_m * omega_nkm + astar * M * astar * M * omega_nkm * omega_nkm - 2;
					P = omega_nkm - i_m * astar * M / (2 * M * rH);
					numer = 256 * pow(2 * M * rH,5) * P * (P*P + 4 * epsilon*epsilon) * (P*P + 16 * epsilon*epsilon) * omega_nkm * omega_nkm * omega_nkm;
					C2 = ((lambda + 2)*(lambda + 2) + 4 * i_m * astar * omega_nkm - 4 * astar * M * astar * M * omega_nkm*omega_nkm) * (lambda*lambda + 36 * i_m * astar * M *omega_nkm - 36 * astar * M * astar * M * omega_nkm*omega_nkm) + (2 * lambda +3) * (96 * astar * M * astar * M * omega_nkm*omega_nkm - 48 * i_m * astar * M * omega_nkm) + 144 * omega_nkm*omega_nkm * (M*M - astar*astar*M*M);
					alphankm = numer / C2;
					//printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg\n", i_n, i_k, i_k, il, omega_nkm, lambda, alphankm);


					Z_down_square = C0[4*il]*C0[4*il] + C0[4*il+1]*C0[4*il+1];
					Z_out_square = C0[4*il+2]*C0[4*il+2] + C0[4*il+3]*C0[4*il+3];
					//printf("%i \t %i \t %i \t %i \t %lg \t %lg \n", i_n, i_k, i_m, il, Z_down_square, Z_out_square);
					//printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n, i_k, i_m, il, C0[4*il], C0[4*il+1], C0[4*il+2], C0[4*il+3], Z_down_square, Z_out_square);
					
					J_dot_sf[0] += -i_n * (Z_out_square + alphankm * Z_down_square) / (2. * omega_nkm*omega_nkm*omega_nkm);
					J_dot_sf[1] += -i_k * (Z_out_square + alphankm * Z_down_square) / (2. * omega_nkm*omega_nkm*omega_nkm);
					J_dot_sf[2] += -i_m * (Z_out_square + alphankm * Z_down_square) / (2. * omega_nkm*omega_nkm*omega_nkm);

				}
			}
		}
	}
	free((char*)E0);
}

/* Computing the J_dot_tidal,i */
void J_dot_tidal(int nl, int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double ra_inner, double rp_inner, double radius_outer, double I_inner, double ra_outer, double rp_outer, double I_outer, double M, double astar, double theta_res_F, double mu_outer ,double *J_dot_td){
	int i, i_n_inner, i_k_inner, i_m_inner, il;
	int i_n_outer, i_k_outer, i_m_outer;
	double EQL_inner[3], J_inner[3], EQL_outer[3], J_outer[3], Minv_inner[9], Minv_outer[9], Omega_inner[3], info[6], info_outer[6], xuorig_inner[6], xuorig_outer[6], cscat[16], aux[4];
	double term, another_term, last_term;
	double Gamma, sgn_Gamma;
	double omega_nkm, omegagw_inner, omegagw_outer;
	double rH;
	double Rtheta, Itheta;
	double epsilon, lambda, numer, C2, P, alphankm;
	double *C0_inner, *E0_inner, *E1_inner, *C0_outer, *E0_outer, *E1_outer;


	/* Converting apocenter, pericenter, and inclination angle into frequencies for: */
	/* Inner body */
	ra_rp_I2EQL(ra_inner, EQL_inner, rp_inner, I_inner, astar, M);
	CKerr_EQL2J(EQL_inner, J_inner, M, astar, NULL);
	CKerr_Minverse(J_inner, Minv_inner, M, astar);
	CKerr_Minv2Omega(Minv_inner, Omega_inner);
	//printf("EQL_inner: %lg %lg %lg \n", EQL_inner[0], EQL_inner[1], EQL_inner[2]);
	//printf("Minv for inner body = \n %lg %lg %lg \n", Minv_inner[0], Minv_inner[1], Minv_inner[2]);
	//printf(" %lg %lg %lg \n", Minv_inner[3], Minv_inner[4], Minv_inner[5]);
	//printf(" %lg %lg %lg \n", Minv_inner[6], Minv_inner[7], Minv_inner[8]);

	//printf("dtau/dt = %lg \n", EQL_inner[0] - Omega_inner[0]*J_inner[0] - Omega_inner[1]*J_inner[1] - Omega_inner[2]*J_inner[2]);
	/* Outer Body */
	// CKerr_FindEQL_IRCirc(0, radius_outer, EQL_outer, M, astar);
	// CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
	// CKerr_Minverse(J_outer, Minv_outer, M, astar);
	// //printf("Minv for outer body = \n %lg %lg %lg \n", Minv_outer[0], Minv_outer[1], Minv_outer[2]);
	// //printf(" %lg %lg %lg \n", Minv_outer[3], Minv_outer[4], Minv_outer[5]);
	// //printf(" %lg %lg %lg \n", Minv_outer[6], Minv_outer[7], Minv_outer[8]);

	// /* Torus origin for inner body orbit */
	// CKerr_TorusOrigin(J_inner, xuorig_inner, M, astar);

	// /* Torus origin for outer body orbit */
	// CKerr_getData_CircEq(M, astar, radius_outer, info_outer);
	// xuorig_outer[0] = radius_outer; xuorig_outer[1] = M_PI/2.; xuorig_outer[2] = 0.;
  	// xuorig_outer[3] = 0.; xuorig_outer[4] = 0.;      xuorig_outer[5] = info_outer[0];

	/* Converting apocenter, pericenter, and inclination angle into frequencies for: */
	/* Outer body, both cases*/
	/* For CIRCULAR orbits */
	if(ra_outer == 0 && rp_outer == 0 && I_outer == 0){
		/* Orbital frequencies for circular equitorial orbits */
		CKerr_FindEQL_IRCirc(0, radius_outer, EQL_outer, M, astar);
		CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
		CKerr_Minverse(J_outer, Minv_outer, M, astar);
		//CKerr_Minv2Omega(Minv_outer, Omega_outer);
		//printf("Minv for circular equitorial case = \n %lg %lg %lg \n", Minv[0], Minv[1], Minv[2]);
		//printf(" %lg %lg %lg \n", Minv[3], Minv[4], Minv[5]);
		//printf(" %lg %lg %lg \n", Minv[6], Minv[7], Minv[8]);
		//printf("Js are: %lg %lg %lg \n", J[0], J[1], J[2]);
		//printf("Omega for circular equatorial is: %lg %lg %lg \n", Omega[0], Omega[1], Omega[2]);

		/* Torus Origin for circular equatorial orbits */
		//CKerr_TorusOrigin(J, xuorig, M, astar);
	}

	/* For GENERIC orbits */
	else{
		/* Routine to find location of resonant apocenter */
		//guess1 = rp + 3.;
		//guess2 = radius_outer - 3.;
		//apo = find_resonance_apo(n_res_inner, k_res_inner, m_res_inner, radius_outer, guess1, guess2, rp, I, astar, M);
		
		/* Get orbital frequencies for inner generic orbit */
		ra_rp_I2EQL(ra_outer, EQL_outer, rp_outer, I_outer, astar, M);
		CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
		CKerr_Minverse(J_outer, Minv_outer, M, astar);
		//CKerr_Minv2Omega(Minv_outer, Omega_outer);
		//printf("Minv for generic case = \n %lg %lg %lg \n", Minv[0], Minv[1], Minv[2]);
		//printf(" %lg %lg %lg \n", Minv[3], Minv[4], Minv[5]);
		//printf(" %lg %lg %lg \n", Minv[6], Minv[7], Minv[8]);
		//printf("Js are: %lg %lg %lg \n", J[0], J[1], J[2]);
		//printf("Omega for generic is: %lg %lg %lg \n", Omega[0], Omega[1], Omega[2]);

		/* Torus origin for inner body orbit */
		//CKerr_TorusOrigin(J, xuorig, M, astar);
		//for(i=0;i<6;i++) printf("xuorig[%d] = %lg\n", i, xuorig[i]);
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


  	E0_inner = (double*)malloc((size_t)(nl*10*sizeof(double)));
  	E1_inner = E0_inner + nl;
  	C0_inner = E1_inner + nl;

  	E0_outer = (double*)malloc((size_t)(nl*10*sizeof(double)));
  	E1_outer = E0_outer + nl;
  	C0_outer = E1_outer + nl;

	rH = M + sqrt(M*M - astar*astar*M*M);
	epsilon = sqrt(M*M - astar*astar*M*M) / (4 * M * rH);
	//Gamma = omega_dot(nl, n_res_inner, k_res_inner, m_res_inner, apo, rp, I, astar, M, radius_outer, 1e-4);
	//sgn_Gamma = Gamma/fabs(Gamma);
	//printf("Gamma and its sign are: %lg %lg \n", Gamma, sgn_Gamma);

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
	
	//printf("About to start for-loops \n");
	//printf("il \t i_n \t i_k \t i_m \t C01_in \t C02_in \t C03_in \t C04_in \t C01_out \t C02_out \t term \t another_term \t last_term \n");
	/* Resonant sum, so it is over integer mulitples of the initial resonance modes for the inner body */
	for (i = -N_res; i <= N_res; i++){
		i_n_inner = i*n_res_inner;
		i_k_inner = i*k_res_inner;
		i_m_inner = i*m_res_inner;
		i_n_outer = i*n_res_outer;
		i_k_outer = i*k_res_outer;
		i_m_outer = i*m_res_outer;
		Rtheta = cos(i * theta_res_F);
		Itheta = sin(i * theta_res_F);
		if((i_n_outer == 0 && i_k_outer == 0 && i_m_outer == 0) || (i_n_inner == 0 && i_k_inner == 0 && i_m_inner == 0))
			continue;
		/* Create the amplitude data for inner body*/
		/* Create the amplitude data for outer body*/
		CKerr_RadialFunc(Minv_outer, xuorig_outer, M, astar, i_n_outer, i_k_outer, i_m_outer, 20, 20, nl, C0_outer, &omegagw_outer, E0_outer);
		CKerr_RadialFunc(Minv_inner, xuorig_inner, M, astar, i_n_inner, i_k_inner, i_m_inner, 20, 20, nl, C0_inner, &omegagw_inner, E0_inner);
		printf("Frequencies from resonance and RadialFunc: %lg \t %lg \t %lg\n", omega_nkm, omegagw_inner, omegagw_outer);
		for (il = 0; (il) < nl; il++){
			omega_nkm = i_n_inner * Omega_inner[0] + i_k_inner * Omega_inner[1] + i_m_inner * Omega_inner[2];
			/* Get the scattering coefficients of the outer body radiation */
			CKerr_GWScatMatrix(M, astar, omega_nkm, i_m_inner, E0_outer[il], cscat, aux);
			/* Define frequency of inner body (omega_nkm) and other functions that depends on n,k,m,l */
			lambda = E0_inner[il] - 2 * astar * M * i_m_inner * omega_nkm + astar * M * astar * M * omega_nkm * omega_nkm - 2;
			P = omega_nkm - i_m_inner * astar * M / (2 * M * rH);
			numer = 256 * pow(2 * M * rH,5) * P * (P*P + 4 * epsilon*epsilon) * (P*P + 16 * epsilon*epsilon) * omega_nkm * omega_nkm * omega_nkm;
			C2 = ((lambda + 2)*(lambda + 2) + 4 * i_m_inner * astar * M * omega_nkm - 4 * astar*astar*M*M * omega_nkm*omega_nkm) * (lambda*lambda + 36 * i_m_inner * astar * M *omega_nkm - 36 * astar*astar*M*M * omega_nkm*omega_nkm) + (2 * lambda +3) * (96 * astar*astar*M*M * omega_nkm*omega_nkm - 48 * i_m_inner * astar * M * omega_nkm) + 144 * omega_nkm*omega_nkm * (M*M - astar*astar*M*M);
			alphankm = numer / C2;
			printf("%i \t %lg \t %lg \t %lg\n", il, omega_nkm, lambda, alphankm);

			printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n_inner, i_k_inner, i_m_inner, il,  C0_inner[4*il], C0_inner[4*il+1], C0_inner[4*il+2], C0_inner[4*il+3], C0_outer[4*il], C0_outer[4*il+1], C0_outer[4*il+2], C0_outer[4*il+3]);
			//printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n, i_k, i_m, il, C0[4*il], C0[4*il+1], C0[4*il+2], C0[4*il+3], Z_down_square, Z_out_square);
			
			/* For INNER/OUTER Body: */
			/* C0_{inner/outer}[4*il    ] = Re[Z_down, inner] */
			/* C0_{inner/outer}[4*il + 1] = Im[Z_down, inner] */
			/* C0_{inner/outer}[4*il + 2] = Re[Z_out, inner] */
			/* C0_{inner/outer}[4*il + 3] = Im[Z_out, inner] */


			/* Expanded out the terms in our expression for \dot{J}_{td} and kept only the real parts*/
			term = C0_inner[4*il+2] * cscat[0] * Rtheta * C0_outer[4*il] + C0_outer[4*il+1] * cscat[0] * C0_inner[4*il + 2] * Itheta + C0_outer[4*il] * cscat[1] * C0_inner[4*il+3] * Rtheta + C0_outer[4*il+1] * cscat[1] * C0_inner[4*il+3] * Itheta;
			another_term = -C0_outer[4*il] * C0_inner[4*il+3] * cscat[0] * Itheta + C0_outer[4*il+1] * C0_inner[4*il+3] * cscat[0] * Rtheta + C0_outer[4*il] * cscat[1] * C0_inner[4*il + 2] * Itheta - C0_outer[4*il+1] * cscat[1] * C0_inner[4*il+2] * Rtheta;
			last_term = C0_outer[4*il] * C0_inner[4*il] * Rtheta + C0_outer[4*il+1] * C0_inner[4*il] * Itheta - C0_outer[4*il] * C0_inner[4*il+1] * Itheta + C0_outer[4*il+1] * C0_inner[4*il+1] * Rtheta;
			printf("%lg \t %lg \t %lg \n", term, another_term, alphankm * last_term);

			/* J_dot of inner body due to tidal field of outer body */
			J_dot_td[0] += -i_n_inner * mu_outer * (term + another_term + alphankm * last_term) / (omega_nkm * omega_nkm * omega_nkm);
			J_dot_td[1] += -i_k_inner * mu_outer * (term + another_term + alphankm * last_term) / (omega_nkm * omega_nkm * omega_nkm);
			J_dot_td[2] += -i_m_inner * mu_outer * (term + another_term + alphankm * last_term) / (omega_nkm * omega_nkm * omega_nkm);
			//printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \n", il, i_n_inner, i_k_inner, i_m_inner, C0_inner[4*il], C0_inner[4*il+1], C0_inner[4*il+2], C0_inner[4*il+3], C0_outer[4*il], C0_outer[4*il+1], term, another_term, alphankm * last_term);
			//printf("%lg \t %lg \t %lg \n", *J_dot_r_tidal, *J_dot_theta_tidal, *J_dot_phi_tidal);
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
	free((char*)E0_inner);
	free((char*)E0_outer);
}

/* Keplerian expression for evolution of L_z, while only keeping the resonant prefactor terms. We will compare this to our general machine for computing J_phi_dot since it equals L_z_dot */
double J_dot_phi_Kepler(double mu_outer, double r_outer, double apo, double peri, double incline, double Theta_res){
//double resonant_angle;
double x_orbitplane_square, y_orbitplane_square;
double a, eccen;
double torque_res;

/* Semi-major axis and eccentricity as a function of apocenter and pericenter */
a = (apo + peri) / 2;
eccen = (apo - peri) / (apo + peri);
//resonant_angle = 2 * (Omega + omega);

/* x^2 and y^2 in the orbit plane after being averaged over mean anomly, M */
x_orbitplane_square = a*a * (1 + 4 * eccen*eccen) / 2;
y_orbitplane_square = a*a * (1 - eccen*eccen) / 2;

//printf("Semi-major axis: %lg \n", a);
//printf("Eccentricity: %lg \n", eccen);
//printf("<x^2>_orbitplane = %lg \n", x_orbitplane_square);
//printf("<y^2>_orbitplane = %lg \n", y_orbitplane_square);

/* Expression for L_z_dot in the observer frame and only keeping the prefactor of the term sin(-resonant_angle) */
torque_res = -3 * mu_outer / (8 * r_outer*r_outer*r_outer) * (x_orbitplane_square - y_orbitplane_square) * (1 + cos(incline))*(1 + cos(incline)) * sin(-2 * Theta_res);
return(torque_res);
}

