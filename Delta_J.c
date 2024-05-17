#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CKerr.h"

/* Computing the Delta J_tidal,i */
/* Input: resonance mode vector for inner and outer body, number of modes to sum over (nl), orbital data (apo, peri, incline), mass and spin of BH, and resonance angle */
/* Output: Delta_J_tidal_{r, \theta, \phi}/(mass_inner * mass_outer) */
void Delta_J_tidal(int nl, int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double ra_inner, double rp_inner, double radius_outer, double I_inner, double ra_outer, double rp_outer, double I_outer, double M, double astar, double theta_res_F, double *Delta_J){
	int i, i_n_inner, i_k_inner, i_m_inner, il;
	int i_n_outer, i_k_outer, i_m_outer;
	double EQL_inner[3], J_inner[3], EQL_outer[3], J_outer[3], Minv_inner[9], Minv_outer[9], Omega_inner[3], Omega_outer[3], info[6], info_outer[6], xuorig_inner[6], xuorig_outer[6], cscat[16], aux[4];
	double term, another_term, last_term;
	double Gamma[2], sgn_sGamma, tot_Gamma, mu_inner = 1, mu_outer = 1;
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
	//printf("Minv for inner body = \n %lg %lg %lg \n", Minv_inner[0], Minv_inner[1], Minv_inner[2]);
	//printf(" %lg %lg %lg \n", Minv_inner[3], Minv_inner[4], Minv_inner[5]);
	//printf(" %lg %lg %lg \n", Minv_inner[6], Minv_inner[7], Minv_inner[8]);

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

	/* Torus origin for inner body orbit */
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
	omega_dot(nl, n_res_inner, k_res_inner, m_res_inner, n_res_outer, k_res_outer, ra_inner, rp_inner, I_inner, ra_outer, rp_outer, I_outer, astar, M, radius_outer, 1.0e-4, Gamma);
	tot_Gamma = -Gamma[0] * mu_inner + Gamma[1] * mu_outer;
	//tot_Gamma = Gamma[0] * mu_inner;
	printf("Angular acceleration near resonance for inner body: %lg \n", Gamma[0]);
	printf("Angular acceleration near resonance for outer body: %lg \n", Gamma[1]);

	Delta_J[0] = 0;
	Delta_J[1] = 0;
	Delta_J[2] = 0;
	
	//printf("About to start for-loops \n");

	/* Resonant sum, so it is over integer mulitples of the initial resonance modes for the inner body */
	for (i = -N_res; i <= N_res; i++){
		i_n_inner = i*n_res_inner;
		i_k_inner = i*k_res_inner;
		i_m_inner = i*m_res_inner;
		i_n_outer = i*n_res_outer;
		i_k_outer = i*k_res_outer;
		i_m_outer = i*m_res_outer;
		sgn_sGamma = (i * tot_Gamma) / fabs(i * tot_Gamma);
		Rtheta = cos(i * theta_res_F + sgn_sGamma * M_PI/4.);
		Itheta = sin(i * theta_res_F + sgn_sGamma * M_PI/4.);
		if((i_n_outer == 0 && i_k_outer == 0 && i_m_outer == 0) || (i_n_inner == 0 && i_k_inner == 0 && i_m_inner == 0))
			continue;
		/* Create the amplitude data for inner body*/
		/* Create the amplitude data for outer body*/
		CKerr_RadialFunc(Minv_outer, xuorig_outer, M, astar, i_n_outer, i_k_outer, i_m_outer, 14, 14, nl, C0_outer, &omegagw_outer, E0_outer);
		CKerr_RadialFunc(Minv_inner, xuorig_inner, M, astar, i_n_inner, i_k_inner, i_m_inner, 14, 14, nl, C0_inner, &omegagw_inner, E0_inner);
		for (il = 0; (il) < nl; il++){
			omega_nkm = i_n_inner * Omega_inner[0] + i_k_inner * Omega_inner[1] + i_m_inner * Omega_inner[2];
			/* Get the scattering coefficients of the outer body radiation */
			CKerr_GWScatMatrix(M, astar, omega_nkm, i_m_inner, E0_outer[il], cscat, aux);
			/* Define frequency of inner body (omega_nkm) and other functions that depends on n,k,m,l */
			lambda = E0_inner[il] - 2 * astar * M * i_m_inner * omega_nkm + astar * M * astar * M * omega_nkm * omega_nkm - 2;
			P = omega_nkm - i_m_inner * astar * M / (2 * M * rH);
			numer = 256 * pow(2 * M * rH,5) * P * (P*P + 4 * epsilon*epsilon) * (P*P + 16 * epsilon*epsilon) * omega_nkm * omega_nkm * omega_nkm;
			C2 = ((lambda + 2)*(lambda + 2) + 4 * i_m_inner * astar * M *omega_nkm - 4 * astar * M * astar * M * omega_nkm*omega_nkm) * (lambda*lambda + 36 * i_m_inner * astar * M * omega_nkm - 36 * astar * M * astar * M * omega_nkm*omega_nkm) + (2 * lambda +3) * (96 * astar*astar * M * M * omega_nkm*omega_nkm - 48 * i_m_inner * astar * M * omega_nkm) + 144 * omega_nkm*omega_nkm * (M*M - astar*astar*M*M);
			alphankm = numer / C2;
			//printf("%i \t %i \t %i \t %i \t %lf \t %lf \t %lf\n", i_n, i_k, i_k, il, omega_nkm, lambda, alphankm);

			//printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n_inner, i_k_inner, i_m_inner, il,  C0_inner[4*il], C0_inner[4*il+1], C0_inner[4*il+2], C0_inner[4*il+3], C0_outer[4*il], C0_outer[4*il+1], C0_outer[4*il+2], C0_outer[4*il+3]);
			//printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n, i_k, i_m, il, C0[4*il], C0[4*il+1], C0[4*il+2], C0[4*il+3], Z_down_square, Z_out_square);

			/* Expanded out the terms in our expression for \dot{J}_{td} and kept only the real parts*/
			term = C0_inner[4*il+2] * cscat[0] * Rtheta * C0_outer[4*il] + C0_outer[4*il+1] * cscat[0] * C0_inner[4*il + 2] * Itheta + C0_outer[4*il] * cscat[1] * C0_inner[4*il+3] * Rtheta + C0_outer[4*il+1] * cscat[1] * C0_inner[4*il+3] * Itheta;
			another_term = -C0_outer[4*il] * C0_inner[4*il+3] * cscat[0] * Itheta + C0_outer[4*il+1] * C0_inner[4*il+3] * cscat[0] * Rtheta + C0_outer[4*il] * cscat[1] * C0_inner[4*il + 2] * Itheta - C0_outer[4*il+1] * cscat[1] * C0_inner[4*il+2] * Rtheta;
			last_term = C0_outer[4*il] * C0_inner[4*il] * Rtheta + C0_outer[4*il+1] * C0_inner[4*il] * Itheta - C0_outer[4*il] * C0_inner[4*il+1] * Itheta + C0_outer[4*il+1] * C0_inner[4*il+1] * Rtheta;
			#if 0
			/* J_dot of inner body due to tidal field of outer body */
			*J_dot_r_tidal += -i_n_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm) ;
			*J_dot_theta_tidal += -i_k_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm);
			*J_dot_phi_tidal += -i_m_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm);
			#endif
			/* \Delta J_td of inner body due to tidal field of outer body, for the real [0] and imaginary [1] parts */
			Delta_J[0] += -i_n_inner * (term + another_term + alphankm * last_term) / (omega_nkm * omega_nkm * omega_nkm)  * sqrt(2. * M_PI / fabs(i * tot_Gamma));
			Delta_J[1] += -i_k_inner * (term + another_term + alphankm * last_term) / (omega_nkm * omega_nkm * omega_nkm) * sqrt(2. * M_PI / fabs(i * tot_Gamma));
			Delta_J[2] += -i_m_inner * (term + another_term + alphankm * last_term) / (omega_nkm * omega_nkm * omega_nkm) * sqrt(2. * M_PI / fabs(i * tot_Gamma));
			//Delta_J_phi_tidal[1] += -i_m_inner * (term + another_term + alphankm * last_term) / (omega_nkm * omega_nkm * omega_nkm) * sqrt(2. * M_PI / fabs(tot_Gamma));
				}
			}
	free((char*)E0_inner);
	free((char*)E0_outer);
}

/* Computes Delta_J/mu_inner (resonant correction to action per unit inner body mass) but with the angular acceleration (Gamma including inner and outer body masses) as an input to the function */
void Delta_J_tidal2(int nl, int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double ra_inner, double rp_inner, double radius_outer, double I_inner, double ra_outer, double rp_outer, double I_outer, double M, double astar, double theta_res_F, double ang_accel, double mu_outer, double *Delta_J){
	int i, i_n_inner, i_k_inner, i_m_inner, il;
	int i_n_outer, i_k_outer, i_m_outer;
	double EQL_inner[3], J_inner[3], EQL_outer[3], J_outer[3], Minv_inner[9], Minv_outer[9], Omega_inner[3], Omega_outer[3], info[6], info_outer[6], xuorig_inner[6], xuorig_outer[6], cscat[16], aux[4];
	double term, another_term, last_term;
	//double Gamma[2], sgn_sGamma, tot_Gamma, mu_inner = 1e-5, mu_outer = 1e-5;
	double sgn_sGamma;
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
	//printf("Minv for inner body = \n %lg %lg %lg \n", Minv_inner[0], Minv_inner[1], Minv_inner[2]);
	//printf(" %lg %lg %lg \n", Minv_inner[3], Minv_inner[4], Minv_inner[5]);
	//printf(" %lg %lg %lg \n", Minv_inner[6], Minv_inner[7], Minv_inner[8]);

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

	/* Torus origin for inner body orbit */
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

	Delta_J[0] = 0;
	Delta_J[1] = 0;
	Delta_J[2] = 0;
	
	//printf("About to start for-loops \n");

	/* Resonant sum, so it is over integer mulitples of the initial resonance modes for the inner body */
	for (i = -N_res; i <= N_res; i++){
		i_n_inner = i*n_res_inner;
		i_k_inner = i*k_res_inner;
		i_m_inner = i*m_res_inner;
		i_n_outer = i*n_res_outer;
		i_k_outer = i*k_res_outer;
		i_m_outer = i*m_res_outer;
		sgn_sGamma = (i * ang_accel) / fabs(i * ang_accel);
		Rtheta = cos(i * theta_res_F + sgn_sGamma * M_PI/4.);
		Itheta = sin(i * theta_res_F + sgn_sGamma * M_PI/4.);
		if((i_n_outer == 0 && i_k_outer == 0 && i_m_outer == 0) || (i_n_inner == 0 && i_k_inner == 0 && i_m_inner == 0))
			continue;
		/* Create the amplitude data for inner body*/
		/* Create the amplitude data for outer body*/
		CKerr_RadialFunc(Minv_outer, xuorig_outer, M, astar, i_n_outer, i_k_outer, i_m_outer, 14, 14, nl, C0_outer, &omegagw_outer, E0_outer);
		CKerr_RadialFunc(Minv_inner, xuorig_inner, M, astar, i_n_inner, i_k_inner, i_m_inner, 14, 14, nl, C0_inner, &omegagw_inner, E0_inner);
		for (il = 0; (il) < nl; il++){
			omega_nkm = i_n_inner * Omega_inner[0] + i_k_inner * Omega_inner[1] + i_m_inner * Omega_inner[2];
			/* Get the scattering coefficients of the outer body radiation */
			CKerr_GWScatMatrix(M, astar, omega_nkm, i_m_inner, E0_outer[il], cscat, aux);
			/* Define frequency of inner body (omega_nkm) and other functions that depends on n,k,m,l */
			lambda = E0_inner[il] - 2 * astar * M * i_m_inner * omega_nkm + astar * M * astar * M * omega_nkm * omega_nkm - 2;
			P = omega_nkm - i_m_inner * astar * M / (2 * M * rH);
			numer = 256 * pow(2 * M * rH,5) * P * (P*P + 4 * epsilon*epsilon) * (P*P + 16 * epsilon*epsilon) * omega_nkm * omega_nkm * omega_nkm;
			C2 = ((lambda + 2)*(lambda + 2) + 4 * i_m_inner * astar * M *omega_nkm - 4 * astar * M * astar * M * omega_nkm*omega_nkm) * (lambda*lambda + 36 * i_m_inner * astar * M * omega_nkm - 36 * astar * M * astar * M * omega_nkm*omega_nkm) + (2 * lambda +3) * (96 * astar*astar * M * M * omega_nkm*omega_nkm - 48 * i_m_inner * astar * M * omega_nkm) + 144 * omega_nkm*omega_nkm * (M*M - astar*astar*M*M);
			alphankm = numer / C2;
			//printf("%i \t %i \t %i \t %i \t %lf \t %lf \t %lf\n", i_n, i_k, i_k, il, omega_nkm, lambda, alphankm);

			//printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n_inner, i_k_inner, i_m_inner, il,  C0_inner[4*il], C0_inner[4*il+1], C0_inner[4*il+2], C0_inner[4*il+3], C0_outer[4*il], C0_outer[4*il+1], C0_outer[4*il+2], C0_outer[4*il+3]);
			//printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n, i_k, i_m, il, C0[4*il], C0[4*il+1], C0[4*il+2], C0[4*il+3], Z_down_square, Z_out_square);

			/* Expanded out the terms in our expression for \dot{J}_{td} and kept only the real parts*/
			term = C0_inner[4*il+2] * cscat[0] * Rtheta * C0_outer[4*il] + C0_outer[4*il+1] * cscat[0] * C0_inner[4*il + 2] * Itheta + C0_outer[4*il] * cscat[1] * C0_inner[4*il+3] * Rtheta + C0_outer[4*il+1] * cscat[1] * C0_inner[4*il+3] * Itheta;
			another_term = -C0_outer[4*il] * C0_inner[4*il+3] * cscat[0] * Itheta + C0_outer[4*il+1] * C0_inner[4*il+3] * cscat[0] * Rtheta + C0_outer[4*il] * cscat[1] * C0_inner[4*il + 2] * Itheta - C0_outer[4*il+1] * cscat[1] * C0_inner[4*il+2] * Rtheta;
			last_term = C0_outer[4*il] * C0_inner[4*il] * Rtheta + C0_outer[4*il+1] * C0_inner[4*il] * Itheta - C0_outer[4*il] * C0_inner[4*il+1] * Itheta + C0_outer[4*il+1] * C0_inner[4*il+1] * Rtheta;
			#if 0
			/* J_dot of inner body due to tidal field of outer body */
			*J_dot_r_tidal += -i_n_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm) ;
			*J_dot_theta_tidal += -i_k_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm);
			*J_dot_phi_tidal += -i_m_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm);
			#endif
			/* \Delta J_td of inner body due to tidal field of outer body */
			Delta_J[0] += -i_n_inner * mu_outer * (term + another_term + alphankm * last_term) / (omega_nkm * omega_nkm * omega_nkm)  * sqrt(2. * M_PI / fabs(i * ang_accel));
			Delta_J[1] += -i_k_inner * mu_outer * (term + another_term + alphankm * last_term) / (omega_nkm * omega_nkm * omega_nkm) * sqrt(2. * M_PI / fabs(i * ang_accel));
			Delta_J[2] += -i_m_inner * mu_outer * (term + another_term + alphankm * last_term) / (omega_nkm * omega_nkm * omega_nkm) * sqrt(2. * M_PI / fabs(i * ang_accel));
			//Delta_J_phi_tidal[1] += -i_m_inner * (term + another_term + alphankm * last_term) / (omega_nkm * omega_nkm * omega_nkm) * sqrt(2. * M_PI / fabs(tot_Gamma));
				}
			}
	free((char*)E0_inner);
	free((char*)E0_outer);
}