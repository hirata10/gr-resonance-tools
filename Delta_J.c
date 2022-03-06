#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CKerr.h"

/* Computing the Delta J_tidal,i */
/* Outputs arrays with two entries: [0] --> real part; [1] --> imaginary part */
double Delta_J_tidal(int nl, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double apo, double rp, double radius_outer, double I, double M, double astar, double theta_res_F, double *Delta_J_r_tidal, double *Delta_J_theta_tidal, double *Delta_J_phi_tidal){
	int i, i_n_inner, i_k_inner, i_m_inner, il;
	int i_n_outer, i_k_outer, i_m_outer;
	double EQL_inner[3], J_inner[3], EQL_outer[3], J_outer[3], Minv_inner[9], Minv_outer[9], Omega_inner[3], Omega_outer, info[6], info_outer[6], xuorig_inner[6], xuorig_outer[6], cscat[16], aux[4];
	double term, another_term, last_term;
	double Gamma, sgn_Gamma;
	double omega_nkm, omegagw_inner, omegagw_outer;
	double rH;
	double Rtheta = 1, Itheta = 0;
	double epsilon, lambda, numer, C2, P, alphankm;
	double *C0_inner, *E0_inner, *E1_inner, *C0_outer, *E0_outer, *E1_outer;


	/* Converting apocenter, pericenter, and inclination angle into frequencies for: */
	/* Inner body */
	ra_rp_I2EQL(apo, EQL_inner, rp, I, astar, M);
	CKerr_EQL2J(EQL_inner, J_inner, M, astar, NULL);
	CKerr_Minverse(J_inner, Minv_inner, M, astar);
	CKerr_Minv2Omega(Minv_inner, Omega_inner);
	printf("Minv for inner body = \n %lg %lg %lg \n", Minv_inner[0], Minv_inner[1], Minv_inner[2]);
	printf(" %lg %lg %lg \n", Minv_inner[3], Minv_inner[4], Minv_inner[5]);
	printf(" %lg %lg %lg \n", Minv_inner[6], Minv_inner[7], Minv_inner[8]);


	/* Outer Body */
	CKerr_FindEQL_IRCirc(I, radius_outer, EQL_outer, M, astar);
	CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
	CKerr_Minverse(J_outer, Minv_outer, M, astar);
	printf("Minv for outer body = \n %lg %lg %lg \n", Minv_outer[0], Minv_outer[1], Minv_outer[2]);
	printf(" %lg %lg %lg \n", Minv_outer[3], Minv_outer[4], Minv_outer[5]);
	printf(" %lg %lg %lg \n", Minv_outer[6], Minv_outer[7], Minv_outer[8]);
	Omega_outer = Omega_outer_direct(radius_outer, M, astar);

	/* Torus origin for inner body orbit */
	CKerr_TorusOrigin(J_inner, xuorig_inner, M, astar);
	//CKerr_TorusOrigin(J_outer, xuorig_outer, M, astar);

	/* Torus origin for outer body orbit */
	CKerr_getData_CircEq(M, astar, radius_outer, info_outer);
	xuorig_outer[0] = radius_outer; xuorig_outer[1] = M_PI/2.; xuorig_outer[2] = 0.;
  	xuorig_outer[3] = 0.; xuorig_outer[4] = 0.;      xuorig_outer[5] = info_outer[0];


  	E0_inner = (double*)malloc((size_t)(nl*10*sizeof(double)));
  	E1_inner = E0_inner + nl;
  	C0_inner = E1_inner + nl;

  	E0_outer = (double*)malloc((size_t)(nl*10*sizeof(double)));
  	E1_outer = E0_outer + nl;
  	C0_outer = E1_outer + nl;

	rH = M + sqrt(M*M - astar*astar);
	epsilon = sqrt(M*M - astar*astar) / (4 * M * rH);
	Gamma = -1 * dOmega_dt(nl, n_res_inner, k_res_inner, m_res_inner, apo, rp, I, astar, M, radius_outer, 1e-4);
	sgn_Gamma = Gamma/fabs(Gamma);
	printf("Gamma and its sign are: %lg %lg \n", Gamma, sgn_Gamma);

	Delta_J_r_tidal[0] = 0;
	Delta_J_r_tidal[1] = 0; 
	Delta_J_theta_tidal[0] = 0;
	Delta_J_theta_tidal[1] = 0;	
	Delta_J_phi_tidal[0] = 0;
	Delta_J_phi_tidal[1] = 0;
	
	printf("About to start for-loops \n");

	/* Resonant sum, so it is over integer mulitples of the initial resonance modes for the inner body */
	for (i = -1; i <= 1; i++){
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
		CKerr_RadialFunc(Minv_outer, xuorig_outer, M, astar, i_n_outer, i_k_outer, i_m_outer, 14, 14, nl, C0_outer, &omegagw_outer, E0_outer);
		CKerr_RadialFunc(Minv_inner, xuorig_inner, M, astar, i_n_inner, i_k_inner, i_m_inner, 14, 14, nl, C0_inner, &omegagw_inner, E0_inner);
		for (il = 0; (il) < nl; il++){
			omega_nkm = i_n_inner * Omega_inner[0] + i_k_inner * Omega_inner[1] + i_m_inner * Omega_inner[2];
			/* Get the scattering coefficients of the outer body radiation */
			CKerr_GWScatMatrix(M, astar, omega_nkm, i_m_inner, E0_outer[il], cscat, aux);
			/* Define frequency of inner body (omega_nkm) and other functions that depends on n,k,m,l */
			lambda = -E0_inner[il] - 2 * astar * i_m_inner * omega_nkm + astar * astar * omega_nkm * omega_nkm - 2;
			P = omega_nkm - i_m_inner * astar / (2 * M * rH);
			numer = 256 * pow(2 * M * rH,5) * P * (P*P + 4 * epsilon*epsilon) * (P*P + 16 * epsilon*epsilon) * omega_nkm * omega_nkm * omega_nkm;
			C2 = ((lambda + 2)*(lambda + 2) + 4 * astar *omega_nkm - 4 * astar*astar * omega_nkm*omega_nkm) * (lambda*lambda + 36 * i_m_inner * astar *omega_nkm - 36 * astar*astar * omega_nkm*omega_nkm) + (2 * lambda +3) * (96 * astar*astar * omega_nkm*omega_nkm - 48 * i_m_inner * astar * omega_nkm) + 144 * omega_nkm*omega_nkm * (M*M - astar*astar);
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
			Delta_J_r_tidal[0] += -i_n_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm) * cos(i * theta_res_F + sgn_Gamma * M_PI/4.) * sqrt(2. * M_PI / fabs(Gamma));
			Delta_J_r_tidal[1] += -i_n_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm) * sin(i * theta_res_F + sgn_Gamma * M_PI/4.) * sqrt(2. * M_PI / fabs(Gamma));
			Delta_J_theta_tidal[0] += -i_k_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm) * cos(i * theta_res_F + sgn_Gamma * M_PI/4.) * sqrt(2. * M_PI / fabs(Gamma));
			Delta_J_theta_tidal[1] += -i_k_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm) * sin(i * theta_res_F + sgn_Gamma * M_PI/4.) * sqrt(2. * M_PI / fabs(Gamma));
			Delta_J_phi_tidal[0] += -i_m_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm) * cos(i * theta_res_F + sgn_Gamma * M_PI/4.) * sqrt(2. * M_PI / fabs(Gamma));
			Delta_J_phi_tidal[1] += -i_m_inner * (term + another_term + alphankm * last_term) / (2 * omega_nkm * omega_nkm * omega_nkm) * sin(i * theta_res_F + sgn_Gamma * M_PI/4.) * sqrt(2. * M_PI / fabs(Gamma));
				}
			}
	free((char*)E0_inner);
	free((char*)E0_outer);
}

int main(){
	int j;
	double J_dot_r, J_dot_theta, J_dot_phi;
	double J_dot_r_tidal, J_dot_theta_tidal, J_dot_phi_tidal;
	double Delta_J_r_tidal[2], Delta_J_theta_tidal[2], Delta_J_phi_tidal[2];
	int nl, nmax, kmax, mmax;
	int n_res_inner, k_res_inner, m_res_inner;
	int n_res_outer, k_res_outer, m_res_outer;
	double apo_res, peri, incline, mass, spin, radius_outer, guess1, guess2;

	printf("Enter pericenter: ");
	scanf("%lf", &peri);
	printf("Enter inlincation angle (radians): ");
	scanf("%lf", &incline);
	printf("Enter central mass: ");
	scanf("%lf", &mass);
	printf("Enter spin parameter of BH: ");
	scanf("%lf", &spin);
	printf("Enter radius of outer orbit: ");
	scanf("%lf", &radius_outer);
	printf("Enter corresponding apocenter at resonance: ");
	scanf("%lf", &apo_res);

	//printf("Number of modes and max n,k,m: ");
	//scanf("%i %i %i %i", &nl, &nmax, &kmax, &mmax);

	printf("Number of modes and mode vector for inner orbit: ");
	scanf("%i %i %i %i", &nl, &n_res_inner, &k_res_inner, &m_res_inner);

	printf("Mode vector for outer orbit: ");
	scanf("%i %i %i", &n_res_outer, &k_res_outer, &m_res_outer);

	/* Now allocate an array with the apropriate size */
	//J_dot_r = (double*)malloc((size_t)(nl*nmax*kmax*mmax*sizeof(double)));
	//J_dot_theta = (double*)malloc((size_t)(nl*nmax*kmax*mmax*sizeof(double)));
	//J_dot_phi = (double*)malloc((size_t)(nl*nmax*kmax*mmax*sizeof(double)));
	
	//J_dot(nl, nmax, kmax, mmax, apo_res, peri, radius_outer, incline, mass, spin, &J_dot_r, &J_dot_theta, &J_dot_phi);
	Delta_J_tidal(nl, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, apo_res, peri, radius_outer, incline, mass, spin, 0, Delta_J_r_tidal, Delta_J_theta_tidal, Delta_J_phi_tidal);
	printf("Delta_J_r_tidal Re and Im = %lg %lg \n", Delta_J_r_tidal[0], Delta_J_r_tidal[1]);
	printf("Delta_J_theta_tidal Re and Im = %lg %lg \n", Delta_J_theta_tidal[0], Delta_J_theta_tidal[1]);
	printf("Delta_J_phi_tidal Re and Im = %lg %lg \n", Delta_J_phi_tidal[0], Delta_J_phi_tidal[1]);
	//for (j=0;j<nl*nmax*kmax*mmax;j++){printf("%2d \t\t %19.12lE \n", j, J_dot_r[j]);}
	return(0);
}
