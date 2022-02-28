#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CKerr.h"

/* Computing the J_dot due to self-force from arXiv:1905.00030v2 Eq. (12) */

double J_dot(int nl, int nmax, int kmax, int mmax, double apo, double rp, double radius_outer, double I, double M, double astar, double *J_dot_r, double *J_dot_theta, double *J_dot_phi){
	int i, i_n, i_k, i_m, il;
	double EQL[3], J[3], Minv[9], Omega_inner[3], info[6], info_outer[6], xuorig[6], xuorig_outer[6];
	double Z_out_square, Z_down_square;
	double omega_nkm, omegagw;
	double rH;
	double J_dot_r_old, J_dot_theta_old, J_dot_phi_old;
	double epsilon, lambda, numer, C2, P, alphankm;
	double *C0, *E0, *E1;

	//printf("IM HERE 1 \n");
	ra_rp_I2EQL(apo, EQL, rp, I, astar, M);
	CKerr_EQL2J(EQL, J, M, astar, NULL);
	CKerr_Minverse(J, Minv, M, astar);
	CKerr_Minv2Omega(Minv, Omega_inner);

	/* Torus origin for inner body orbit */
	CKerr_TorusOrigin(J, xuorig, M, astar);
	for(i=0;i<6;i++) printf("xuorig[%d] = %lg\n", i, xuorig[i]);

	//printf("IM HERE 2 \n");
	/* Torus origin for outer body orbit */
	CKerr_getData_CircEq(M, astar, radius_outer, info_outer);
	xuorig_outer[0] = radius_outer; xuorig_outer[1] = M_PI/2.; xuorig_outer[2] = 0.;
  	xuorig_outer[3] = 0.; xuorig_outer[4] = 0.;      xuorig_outer[5] = info[0];

  	/* This takes the M_inv, the origin in phase space, mass and spin of BH, {n,k,m}, 
	number of steps in r and theta, number of modes (nl), and the angular frequency of the perturber orbit
	Returns the Z-coeff. and energy eigenvalues E[nl]. This is also
	computed at the specific resonance condition (1, 2, -2). */
	E0 = (double*)malloc((size_t)(nl*10*sizeof(double)));
  	E1 = E0 + nl;
  	C0 = E1 + nl;
	
	//printf("IM HERE 3 \n");

	rH = M + sqrt(M*M - astar*astar);
	epsilon = sqrt(M*M - astar*astar) / (4 * M * rH);

	*J_dot_r = 0;
	*J_dot_theta = 0; 
	*J_dot_phi = 0;
	//printf("IM HERE 4 \n");
	/* Compute the J_dots for the inner orbit */
	for (i_n = -nmax; i_n <= nmax; i_n++){
		for (i_k = -kmax; i_k <= kmax; i_k++){
			for (i_m = -mmax; i_m <= mmax; i_m++){
				if(i_n == 0 && i_k == 0 && i_m == 0)
					continue;
				CKerr_RadialFunc(Minv, xuorig, M, astar, i_n, i_k, i_m, 14, 14, nl, C0, &omegagw, E0);
				for (il = 0; (il) < nl; il++){
					
					omega_nkm = i_n * Omega_inner[0] + i_k * Omega_inner[1] + i_m * Omega_inner[2];
					lambda = -E0[il] - 2 * astar * i_m * omega_nkm;
					P = omega_nkm - i_m * astar / (2 * M * rH);
					numer = 256 * pow(2 * M * rH,5) * P * (P*P + 4 * epsilon*epsilon) * (P*P + 16 * epsilon*epsilon) * omega_nkm;
					C2 = ((lambda + 2)*(lambda + 2) + 4 * astar *omega_nkm - 4 * astar*astar * omega_nkm*omega_nkm) * (lambda*lambda + 36 * i_m * astar *omega_nkm - 36 * astar*astar * omega_nkm*omega_nkm) + (2 * lambda +3) * (96 * astar*astar * omega_nkm*omega_nkm - 48 * i_m * astar * omega_nkm) + 144 * omega_nkm*omega_nkm * (M*M - astar*astar);
					alphankm = numer / C2;
					//printf("%i \t %i \t %i \t %i \t %lf \t %lf \t %lf\n", i_n, i_k, i_k, il, omega_nkm, lambda, alphankm);


					Z_out_square = C0[4*il]*C0[4*il] + C0[4*il+1]*C0[4*il+1];
					Z_down_square = C0[4*il+2]*C0[4*il+2] + C0[4*il+3]*C0[4*il+3];
					//printf("%i \t %i \t %i \t %i \t %lf \t %lf\n", i_n, i_k, i_k, il, Z_down_square, Z_out_square);
					printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n, i_k, i_m, il, C0[4*il], C0[4*il+1], C0[4*il+2], C0[4*il+3], Z_down_square, Z_out_square);
					
					*J_dot_r += -i_n * (Z_out_square + alphankm * Z_down_square) / (2. * omega_nkm*omega_nkm*omega_nkm);
					*J_dot_theta += -i_k * (Z_out_square + alphankm * Z_down_square) / (2. * omega_nkm*omega_nkm*omega_nkm);
					*J_dot_phi += -i_m * (Z_out_square + alphankm * Z_down_square) / (2. * omega_nkm*omega_nkm*omega_nkm);

				}
			}
		}
	}
	free((char*)E0);
}

/* Computing the J_dot due to tidal field */
double J_dot_tidal(int nl, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double apo, double rp, double radius_outer, double I, double M, double astar, double *J_dot_r_tidal, double *J_dot_theta_tidal, double *J_dot_phi_tidal){
	int i, i_n, i_k, i_m, il;
	double EQL_inner[3], J_inner[3], EQL_outer[3], J_outer[3], Minv_inner[9], Minv_outer[9], Omega_inner[3], Omega_outer, info[6], info_outer[6], xuorig_inner[6], xuorig_outer[6], cscat[16], aux[4];
	double Z_out_inner, Z_out_inner_cc, Z_down_inner, Z_down_inner_cc;
	double Z_out_outer, Z_out_outer_cc, Z_down_outer, Z_down_outer_cc, c13, c13_cc;
	double term1, term2;
	double omega_nkm, omegagw_inner, omegagw_outer;
	double rH;
	double J_dot_r_old, J_dot_theta_old, J_dot_phi_old;
	double epsilon, lambda, numer, C2, P, alphankm;
	double *C0_inner, *E0_inner, *E1_inner, *C0_outer, *E0_outer, *E1_outer;


	/* Converting apocenter, pericenter, and inclination angle into frequencies for: */
	/* Inner body */
	ra_rp_I2EQL(apo, EQL_inner, rp, I, astar, M);
	CKerr_EQL2J(EQL_inner, J_inner, M, astar, NULL);
	CKerr_Minverse(J_inner, Minv_inner, M, astar);
	CKerr_Minv2Omega(Minv_inner, Omega_inner);

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
  	xuorig_outer[3] = 0.; xuorig_outer[4] = 0.;      xuorig_outer[5] = info[0];


  	E0_inner = (double*)malloc((size_t)(nl*10*sizeof(double)));
  	E1_inner = E0_inner + nl;
  	C0_inner = E1_inner + nl;

  	E0_outer = (double*)malloc((size_t)(nl*10*sizeof(double)));
  	E1_outer = E0_outer + nl;
  	C0_outer = E1_outer + nl;

	rH = M + sqrt(M*M - astar*astar);
	epsilon = sqrt(M*M - astar*astar) / (4 * M * rH);

	*J_dot_r_tidal = 0;
	*J_dot_theta_tidal = 0; 
	*J_dot_phi_tidal = 0;  	

	/* Create the amplitude data for outer body, still need to compute Minv_outer, though I think I can write this analytically */
	CKerr_RadialFunc(Minv_outer, xuorig_outer, M, astar, n_res_outer, k_res_outer, m_res_outer, 14, 14, nl, C0_outer, &omegagw_outer, E0_outer);

	for (i_n = -n_res_inner;  i_n <= 2 * n_res_inner; i_n = i_n + n_res_inner){
		for (i_k = -k_res_inner; i_k <= 2 * k_res_inner; i_k = i_k + k_res_inner){
			for (i_m = -m_res_inner; i_m <= 2 * m_res_inner; i_m = i_m + m_res_inner){
				if(i_n == 0 && i_k == 0 && i_m == 0)
					continue;
				/* Create the amplitude data for inner body*/
				CKerr_RadialFunc(Minv_inner, xuorig_inner, M, astar, i_n, i_k, i_m, 14, 14, nl, C0_inner, &omegagw_inner, E0_inner);
				for (il = 0; (il) < nl; il++){
					/* Get the scattering coefficients of the outer body radiation */
					CKerr_GWScatMatrix(M, astar, Omega_outer, i_m, E0_outer[il], cscat, aux);
					/* Define frequency of inner body (omega_nkm) and other functions that depends on n,k,m,l */
					omega_nkm = i_n * Omega_inner[0] + i_k * Omega_inner[1] + i_m * Omega_inner[2];
					lambda = -E0_inner[il] - 2 * astar * i_m * omega_nkm;
					P = omega_nkm - i_m * astar / (2 * M * rH);
					numer = 256 * pow(2 * M * rH,5) * P * (P*P + 4 * epsilon*epsilon) * (P*P + 16 * epsilon*epsilon) * omega_nkm;
					C2 = ((lambda + 2)*(lambda + 2) + 4 * astar *omega_nkm - 4 * astar*astar * omega_nkm*omega_nkm) * (lambda*lambda + 36 * i_m * astar *omega_nkm - 36 * astar*astar * omega_nkm*omega_nkm) + (2 * lambda +3) * (96 * astar*astar * omega_nkm*omega_nkm - 48 * i_m * astar * omega_nkm) + 144 * omega_nkm*omega_nkm * (M*M - astar*astar);
					alphankm = numer / C2;
					//printf("%i \t %i \t %i \t %i \t %lf \t %lf \t %lf\n", i_n, i_k, i_k, il, omega_nkm, lambda, alphankm);

					/* Amplitudes for inner bodyʻs gravitational radiation */
					Z_out_inner = C0_inner[4*il] + C0_inner[4*il+1];
					Z_down_inner = C0_inner[4*il+2] + C0_inner[4*il+3];
					Z_out_inner_cc = C0_inner[4*il] - C0_inner[4*il+1];
					Z_down_inner_cc = C0_inner[4*il+2] - C0_inner[4*il+3];

					/* Amplitudes for outer bodyʻs gravitational radiation */
					Z_out_outer = C0_outer[4*il] + C0_outer[4*il+1];
					Z_down_outer = C0_outer[4*il+2] + C0_outer[4*il+3];
					Z_out_outer_cc = C0_outer[4*il] - C0_outer[4*il+1];
					Z_down_outer_cc = C0_outer[4*il+2] - C0_outer[4*il+3];
					printf("%i \t %i \t %i \t %i \t %lf \t %lf \t %lf \t %lf\n", i_n, i_k, i_k, il, Z_out_inner, Z_down_inner, Z_out_outer, Z_down_outer);
					//printf("%i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg\n", i_n, i_k, i_m, il, C0[4*il], C0[4*il+1], C0[4*il+2], C0[4*il+3], Z_down_square, Z_out_square);
					
					/* Defining the c.c of the c13 scattering coefficient for the outer body */
					c13 = cscat[0] + cscat[1];
					c13_cc = cscat[0] - cscat[1];

					term1 = Z_down_outer_cc * (c13_cc * Z_out_inner + alphankm * Z_down_inner);
					term2 = Z_down_outer * (c13 * Z_out_inner_cc + alphankm * Z_down_inner_cc);

					/* J_dot of inner body due to tidal field of outer body */
					*J_dot_r_tidal += -i_n * (term1 + term2) / (2 * omega_nkm * omega_nkm * omega_nkm);
					*J_dot_theta_tidal += -i_k * (term1 + term2) / (2 * omega_nkm * omega_nkm * omega_nkm);
					*J_dot_phi_tidal += -i_m * (term1 + term2) / (2 * omega_nkm * omega_nkm * omega_nkm);
				}
			}
		}
	}
	free((char*)E0_inner);
	free((char*)E0_outer);
}



int main(){
	int j;
	double J_dot_r, J_dot_theta, J_dot_phi;
	double J_dot_r_tidal, J_dot_theta_tidal, J_dot_phi_tidal;
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
	J_dot_tidal(nl, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, apo_res, peri, radius_outer, incline, mass, spin, &J_dot_r_tidal, &J_dot_theta_tidal, &J_dot_phi_tidal);
	printf("J_dot_r_tidal = %lg \n", J_dot_r_tidal);
	printf("J_dot_theta_tidal = %lg \n", J_dot_theta_tidal);
	printf("J_dot_phi_tidal = %lg \n", J_dot_phi_tidal);
	//for (j=0;j<nl*nmax*kmax*mmax;j++){printf("%2d \t\t %19.12lE \n", j, J_dot_r[j]);}
	return(0);
}
