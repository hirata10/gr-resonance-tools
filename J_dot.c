#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CKerr.h"

/* Computing the J_dot from arXiv:1905.00030v2 Eq. (12) */

double J_dot(int nl, int nmax, int kmax, int mmax, double apo, double rp, double radius_outer, double I, double M, double astar, double *J_dot_r, double *J_dot_theta, double *J_dot_phi){
	int i, i_n, i_k, i_m, il;
	double EQL[3], J[3], Minv[9], Omega_inner[3], info[6], info_outer[6], xuorig[6], xuorig_outer[6];
	double Z_out_square, Z_down_square;
	double omega_nkm, omegagw;
	double rH;
	double J_dot_r_old, J_dot_theta_old, J_dot_phi_old;
	double epsilon, lambda, numer, C2, P, alphankm;
	double *C0, *E0, *E1;

	printf("IM HERE 1 \n");
	ra_rp_I2EQL(apo, EQL, rp, I, astar, M);
	CKerr_EQL2J(EQL, J, M, astar, NULL);
	CKerr_Minverse(J, Minv, M, astar);
	CKerr_Minv2Omega(Minv, Omega_inner);

	/* Torus origin for inner body orbit */
	CKerr_TorusOrigin(J, xuorig, M, astar);
	for(i=0;i<6;i++) printf("xuorig[%d] = %lg\n", i, xuorig[i]);

	printf("IM HERE 2 \n");
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
	
	printf("IM HERE 3 \n");

	rH = M + sqrt(M*M - astar*astar);
	epsilon = sqrt(M*M - astar*astar) / (4 * M * rH);

	*J_dot_r = 0;
	*J_dot_theta = 0; 
	*J_dot_phi = 0;
	printf("IM HERE 4 \n");
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

int main(){
	int j;
	double J_dot_r, J_dot_theta, J_dot_phi;
	int nl, nmax, kmax, mmax;
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

	printf("Number of modes and max n,k,m : ");
	scanf("%i %i %i %i", &nl, &nmax, &kmax, &mmax);

	/* Now allocate an array with the apropriate size */
	//J_dot_r = (double*)malloc((size_t)(nl*nmax*kmax*mmax*sizeof(double)));
	//J_dot_theta = (double*)malloc((size_t)(nl*nmax*kmax*mmax*sizeof(double)));
	//J_dot_phi = (double*)malloc((size_t)(nl*nmax*kmax*mmax*sizeof(double)));
	
	J_dot(nl, nmax, kmax, mmax, apo_res, peri, radius_outer, incline, mass, spin, &J_dot_r, &J_dot_theta, &J_dot_phi);
	printf("IM HERE 5 \n");
	printf("J_dot_r = %lg \n", J_dot_r);
	printf("J_dot_theta = %lg \n", J_dot_theta);
	printf("J_dot_phi = %lg \n", J_dot_phi);
	//for (j=0;j<nl*nmax*kmax*mmax;j++){printf("%2d \t\t %19.12lE \n", j, J_dot_r[j]);}
	return(0);
}
