#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CKerr.h"

/* \dot{\Delta \omega} = \vec{N}_{outer} \cdot \vec{\dot{\Omega}}_{outer} - \vec{N}_{inner} \cdot \vec{\dot{\Omega}}_{inner} */
/* \dot{\Omega}_j = dOmega_j_dt = dOmega_j/dJ_i * J_dot_i */
/* dOmega/dJ_i will be numerically computed using symmetric difference quotient */
/* N = (n, k, m) */
/* Outputs the inner and outer body's change in angular velocity as an array Gamma[0] = Gamma_inner, Gamma[1] = Gamma_outer */
/* m_res_inner = m_res_outer by selection rules */

int omega_dot(int nl, int n_res_inner, int k_res_inner, int m_res_inner, int n_res_outer, int k_res_outer, double ra_inner, double rp_inner, double I_inner, double ra_outer, double rp_outer, double I_outer, double astar, double M, double radius_outer, double delta_t, double *Gamma){
	int k;
	double Minv_inner_plus[9], Minv_inner_minus[9];
	double Minv_outer_plus[9], Minv_outer_minus[9];
	double J_step_inner[3], Jplus_inner[3], Jminus_inner[3];
	double J_step_outer[3], Jplus_outer[3], Jminus_outer[3];
	double EQL_inner[3], J_inner[3], Omega_inner_plus[3], Omega_inner_minus[3], omega_dot_inner[3];
	double EQL_outer[3], J_outer[3], Omega_outer_plus[3], Omega_outer_minus[3], omega_dot_outer[3];
	double dOmegai_dJj[9];
	int nmax = 1, kmax = 1, mmax = 1;
	//double Jdotr_inner, Jdottheta_inner, Jdotphi_inner;
	double J_dot_inner[3], J_dot_outer[3];
	//double Jdotr_outer, Jdottheta_outer, Jdotphi_outer;

	/* Inner Body: Generic Orbit */
	ra_rp_I2EQL(ra_inner, EQL_inner, rp_inner, I_inner, astar, M);
	CKerr_EQL2J(EQL_inner, J_inner, M, astar, NULL);
	printf("Inner J's: %lg %lg %lg \n", J_inner[0], J_inner[1], J_inner[2]);

	/* Outer body: Generic Orbit */
	ra_rp_I2EQL(ra_outer, EQL_outer, rp_outer, I_outer, astar, M);
	CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
	printf("Inner J's: %lg %lg %lg \n", J_outer[0], J_outer[1], J_outer[2]);
	
	#if 0
	/* Outer Body: Circular Equitorial Orbit */
	CKerr_FindEQL_IRCirc(0, radius_outer, EQL_outer, M, astar);
	CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
	printf("Outer J's: %lg %lg %lg \n", J_outer[0], J_outer[1], J_outer[2]);
	#endif

	//WRITE A ROUTINE FOR BOTH INNER AND OUTER BODY INSPIRALS
	J_dot_selfforce(nl, nmax, kmax, mmax, ra_inner, rp_inner, radius_outer, I_inner, M, astar, J_dot_inner);
	J_dot_selfforce(nl, nmax, kmax, mmax, ra_outer, rp_outer, radius_outer, I_outer, M, astar, J_dot_outer);
	printf("J dot for inner: %lg %lg %lg \n", J_dot_inner[0], J_dot_inner[1], J_dot_inner[2]);
	printf("J dot for outer: %lg %lg %lg \n", J_dot_outer[0], J_dot_outer[1], J_dot_outer[2]);

	/* Step size in Jr, Jtheta, and Jphi */
	/* We set J_step_outer for Jr and Jtheta to zero since they should be identically zero for the case of a circular equatorial orbit */
	J_step_inner[0] = J_dot_inner[0] * delta_t;
	J_step_inner[1] = J_dot_inner[1] * delta_t;
	J_step_inner[2] = J_dot_inner[2] * delta_t;
	J_step_outer[0] = J_dot_outer[0] * delta_t;
	J_step_outer[1] = J_dot_outer[1] * delta_t;
	J_step_outer[2] = J_dot_outer[2] * delta_t;

	//J_step_outer[0] = Jdotr_outer * delta_t;
	//J_step_outer[1] = Jdottheta_outer * delta_t;
	//J_step_outer[2] = Jdotphi_outer * delta_t;
	printf("Steps in J for outer: %lg %lg %lg \n", J_step_outer[0], J_step_outer[1], J_step_outer[2]);

	/* Defining a small increment in J in both +/- */
	for(k=0;k<3;k++){Jplus_inner[k] = J_inner[k] + J_step_inner[k]; Jminus_inner[k] = J_inner[k] - J_step_inner[k];}
	for(k=0;k<3;k++){Jplus_outer[k] = J_outer[k] + J_step_outer[k]; Jminus_outer[k] = J_outer[k] - J_step_outer[k];}
	//printf("J plus outer: %lg %lg %lg \n", Jplus_outer[0], Jplus_outer[1], Jplus_outer[2]);
	//printf("J minus outer: %lg %lg %lg \n", Jminus_outer[0], Jminus_outer[1], Jminus_outer[2]);
	//if(Jminus_outer[0] < 0){Jminus_outer[0] = fabs(Jminus_outer[0]);}
	//printf("The new J minus outer if Jr<0: %lg %lg %lg \n", Jminus_outer[0], Jminus_outer[1], Jminus_outer[2]);


	/* Get frequencies for incremented J's of inner body */
	CKerr_Minverse(Jplus_inner, Minv_inner_plus, M, astar);
	CKerr_Minverse(Jminus_inner, Minv_inner_minus, M, astar);
	CKerr_Minv2Omega(Minv_inner_plus, Omega_inner_plus);
	CKerr_Minv2Omega(Minv_inner_minus, Omega_inner_minus);
	
	/* Get frequencies for incremented J's of outer body */
	CKerr_Minverse(Jplus_outer, Minv_outer_plus, M, astar);
	CKerr_Minverse(Jminus_outer, Minv_outer_minus, M, astar);
	CKerr_Minv2Omega(Minv_outer_plus, Omega_outer_plus);
	CKerr_Minv2Omega(Minv_outer_minus, Omega_outer_minus);

	#if 0
	printf("Minv_plus for outer body = \n %lg %lg %lg \n", Minv_outer_plus[0], Minv_outer_plus[1], Minv_outer_plus[2]);
	printf(" %lg %lg %lg \n", Minv_outer_plus[3], Minv_outer_plus[4], Minv_outer_plus[5]);
	printf(" %lg %lg %lg \n",  Minv_outer_plus[6], Minv_outer_plus[7], Minv_outer_plus[8]);

	printf("Minv_minus for outer body = \n %lg %lg %lg \n", Minv_outer_minus[0], Minv_outer_minus[1], Minv_outer_minus[2]);
	printf(" %lg %lg %lg \n", Minv_outer_minus[3], Minv_outer_minus[4], Minv_outer_minus[5]);
	printf(" %lg %lg %lg \n",  Minv_outer_minus[6], Minv_outer_minus[7], Minv_outer_minus[8]);

	printf("Omega_outer for plus: %lg %lg %lg \n", Omega_outer_plus[0], Omega_outer_plus[1], Omega_outer_plus[2]);
	printf("Omega_outer for minus: %lg %lg %lg \n", Omega_outer_minus[0], Omega_outer_minus[1], Omega_outer_minus[2]);
	#endif

	/* Symmetric derivative: d\Omega_i/dJ_j \approx (\Omega_i(J_j+J_step_j) - \Omega_i(J_j-J_step_j))/(2 * J_step_j) */
	/* dOmegai/dJj[0] = dOmegar/dJr, dOmegai/dJj[1] = dOmegar/dJtheta, ... */
	/* CHANGE WITH NEW J_step_j */
	omega_dot_inner[0] = (Omega_inner_plus[0] - Omega_inner_minus[0])/(2 * delta_t);
	omega_dot_inner[1] = (Omega_inner_plus[1] - Omega_inner_minus[1])/(2 * delta_t);
	omega_dot_inner[2] = (Omega_inner_plus[2] - Omega_inner_minus[2])/(2 * delta_t);

	omega_dot_outer[0] = (Omega_outer_plus[0] - Omega_outer_minus[0])/(2 * delta_t);
	omega_dot_outer[1] = (Omega_outer_plus[1] - Omega_outer_minus[1])/(2 * delta_t);
	omega_dot_outer[2] = (Omega_outer_plus[2] - Omega_outer_minus[2])/(2 * delta_t);

	#if 0
	printf("%lg %lg %lg \n", omega_dot_inner[0], omega_dot_inner[1], omega_dot_inner[2]);
	printf("%lg %lg %lg \n", omega_dot_outer[0], omega_dot_outer[1], omega_dot_outer[2]);
	#endif

	/* omega_dot[0] = omega_r_dot, omega_dot[1] = omega_theta_dot, omega_dot[2] = omega_phi_dot */
	/*omega_dot[0] = dOmegai_dJj[0] * Jdotr + dOmegai_dJj[1] * Jdottheta + dOmegai_dJj[2] * Jdotphi;
	omega_dot[1] = dOmegai_dJj[3] * Jdotr + dOmegai_dJj[4] * Jdottheta + dOmegai_dJj[5] * Jdotphi;
	omega_dot[2] = dOmegai_dJj[6] * Jdotr + dOmegai_dJj[7] * Jdottheta + dOmegai_dJj[8] * Jdotphi;*/

	/* Gamma = \vec{N}_{outer} \cdot \vec{\dot{\Omega}}_{outer} - \vec{N}_{inner} \cdot \vec{\dot{\Omega}}_{inner}  */
	/* For circular equatorial orbits, we set n,k for outer obrit modes to zero. And the m mode matches the inner orbit */
	Gamma[0] = n_res_inner * omega_dot_inner[0] + k_res_inner * omega_dot_inner[1] + m_res_inner * omega_dot_inner[2];
	Gamma[1] = n_res_outer * omega_dot_outer[0] + k_res_outer * omega_dot_outer[1] + m_res_inner * omega_dot_outer[2];
	return(0);
}

/* Eq. (17) of arXiv:2207.07733 with the analytic dOmega/dJ from the contour integrals (J2dOmega_dJ in kerrphase.c) */
/* Gamma = \vec{N}_{outer} \cdot (dOmega_outer/dJ) \cdot \dot{J}_{outer} - \vec{N}_{inner} \cdot (dOmega_inner/dJ) \cdot \dot{J}_{inner} */
/* Input: Js and J_dots (per unit mass ratio, as returned by J_dot_selfforce) of both bodies, mass ratios, mass and spin of BH */
/* With mu_inner = mu_outer = 1, Gamma is in the same (per mass ratio) units as gamma in potential_resonances_*.txt */
/* Output: Gamma[0] = N_inner.Omega_dot_inner, Gamma[1] = N_outer.Omega_dot_outer, Gamma[2] = Gamma[1] - Gamma[0] */
/* m_res_inner = m_res_outer = m_res by selection rules */
int Gamma_dOmegadJ(int n_inner, int k_inner, int n_outer, int k_outer, int m_res, double *J_inner, double *J_outer, double *J_dot_inner, double *J_dot_outer, double mu_inner, double mu_outer, double M, double astar, double *Gamma){
	int i, k;
	double dOmega_dJ_inner[9], dOmega_dJ_outer[9];
	double Omega_dot_inner[3], Omega_dot_outer[3];

	J2dOmega_dJ(J_inner, M, astar, dOmega_dJ_inner);
	J2dOmega_dJ(J_outer, M, astar, dOmega_dJ_outer);

	/* \dot{\Omega}_i = dOmega_i/dJ_k * mu * \dot{J}_k */
	for (i = 0; i < 3; i++){
		Omega_dot_inner[i] = 0.;
		Omega_dot_outer[i] = 0.;
		for (k = 0; k < 3; k++){
			Omega_dot_inner[i] += dOmega_dJ_inner[3 * i + k] * mu_inner * J_dot_inner[k];
			Omega_dot_outer[i] += dOmega_dJ_outer[3 * i + k] * mu_outer * J_dot_outer[k];
		}
	}

	Gamma[0] = n_inner * Omega_dot_inner[0] + k_inner * Omega_dot_inner[1] + m_res * Omega_dot_inner[2];
	Gamma[1] = n_outer * Omega_dot_outer[0] + k_outer * Omega_dot_outer[1] + m_res * Omega_dot_outer[2];
	Gamma[2] = Gamma[1] - Gamma[0];
	return(0);
}

/* Same as Gamma_dOmegadJ, but computes the J_dots from the Js with the OpenMP self-force (J_dot_selfforce_openmp) */
/* Output: J_dot_inner[3], J_dot_outer[3] (so they can be reused), and Gamma[3] as in Gamma_dOmegadJ */
int Gamma_dOmegadJ_selfforce(int nl, int nmax, int kmax, int mmax, int Nr, int Nt, int n_inner, int k_inner, int n_outer, int k_outer, int m_res, double *J_inner, double *J_outer, double mu_inner, double mu_outer, double M, double astar, double *J_dot_inner, double *J_dot_outer, double *Gamma){
	double J_local[3], EQL[3], anc[3];
	int i;

	/* Inner body: J -> EQL -> (I, rp, ra) -> J_dot */
	for (i = 0; i < 3; i++){J_local[i] = J_inner[i];}
	CKerr_J2EQL(J_local, EQL, M, astar);
	CKerr_EQL2J(EQL, J_local, M, astar, anc);
	J_dot_selfforce_openmp(nl, nmax, kmax, mmax, Nr, Nt, anc[2], anc[1], anc[0], M, astar, J_dot_inner);

	/* Outer body */
	for (i = 0; i < 3; i++){J_local[i] = J_outer[i];}
	CKerr_J2EQL(J_local, EQL, M, astar);
	CKerr_EQL2J(EQL, J_local, M, astar, anc);
	J_dot_selfforce_openmp(nl, nmax, kmax, mmax, Nr, Nt, anc[2], anc[1], anc[0], M, astar, J_dot_outer);

	return(Gamma_dOmegadJ(n_inner, k_inner, n_outer, k_outer, m_res, J_inner, J_outer, J_dot_inner, J_dot_outer, mu_inner, mu_outer, M, astar, Gamma));
}