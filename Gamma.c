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