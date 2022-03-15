#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CKerr.h"

/* \dot{\Delta \omega} = \vec{N}_{outer} \cdot \vec{\dot{\Omega}}_{outer} - \vec{N}_{inner} \cdot \vec{\dot{\Omega}}_{inner} */
/* \dot{\Omega}_j = dOmega_j_dt = dOmega_j/dJ_i * J_dot_i */
/* dOmega/dJ_i will be numerically computed using symmetric difference quotient */
/* N = (n, k, m) */

double omega_dot(int nl, int n_res_inner, int k_res_inner, int m_res_inner, double ra, double rp, double I, double astar, double M, double radius_outer, double delta_t){
	int k;
	double Minv_inner_plus[9], Minv_inner_minus[9];
	double Minv_outer_plus[9], Minv_outer_minus[9];
	double J_step_inner[3], Jplus_inner[3], Jminus_inner[3];
	double J_step_outer[3], Jplus_outer[3], Jminus_outer[3];
	double EQL_inner[3], J_inner[3], Omega_inner_plus[3], Omega_inner_minus[3], omega_dot_inner[3];
	double EQL_outer[3], J_outer[3], Omega_outer_plus[3], Omega_outer_minus[3], omega_dot_outer[3];
	double dOmegai_dJj[9];
	double Gamma;
	int nmax = 1, kmax = 1, mmax = 1;
	double Jdotr_inner, Jdottheta_inner, Jdotphi_inner;
	double Jdotr_outer, Jdottheta_outer, Jdotphi_outer;

	/* Inner Body: Generic Orbit */
	ra_rp_I2EQL(ra, EQL_inner, rp, I, astar, M);
	CKerr_EQL2J(EQL_inner, J_inner, M, astar, NULL);
	printf("Inner J's: %lg %lg %lg \n", J_inner[0], J_inner[1], J_inner[2]);

	/* Outer Body: Circular Equitorial Orbit */
	CKerr_FindEQL_IRCirc(0, radius_outer, EQL_outer, M, astar);
	CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
	printf("Outer J's: %lg %lg %lg \n", J_outer[0], J_outer[1], J_outer[2]);
	

	//WRITE A ROUTINE FOR BOTH INNER AND OUTER BODY INSPIRALS
	J_dot(nl, nmax, kmax, mmax, ra, rp, radius_outer, I, M, astar, &Jdotr_inner, &Jdottheta_inner, &Jdotphi_inner);
	J_dot(nl, nmax, kmax, mmax, 0, 0, radius_outer, 0, M, astar, &Jdotr_outer, &Jdottheta_outer, &Jdotphi_outer);
	printf("J dot for inner: %lg %lg %lg \n", Jdotr_inner, Jdottheta_inner, Jdotphi_inner);
	printf("J dot for outer: %lg %lg %lg \n", Jdotr_outer, Jdottheta_outer, Jdotphi_outer);

	/* Step size in Jr, Jtheta, and Jphi */
	J_step_inner[0] = Jdotr_inner * delta_t;
	J_step_inner[1] = Jdottheta_inner * delta_t;
	J_step_inner[2] = Jdotphi_inner * delta_t;
	J_step_outer[0] = Jdotr_outer * delta_t;
	J_step_outer[1] = Jdottheta_outer * delta_t;
	J_step_outer[2] = Jdotphi_outer * delta_t;
	printf("Steps in J for outer: %lg %lg %lg \n", J_step_outer[0], J_step_outer[1], J_step_outer[2]);

	/* Defining a small increment in J in both +/- */
	for(k=0;k<3;k++){Jplus_inner[k] = J_inner[k] + J_step_inner[k]; Jminus_inner[k] = J_inner[k] - J_step_inner[k];}
	for(k=0;k<3;k++){Jplus_outer[k] = J_outer[k] + J_step_outer[k]; Jminus_outer[k] = J_outer[k] - J_step_outer[k];}
	printf("J plus outer: %lg %lg %lg \n", Jplus_outer[0], Jplus_outer[1], Jplus_outer[2]);
	printf("J minus outer: %lg %lg %lg \n", Jminus_outer[0], Jminus_outer[1], Jminus_outer[2]);
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

	printf("Minv_plus for outer body = \n %lg %lg %lg \n", Minv_outer_plus[0], Minv_outer_plus[1], Minv_outer_plus[2]);
	printf(" %lg %lg %lg \n", Minv_outer_plus[3], Minv_outer_plus[4], Minv_outer_plus[5]);
	printf(" %lg %lg %lg \n",  Minv_outer_plus[6], Minv_outer_plus[7], Minv_outer_plus[8]);

	printf("Minv_minus for outer body = \n %lg %lg %lg \n", Minv_outer_minus[0], Minv_outer_minus[1], Minv_outer_minus[2]);
	printf(" %lg %lg %lg \n", Minv_outer_minus[3], Minv_outer_minus[4], Minv_outer_minus[5]);
	printf(" %lg %lg %lg \n",  Minv_outer_minus[6], Minv_outer_minus[7], Minv_outer_minus[8]);

	printf("Omega_outer for plus: %lg %lg %lg \n", Omega_outer_plus[0], Omega_outer_plus[1], Omega_outer_plus[2]);
	printf("Omega_outer for minus: %lg %lg %lg \n", Omega_outer_minus[0], Omega_outer_minus[1], Omega_outer_minus[2]);

	/* Symmetric derivative: d\Omega_i/dJ_j \approx (\Omega_i(J_j+J_step_j) - \Omega_i(J_j-J_step_j))/(2 * J_step_j) */
	/* dOmegai/dJj[0] = dOmegar/dJr, dOmegai/dJj[1] = dOmegar/dJtheta, ... */
	/* CHANGE WITH NEW J_step_j */
	omega_dot_inner[0] = (Omega_inner_plus[0] - Omega_inner_minus[0])/(2 * delta_t);
	omega_dot_inner[1] = (Omega_inner_plus[1] - Omega_inner_minus[1])/(2 * delta_t);
	omega_dot_inner[2] = (Omega_inner_plus[2] - Omega_inner_minus[2])/(2 * delta_t);

	omega_dot_outer[0] = (Omega_outer_plus[0] - Omega_outer_minus[0])/(2 * delta_t);
	omega_dot_outer[1] = (Omega_outer_plus[1] - Omega_outer_minus[1])/(2 * delta_t);
	omega_dot_outer[2] = (Omega_outer_plus[2] - Omega_outer_minus[2])/(2 * delta_t);

	printf("%lg %lg %lg \n", omega_dot_inner[0], omega_dot_inner[1], omega_dot_inner[2]);
	printf("%lg %lg %lg \n", omega_dot_outer[0], omega_dot_outer[1], omega_dot_outer[2]);

	/* omega_dot[0] = omega_r_dot, omega_dot[1] = omega_theta_dot, omega_dot[2] = omega_phi_dot */
	/*omega_dot[0] = dOmegai_dJj[0] * Jdotr + dOmegai_dJj[1] * Jdottheta + dOmegai_dJj[2] * Jdotphi;
	omega_dot[1] = dOmegai_dJj[3] * Jdotr + dOmegai_dJj[4] * Jdottheta + dOmegai_dJj[5] * Jdotphi;
	omega_dot[2] = dOmegai_dJj[6] * Jdotr + dOmegai_dJj[7] * Jdottheta + dOmegai_dJj[8] * Jdotphi;*/

	/* Gamma = \vec{N}_{outer} \cdot \vec{\dot{\Omega}}_{outer} - \vec{N}_{inner} \cdot \vec{\dot{\Omega}}_{inner}  */
	/* For circular equatorial orbits, we set n,k for outer obrit modes to zero. And the m mode matches the inner orbit */
	Gamma = -1 * n_res_inner * omega_dot_inner[0] - k_res_inner * omega_dot_inner[1] - m_res_inner * (omega_dot_inner[2] - omega_dot_outer[2]);
	return(Gamma);
}


int main(){
	double jstep[3];
	int n, k, m, nl;
	double apo_res, peri, incline, mass, spin, radius_outer, guess1, guess2;
	double Gam;

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

	printf("Resonance condition for nl,n,k,m : ");
	scanf("%i %i %i %i", &nl, &n, &k, &m);

	
	Gam = omega_dot(nl, n, k, m, apo_res, peri, incline, spin, mass, radius_outer, 1e-4);
	printf("Gamma is  = %lg \n", Gam);
	//for (j=0;j<nl*nmax*kmax*mmax;j++){printf("%2d \t\t %19.12lE \n", j, J_dot_r[j]);}
	return(0);
}