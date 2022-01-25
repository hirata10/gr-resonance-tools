#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CKerr.h"

/* dOmega_dt = dOmega/dJ_i * J_dot_i */
/* dOmega/dJ_i will be numerically computed using symmetric difference quotient */
/* N = (n, k, m) */

double dOmega_dt(int n_res, int k_res, int m_res, double ra, double rp, double I, double astar, double M, double radius_outer, double delta_t){
	int k;
	double Minv_plus[9], Minv_minus[9];
	double J_step[3], Jplus[3], Jminus[3];
	double EQL[3], J[3], Omega_inner_plus[3], Omega_inner_minus[3], omega_dot[3];
	double dOmegai_dJj[9];
	double Gamma;
	int nl = 3, nmax = 3, kmax = 3, mmax = 3;
	double Jdotr, Jdottheta, Jdotphi;

	ra_rp_I2EQL(ra, EQL, rp, I, astar, M);
	CKerr_EQL2J(EQL, J, M, astar, NULL);

	
	J_dot(nl, nmax, kmax, mmax, ra, rp, radius_outer, I, M, astar, &Jdotr, &Jdottheta, &Jdotphi);

	/* Step size in Jr, Jtheta, and Jphi */
	J_step[0] = Jdotr * delta_t;
	J_step[1] = Jdottheta * delta_t;
	J_step[2] = Jdotphi * delta_t;

	/* Defining a small increment in J in both +/- */
	for(k=0;k<3;k++){Jplus[k] = J[k] + J_step[k]; Jminus[k] = J[k] - J_step[k];}

	CKerr_Minverse(Jplus, Minv_plus, M, astar);
	CKerr_Minverse(Jminus, Minv_minus, M, astar);
	CKerr_Minv2Omega(Minv_plus, Omega_inner_plus);
	CKerr_Minv2Omega(Minv_minus, Omega_inner_minus);

	/* Symmetric derivative: d\Omega_i/dJ_j \approx (\Omega_i(J_j+J_step_j) - \Omega_i(J_j-J_step_j))/(2 * J_step_j) */
	/* dOmegai/dJj[0] = dOmegar/dJr, dOmegai/dJj[1] = dOmegar/dJtheta, ... */
	/* CHANGE WITH NEW J_step_j */
	omega_dot[0] = (Omega_inner_plus[0] - Omega_inner_minus[0])/(2 * delta_t);
	omega_dot[1] = (Omega_inner_plus[0] - Omega_inner_minus[0])/(2 * delta_t);
	omega_dot[2] = (Omega_inner_plus[0] - Omega_inner_minus[0])/(2 * delta_t);

	/* omega_dot[0] = omega_r_dot, omega_dot[1] = omega_theta_dot, omega_dot[2] = omega_phi_dot */
	/*omega_dot[0] = dOmegai_dJj[0] * Jdotr + dOmegai_dJj[1] * Jdottheta + dOmegai_dJj[2] * Jdotphi;
	omega_dot[1] = dOmegai_dJj[3] * Jdotr + dOmegai_dJj[4] * Jdottheta + dOmegai_dJj[5] * Jdotphi;
	omega_dot[2] = dOmegai_dJj[6] * Jdotr + dOmegai_dJj[7] * Jdottheta + dOmegai_dJj[8] * Jdotphi;*/
	/* Gamma = N \dot \Omega_dot  */
	Gamma = n_res * omega_dot[0] + k_res * omega_dot[1] + m_res * omega_dot[2];
	return(Gamma);
}

int main(){
	double jstep[3];
	int n, k, m;
	double apo_res, peri, incline, mass, spin, radius_outer, guess1, guess2;
	double Gamma;

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

	printf("Resonance condition for n,k,m : ");
	scanf("%i %i %i", &n, &k, &m);

	
	Gamma = dOmega_dt(n, k, m, apo_res, peri, incline, spin, mass, radius_outer, 1e-4);
	printf("Gamma is  = %lg \n", Gamma);
	//for (j=0;j<nl*nmax*kmax*mmax;j++){printf("%2d \t\t %19.12lE \n", j, J_dot_r[j]);}
	return(0);
}