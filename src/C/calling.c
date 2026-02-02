#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "CKerr.h"
#include "globalpars_c.h"

#define N_LINES_MAX 30000

#ifdef IS_GAMMA
int main(){
	double jstep[3];
	int  n_res_inner, k_res_inner, m_res_inner, n_res_outer, k_res_outer, nl;
	double ra_inner, rp_inner, I_inner, ra_outer, rp_outer, I_outer;
	double astar, M, radius_outer, delta_t = 1e-4;
	//double apo_res, peri, incline, mass, spin, radius_outer, guess1, guess2;
	double Gamma[2];

	printf("Enter inner body pericenter: ");
	scanf("%lf", &rp_inner);
	printf("Enter inner body apocenter: ");
	scanf("%lf", &ra_inner);
	printf("Enter inner body inlincation angle (radians): ");
	scanf("%lf", &I_inner);

	printf("Enter outer pericenter: ");
	scanf("%lf", &rp_outer);
	printf("Enter outer apocenter: ");
	scanf("%lf", &ra_outer);
	printf("Enter outer inlincation angle (radians): ");
	scanf("%lf", &I_outer);

	printf("Enter central mass: ");
	scanf("%lf", &M);
	printf("Enter spin parameter of BH: ");
	scanf("%lf", &astar);
	printf("Enter radius of outer orbit (if applicable): ");
	scanf("%lf", &radius_outer);
	

	printf("Number of modes and mode vectors for inner and outer body: ");
	scanf("%i %i %i %i %i %i", &nl, &n_res_inner, &k_res_inner, &m_res_inner, &n_res_outer, &k_res_outer);

	omega_dot(nl, n_res_inner, k_res_inner, m_res_inner, n_res_outer, k_res_outer, ra_inner, rp_inner, I_inner, ra_outer, rp_outer, I_outer, astar, M, radius_outer, delta_t, Gamma);
	//omega_dot(nl, n, k, m, apo_res, peri, incline, spin, mass, radius_outer, 1e-4, Gamma);
	printf("Gamma for the inner and outer bodies are: %lg %lg \n", Gamma[0], Gamma[1]);
	//for (j=0;j<nl*nmax*kmax*mmax;j++){printf("%2d \t\t %19.12lE \n", j, J_dot_r[j]);}
	return(0);
}
#endif

#ifdef IS_RESFIND
int main(){
	int j;
	int n, k, m;
	double step, next_step, Omega_condition;
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
	//printf("Enter corresponding apocenter at resonance: ");
	//scanf("%lf", &apo_res);

	//printf("Number of modes and max n,k,m: ");
	//scanf("%i %i %i %i", &nl, &nmax, &kmax, &mmax);

	printf("Enter mode vector: ");
	scanf("%i %i %i", &n, &k, &m);

	//printf("Mode vector for outer orbit: ");
	//scanf("%i %i %i", &n_res_outer, &k_res_outer, &m_res_outer);

	guess1 = peri + 1.5;
	guess2 = radius_outer - 2.;

	apo_res = find_resonance_apo_OuterCirc(n, k, m, radius_outer, guess1, guess2, peri, incline, spin, mass);
	printf("Apocenter for this resonance is = %lg \n", apo_res);
	#if 0
	for (j = 0; j < 1000; j++){
		step = (guess2 - guess1)/1000;
		next_step = guess1 + step * j;
		Omega_condition = ra_rp_I2Omega_OuterCirc(n, k, m, radius_outer, next_step, peri, incline, spin, mass);
		printf("%i %lg %lg \n", j, next_step, Omega_condition);
	}
	#endif
	
	//for (j=0;j<nl*nmax*kmax*mmax;j++){printf("%2d \t\t %19.12lE \n", j, J_dot_r[j]);}
	return(0);
}
#endif

#ifdef IS_J_DOT_SF
int main(int argc, char **argv){
	int j;
	double J[3], EQL[3], anc[3];
	double J_dot_sf[3];
	double angle_space, angle_step, L_dot[50], Kepler_torque[50];
	int nmax = GLOBALPAR_nmax, kmax = GLOBALPAR_kmax, mmax = GLOBALPAR_mmax, nl_self = GLOBALPAR_nl_self;
	double mass = GLOBALPAR_M, spin = GLOBALPAR_astar, radius_outer = 0.;

	/* Input Js on command line */
	sscanf(argv[1], "%lg", &J[0]);
	sscanf(argv[2], "%lg", &J[1]);
	sscanf(argv[3], "%lg", &J[2]);
	

	CKerr_J2EQL(J, EQL, mass, spin);
	CKerr_EQL2J(EQL, J, mass, spin, anc);

	double I = anc[0];
	double rp = anc[1];
	double ra = anc[2];


	double eccen = (ra - rp) / (ra + rp);

	printf("I, rp, ra, eccentricity: %lg %lg %lg %lg \n", I, rp, ra, eccen);

	printf("EQL: %lg %lg %lg \n", EQL[0], EQL[1], EQL[2]);

	printf("About to compute self-force J_dots \n");
	J_dot_selfforce(nl_self, nmax, kmax, mmax, ra, rp, radius_outer, I, mass, spin, J_dot_sf);
	printf("J_dot_r_sf = %lg \n", J_dot_sf[0]);
	printf("J_dot_theta_sf = %lg \n", J_dot_sf[1]);
	printf("J_dot_phi_sf = %lg \n", J_dot_sf[2]);

	return(0);
}
#endif

#ifdef IS_J_DOT_TIDAL
int main(int argc, char **argv){
	int j;
	double J_dot_r, J_dot_theta, J_dot_phi;
	double J_dot_r_tidal, J_dot_theta_tidal, J_dot_phi_tidal;
	double J_dot_sf[3], J_dot_td[3];
	double J_inner[3], J_outer[3], EQL_inner[3], EQL_outer[3], anc_inner[3], anc_outer[3];
	double angle_space, angle_step, L_dot[50], Kepler_torque[50];
	int nl = GLOBALPAR_nl_res, nmax = GLOBALPAR_nmax, kmax = GLOBALPAR_kmax, mmax = GLOBALPAR_mmax, nl_self = GLOBALPAR_nl_self, N_res = GLOBALPAR_N_res;
	int n_res_inner, k_res_inner, m_res_inner;
	int n_res_outer, k_res_outer, m_res_outer;
	double ra_inner, rp_inner, I_inner, ra_outer, rp_outer, I_outer, radius_outer=0., guess1, guess2, angle_torus;
	double mass = GLOBALPAR_M, spin = GLOBALPAR_astar, mu_outer = GLOBALPAR_mu_outer;

	sscanf(argv[1], "%lg", &J_inner[0]);
	sscanf(argv[2], "%lg", &J_inner[1]);
	sscanf(argv[3], "%lg", &J_inner[2]);
	sscanf(argv[4], "%lg", &J_outer[0]);
	sscanf(argv[5], "%lg", &J_outer[1]);
	sscanf(argv[6], "%lg", &J_outer[2]);
	sscanf(argv[7], "%d", &n_res_inner);
	sscanf(argv[8], "%d", &k_res_inner);
	sscanf(argv[9], "%d", &m_res_inner);
	sscanf(argv[10], "%d", &n_res_outer);
	sscanf(argv[11], "%d", &k_res_outer);
	sscanf(argv[12], "%d", &m_res_outer);
	sscanf(argv[13], "%lg", &angle_torus);

	#if 0
	// printf("Enter inner pericenter: ");
	// scanf("%lf", &rp_inner);
	// printf("Enter inner inlincation angle (radians): ");
	// scanf("%lf", &I_inner);
	// printf("Enter central mass: ");
	// scanf("%lf", &mass);
	// printf("Enter spin parameter of BH: ");
	// scanf("%lf", &spin);
	printf("Enter radius of outer orbit (if applicable): ");
	scanf("%lf", &radius_outer);
	// printf("Enter inner apocenter: ");
	// scanf("%lf", &ra_inner);

	// printf("Enter outer pericenter: ");
	// scanf("%lf", &rp_outer);
	// printf("Enter outer inlincation angle (radians): ");
	// scanf("%lf", &I_outer);
	// printf("Enter outer apocenter: ");
	// scanf("%lf", &ra_outer);

	printf("Enter J_inner components: ");
	scanf("%lg %lg %lg", &J_inner[0], &J_inner[1], &J_inner[2]);

	printf("Enter J_outer components: ");
	scanf("%lg %lg %lg", &J_outer[0], &J_outer[1], &J_outer[2]);

	#endif

	CKerr_J2EQL(J_inner, EQL_inner, mass, spin);
	CKerr_EQL2J(EQL_inner, J_inner, mass, spin, anc_inner);

	CKerr_J2EQL(J_outer, EQL_outer, mass, spin);
	CKerr_EQL2J(EQL_outer, J_outer, mass, spin, anc_outer);

	I_inner = anc_inner[0];
	rp_inner = anc_inner[1];
	ra_inner = anc_inner[2];

	I_outer = anc_outer[0];
	rp_outer = anc_outer[1];
	ra_outer = anc_outer[2];


	// printf("Mode vector for inner orbit: ");
	// scanf("%i %i %i", &n_res_inner, &k_res_inner, &m_res_inner);

	// printf("Number of resonance terms: ");
	// scanf("%i", &N_res);

	// printf("Mode vector for outer orbit: ");
	// scanf("%i %i %i", &n_res_outer, &k_res_outer, &m_res_outer);

	double e_inner = (ra_inner - rp_inner) / (ra_inner + rp_inner);
	double e_outer = (ra_outer - rp_outer) / (ra_outer + rp_outer);

	printf("I_inner, rp_inner, ra_inner, eccen_inner: %lg %lg %lg %lg \n", I_inner, rp_inner, ra_inner, e_inner);
	printf("I_outer, rp_outer, ra_outer, eccen_outer: %lg %lg %lg %lg \n", I_outer, rp_outer, ra_outer, e_outer);
	printf("Inner and Outer EQL: %lg %lg %lg, %lg %lg %lg \n", EQL_inner[0], EQL_inner[1], EQL_inner[2], EQL_outer[0], EQL_outer[1], EQL_outer[2]);
	
	#if 1
	printf("About to compute tidal J_dots \n");
	J_dot_tidal(nl, N_res, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, ra_inner, rp_inner, radius_outer, I_inner, ra_outer, rp_outer, I_outer, mass, spin, angle_torus, mu_outer, J_dot_td);
	printf("J_dot_r_tidal = %lg \n", J_dot_td[0]);
	printf("J_dot_theta_tidal = %lg \n", J_dot_td[1]);
	printf("J_dot_phi_tidal = %lg \n", J_dot_td[2]);
	// printf("Keplerian J_dot_phi at resonance = %lg \n", J_dot_phi_Kepler(1.0, radius_outer, apo_res, peri, incline));
	#endif
	
	// printf("About to start L_dot \n");
	// for (j = 0; j < 50; j++){
	// 	angle_space = 2 * M_PI / 50;
	// 	angle_step = j * angle_space;
	// 	J_dot_tidal(nl, N_res, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, apo_res, peri, radius_outer, incline, mass, spin, angle_step, J_dot_td);
	// 	L_dot[j] = J_dot_td[2];
	// 	Kepler_torque[j] = J_dot_phi_Kepler(1.0, radius_outer, apo_res, peri, incline, angle_step);
	// 	printf("%lg %lg %lg \n", angle_step, L_dot[j], Kepler_torque[j]);
	// }
	//for (j=0;j<nl*nmax*kmax*mmax;j++){printf("%2d \t\t %19.12lE \n", j, J_dot_r[j]);}
	return(0);
}
#endif

#ifdef IS_J2EQL
int main(int argc, char **argv){
	double J[3], EQL[3], anc[3];
	double mass = GLOBALPAR_M, spin = GLOBALPAR_astar;

	/* Look for Jr, J_theta, and J_phi on command line */
	sscanf(argv[1], "%lg", &J[0]);
	sscanf(argv[2], "%lg", &J[1]);
	sscanf(argv[3], "%lg", &J[2]);

	CKerr_J2EQL(J, EQL, mass, spin);
	CKerr_EQL2J(EQL, J, mass, spin, anc);

	double I = anc[0];
	double rp = anc[1];
	double ra = anc[2];


	double eccen = (ra - rp) / (ra + rp);

	printf("I, rp, ra, eccentricity: %lg %lg %lg %lg \n", I, rp, ra, eccen);

	printf("EQL: %lg %lg %lg \n", EQL[0], EQL[1], EQL[2]);

	return(0);
}
#endif

#ifdef IS_DELTA_J
int main(){
	int j;
	double Delta_J_r_tidal, Delta_J_theta_tidal, Delta_J_phi_tidal;
	double Delta_J_tidal_value[3];
	int nl = GLOBALPAR_nl_res, N_res = GLOBALPAR_N_res;
	int n_res_inner, k_res_inner, m_res_inner;
	int n_res_outer, k_res_outer, m_res_outer;
	double EQL_inner[3], anc_inner[3], J_inner[3];
	double EQL_outer[3], anc_outer[3], J_outer[3];
	double apo_res, peri, incline, mass, spin, mu_outer = GLOBALPAR_mu_outer, radius_outer, guess1, guess2;
	double ra_inner, I_inner, rp_inner, angle_torus, ang_accel, eccentricity_inner;
	double ra_outer, I_outer, rp_outer, eccentricity_outer;

	// printf("Enter radius of outer orbit (if applicable): ");
	// scanf("%lf", &radius_outer);
	// printf("Enter inner apocenter: ");
	// scanf("%lf", &ra_inner);

	// printf("Enter outer pericenter: ");
	// scanf("%lf", &rp_outer);
	// printf("Enter outer inlincation angle (radians): ");
	// scanf("%lf", &I_outer);
	// printf("Enter outer apocenter: ");
	// scanf("%lf", &ra_outer);

	printf("Enter Mass and Spin of SMBH: ");
	scanf("%lg %lg", &mass, &spin);
	printf("Enter J_inner components: ");
	scanf("%lg %lg %lg", &J_inner[0], &J_inner[1], &J_inner[2]);
	printf("Enter J_outer components: ");
	scanf("%lg %lg %lg", &J_outer[0], &J_outer[1], &J_outer[2]);

	CKerr_J2EQL(J_inner, EQL_inner, mass, spin);
	CKerr_EQL2J(EQL_inner, J_inner, mass, spin, anc_inner);

	CKerr_J2EQL(J_outer, EQL_outer, mass, spin);
	CKerr_EQL2J(EQL_outer, J_outer, mass, spin, anc_outer);

	I_inner = anc_inner[0];
	rp_inner = anc_inner[1];
	ra_inner = anc_inner[2];

	I_outer = anc_outer[0];
	rp_outer = anc_outer[1];
	ra_outer = anc_outer[2];


	// printf("Location on torus (in radians): ");
	// scanf("%lg", &angle_torus);

	// printf("Angular acceleration: ");
	// scanf("%lg", &ang_accel);


	// printf("Mode vector for inner orbit: ");
	// scanf("%i %i %i", &n_res_inner, &k_res_inner, &m_res_inner);

	// // printf("Number of resonance terms: ");
	// // scanf("%i", &N_res);

	// printf("Mode vector for outer orbit: ");
	// scanf("%i %i %i", &n_res_outer, &k_res_outer, &m_res_outer);

	eccentricity_inner = (ra_inner - rp_inner)/(ra_inner + rp_inner);
	eccentricity_outer = (ra_outer - rp_outer)/(ra_outer + rp_outer);

	printf("I_inner, rp_inner, ra_inner, eccentricity_inner: %lg %lg %lg %lg \n", I_inner, rp_inner, ra_inner, eccentricity_inner);
	printf("I_outer, rp_outer, ra_outer, eccentricity_outer: %lg %lg %lg %lg \n", I_outer, rp_outer, ra_outer, eccentricity_outer);
	printf("Inner and Outer EQL: %lg %lg %lg, %lg %lg %lg \n", EQL_inner[0], EQL_inner[1], EQL_inner[2], EQL_outer[0], EQL_outer[1], EQL_outer[2]);
	#if 0
	/* Now allocate an array with the apropriate size */
	//J_dot_r = (double*)malloc((size_t)(nl*nmax*kmax*mmax*sizeof(double)));
	//J_dot_theta = (double*)malloc((size_t)(nl*nmax*kmax*mmax*sizeof(double)));
	//J_dot_phi = (double*)malloc((size_t)(nl*nmax*kmax*mmax*sizeof(double)));
	
	//J_dot(nl, nmax, kmax, mmax, apo_res, peri, radius_outer, incline, mass, spin, &J_dot_r, &J_dot_theta, &J_dot_phi);
	// Delta_J_tidal(nl, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, apo_inner, peri_inner, radius_outer, I_inner, apo_outer, peri_outer, I_outer, mass, spin, theta_res_F, Delta_J_tidal_value);
	Delta_J_tidal2(nl, N_res, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, ra_inner, rp_inner, radius_outer, I_inner, ra_outer, rp_outer, I_outer, mass, spin, angle_torus, ang_accel, mu_outer, Delta_J_tidal_value);
	printf("Delta_J_r_tidal = %lg \n", Delta_J_tidal_value[0]);
	printf("Delta_J_theta_tidal = %lg \n", Delta_J_tidal_value[1]);
	printf("Delta_J_phi_tidal = %lg \n", Delta_J_tidal_value[2]);
	//for (j=0;j<nl*nmax*kmax*mmax;j++){printf("%2d \t\t %19.12lE \n", j, J_dot_r[j]);}
	#endif
	return(0);
}
#endif

#ifdef IS_DELTA_J_SINGLE
int main(int argc, char **argv){
	int j;
	long i;
	double Delta_J_r_tidal, Delta_J_theta_tidal, Delta_J_phi_tidal;
	double Delta_J_tidal_value[3];
	int nl = GLOBALPAR_nl_res; /* Number of modes */
	int n_res_inner, k_res_inner, m_res_inner;
	int n_res_outer, k_res_outer, m_res_outer;
	int system_label, res_label;
	int N_res = GLOBALPAR_N_res;
	double ra_inner, I_inner, rp_inner;
	double ra_outer, I_outer, rp_outer;
	double Omega_res_inner_r, Omega_res_inner_theta, Omega_res_inner_phi;
	double Omega_res_outer_r, Omega_res_outer_theta, Omega_res_outer_phi;
	double J_res_inner_r, J_res_inner_theta, J_res_inner_phi;
	double J_res_outer_r, J_res_outer_theta, J_res_outer_phi;
	double ratio_Delta_J_J_r, ratio_Delta_J_J_theta, ratio_Delta_J_J_phi;
	double theta_res_F;
	double radius_outer = 0;
	double M = GLOBALPAR_M, astar = GLOBALPAR_astar; //Black hole mass set to 1 for units, BH spin wrt to BH mass = 1
	double ang_accel, mu_outer; //Value of angular acceleration (including body mass)

	/* Inputs to be given in command line */
	//sscanf(argv[1], "%ld", &nl);
	sscanf(argv[1], "%d", &n_res_inner);
	sscanf(argv[2], "%d", &n_res_outer);
	sscanf(argv[3], "%d", &k_res_inner);
	sscanf(argv[4], "%d", &k_res_outer);
	sscanf(argv[5], "%d", &m_res_inner);
	sscanf(argv[6], "%d", &m_res_outer);
	sscanf(argv[7], "%lg", &ra_inner);
	sscanf(argv[8], "%lg", &rp_inner);
	//sscanf(argv[10], "%lg", &radius_outer);
	sscanf(argv[9], "%lg", &I_inner);
	sscanf(argv[10], "%lg", &ra_outer);
	sscanf(argv[11], "%lg", &rp_outer);
	sscanf(argv[12], "%lg", &I_outer);
	//sscanf(argv[15], "%lg", &M);
	//sscanf(argv[16], "%lg", &astar);
	//sscanf(argv[17], "%lg", &theta_res_F);
	sscanf(argv[13], "%lg", &ang_accel);
	sscanf(argv[14], "%lg", &theta_res_F);
	sscanf(argv[15], "%lg", &mu_outer);
	sscanf(argv[16], "%d", &system_label);
	sscanf(argv[17], "%d", &res_label);
	sscanf(argv[18], "%lg", &Omega_res_inner_r);
	sscanf(argv[19], "%lg", &Omega_res_inner_theta);
	sscanf(argv[20], "%lg", &Omega_res_inner_phi);
	sscanf(argv[21], "%lg", &Omega_res_outer_r);
	sscanf(argv[22], "%lg", &Omega_res_outer_theta);
	sscanf(argv[23], "%lg", &Omega_res_outer_phi);
	sscanf(argv[24], "%lg", &J_res_inner_r);
	sscanf(argv[25], "%lg", &J_res_inner_theta);
	sscanf(argv[26], "%lg", &J_res_inner_phi);
	sscanf(argv[27], "%lg", &J_res_outer_r);
	sscanf(argv[28], "%lg", &J_res_outer_theta);
	sscanf(argv[29], "%lg", &J_res_outer_phi);

	//theta_res_F = (double)rand()/(double)RAND_MAX * 2*M_PI; //Random number generator from 0 to 2*pi for resonant angle value

	//printf("system label \t n_inner \t k_inner \t n_outer \t k_outer \t m \t Gamma \t Delta_J_r \t Delta_J_theta \t Delta_J_phi \n------------\n");
	Delta_J_tidal2(nl, N_res, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, ra_inner, rp_inner, radius_outer, I_inner, ra_outer, rp_outer, I_outer, M, astar, theta_res_F, ang_accel, mu_outer, Delta_J_tidal_value);
	ratio_Delta_J_J_r = Delta_J_tidal_value[0]/J_res_inner_r;
	ratio_Delta_J_J_theta = Delta_J_tidal_value[1]/J_res_inner_theta;
	ratio_Delta_J_J_phi = Delta_J_tidal_value[2]/J_res_inner_phi;
	printf("%i \t %i \t %i \t %i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \n", res_label, system_label, n_res_inner, k_res_inner, n_res_outer, k_res_outer, m_res_outer, theta_res_F, ang_accel, Omega_res_inner_r, Omega_res_inner_theta, Omega_res_inner_phi, Omega_res_outer_r, Omega_res_outer_theta, Omega_res_outer_phi, J_res_inner_r, J_res_inner_theta, J_res_inner_phi, J_res_outer_r, J_res_outer_theta, J_res_outer_phi, Delta_J_tidal_value[0], Delta_J_tidal_value[1], Delta_J_tidal_value[2], ratio_Delta_J_J_r, ratio_Delta_J_J_theta, ratio_Delta_J_J_phi);
	/* for(i = 0; i < MAX_ROW; i++){

		//ang_accel = data[i].data_col5 * data[i].data_col7; 

		Delta_J_tidal2(nl, data[i].data_col9, data[i].data_col11, data[i].data_col10, data[i].data_col12, data[i].data_col8, data[i].data_col8, data[i].data_col21, data[i].data_col20, 0, data[i].data_col19, data[i].data_col24, data[i].data_col23, data[i].data_col22, mass, spin, theta_res_F, data[i].data_col7, data[i].data_col6, Delta_J_tidal_value);
		printf("%i \t %i \t %i \t %i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \n", i, data[i].data_col1, data[i].data_col9, data[i].data_col10, data[i].data_col11, data[i].data_col12, data[i].data_col8, data[i].data_col7, Delta_J_tidal_value[0], Delta_J_tidal_value[1], Delta_J_tidal_value[2]);

	} */
	// Delta_J_tidal2(nl, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, apo_inner, peri_inner, radius_outer, I_inner, apo_outer, peri_outer, I_outer, mass, spin, theta_res_F, ang_accel, Delta_J_tidal_value);
	// printf("Delta_J_r_tidal = %lg \n", Delta_J_tidal_value[0]);
	// printf("Delta_J_theta_tidal = %lg \n", Delta_J_tidal_value[1]);
	// printf("Delta_J_phi_tidal = %lg \n", Delta_J_tidal_value[2]);
	//free((char*)data);
	return(0);
}
#endif

#ifdef IS_DELTA_OMEGA
int main(){

	int n_inner, k_inner, m_inner;
	int n_outer, k_outer;
	double value;

	n_inner=1;
	k_inner=2;
	m_inner=-3;

	n_outer = 1;
	k_outer = 1;

	value = ra_rp_I2Omega_generic(n_inner, k_inner, m_inner, n_outer, k_outer, 6, 3, 0.1, 20, 0, 0, 0.9, 1.0);
	printf("%lg \n", value);
	return(0);
}
#endif

#if IS_RK4_J_DOT
int main(int argc, char **argv)
{
	//double h=100, t, t0 = 1.; //Steps and initial start time
	double *J_r_final, *J_theta_final, *J_phi_final, *t;
	double J_r_ini, J_theta_ini, J_phi_ini, t0;
	double mu_body, M = GLOBALPAR_M, astar = GLOBALPAR_astar; //Mass of body and BH parameters
	long n; //Number of time steps
	long label; //System label from data file
	char sys_type[100]; //Inner or outer body label
	

	/* Inputs to be given in command line */
	/* Initial action variables (J_i), initial starting time (t0), number of steps (n), mass ratio of orbiting body, system label, and system type */
	sscanf(argv[1], "%lg", &J_r_ini);
	sscanf(argv[2], "%lg", &J_theta_ini);
	sscanf(argv[3], "%lg", &J_phi_ini);
	sscanf(argv[4], "%lg", &t0);
	sscanf(argv[5], "%ld", &n);
	sscanf(argv[6], "%lg", &mu_body);
	sscanf(argv[7], "%ld", &label);
	sscanf(argv[8], "%s", &sys_type);

	printf("Total number of arguments is %ld \n", argc);
	printf("Number of time steps is n = %ld \n", n);
	printf("Initial time is t0 = %lg \n", t0);
	printf("Mass ratio of the inspiral body = %lg \n", mu_body);
	printf("Mass and spin of central BH are = %lg %lg \n", M, astar);
	printf("Inital Js are: %lg %lg %lg\n", J_r_ini, J_theta_ini, J_phi_ini);
	printf("System label and type are: %ld %s \n", label, sys_type);

	
	J_r_final = (double *)malloc(sizeof(double) * n);
	J_theta_final = (double *)malloc(sizeof(double) * n);
	J_phi_final = (double *)malloc(sizeof(double) * n);

	/* Define time step array, starts at inital time (input) */
    t = (double *)malloc(sizeof(double) * (n));
    t[0] = t0;
	

	/* Make a .txt file with the date in the name or values of J_i_ini */
	char *filename[50];
	//struct tm *timenow;

	//time_t now = time(NULL);
	//timenow = gmtime(&now);

	//strftime(filename, sizeof(filename), "outputs_data/J_evolv_testrun_%Y-%m-%d_%H:%M:%S.txt", timenow);
	//sprintf(filename, "outputs_data/J_evolve_%ld.txt", label);
	sprintf(filename, "outputs_data/J_evolve_%s_%ld.txt", sys_type, label);
	

	FILE *fptr = fopen(filename,"w");


	rk4_J2Jdot(t, n, J_r_ini, J_theta_ini, J_phi_ini, J_r_final, J_theta_final, J_phi_final, fptr, mu_body, M, astar);
	
	fclose(fptr);
		 
 
	//printf("t \t J_r \t J_theta \t J_phi \n------------\n");
	/* for (i = 0; i < n; i++) {
		t = t0 + h * i;
		printf("%i \t %g \t %e \t %e \t %e\n", i, t, J_r_final[i], J_theta_final[i], J_phi_final[i]);
	} */
 
 	free((char*)J_r_final);
 	free((char*)J_theta_final);
	free((char*)J_phi_final);
	free((char*)t);
	return(0);

}
#endif

#if IS_RK4_J_DOT_start
int main(int argc, char **argv)
{
	//double h=100, t, t0 = 1.; //Steps and initial start time
	double *J_r_final, *J_theta_final, *J_phi_final, *t;
	double J_r_ini, J_theta_ini, J_phi_ini, t0;
	double mu_body, M = GLOBALPAR_M, astar = GLOBALPAR_astar; //Mass of body and BH parameters
	int n, i_start; //Number of time steps
	long label; //System label from data file
	char sys_type[100]; //Inner or outer body label
	int RESTART = 1; //Start from the beginning
	

	/* Inputs to be given in command line */
	/* Initial action variables (J_i), initial starting time (t0), number of steps (n), mass ratio of orbiting body, system label, and system type */
	sscanf(argv[1], "%lg", &J_r_ini);
	sscanf(argv[2], "%lg", &J_theta_ini);
	sscanf(argv[3], "%lg", &J_phi_ini);
	sscanf(argv[4], "%lg", &t0);
	sscanf(argv[5], "%ld", &n);
	sscanf(argv[6], "%lg", &mu_body);
	sscanf(argv[7], "%ld", &label);
	sscanf(argv[8], "%s", &sys_type);
	sscanf(argv[9], "%i", &i_start);

	printf("Total number of arguments is %ld \n", argc);
	printf("Number of time steps is n = %ld \n", n);
	printf("Initial time is t0 = %lg \n", t0);
	printf("Mass ratio of the inspiral body = %lg \n", mu_body);
	printf("Mass and spin of central BH are = %lg %lg \n", M, astar);
	printf("Inital Js are: %lg %lg %lg\n", J_r_ini, J_theta_ini, J_phi_ini);
	printf("System label and type are: %ld %s \n", label, sys_type);

	
	J_r_final = (double *)malloc(sizeof(double) * n);
	J_theta_final = (double *)malloc(sizeof(double) * n);
	J_phi_final = (double *)malloc(sizeof(double) * n);

	/* Define time step array, starts at inital time (input) */
    t = (double *)malloc(sizeof(double) * (n));
    // t[0] = t0;


	rk4_J2Jdot_restart(t0, t, i_start, n, J_r_ini, J_theta_ini, J_phi_ini, J_r_final, J_theta_final, J_phi_final, RESTART, mu_body, M, astar);
	
		 
 
	//printf("t \t J_r \t J_theta \t J_phi \n------------\n");
	/* for (i = 0; i < n; i++) {
		t = t0 + h * i;
		printf("%i \t %g \t %e \t %e \t %e\n", i, t, J_r_final[i], J_theta_final[i], J_phi_final[i]);
	} */
 
 	free((char*)J_r_final);
 	free((char*)J_theta_final);
	free((char*)J_phi_final);
	free((char*)t);
	return(0);

}
#endif

#if IS_RK4_J_DOT_restart
int main(int argc, char **argv)
{
	//double h=100, t, t0 = 1.; //Steps and initial start time
	double *J_r_final, *J_theta_final, *J_phi_final, *t;
	double J_r_ini, J_theta_ini, J_phi_ini, t_start;
	double mu_body, M = GLOBALPAR_M, astar = GLOBALPAR_astar; //Mass of body and BH parameters
	int n, i_start; //Number of time steps
	long label; //System label from data file
	char sys_type[100]; //Inner or outer body label
	int RESTART = 0; // Restart from intermediate time step
	

	/* Inputs to be given in command line */
	/* Initial action variables (J_i), initial starting time (t0), number of steps (n), mass ratio of orbiting body, system label, and system type */
	sscanf(argv[1], "%lg", &J_r_ini);
	sscanf(argv[2], "%lg", &J_theta_ini);
	sscanf(argv[3], "%lg", &J_phi_ini);
	sscanf(argv[4], "%lg", &t_start);
	sscanf(argv[5], "%i", &n);
	sscanf(argv[6], "%lg", &mu_body);
	sscanf(argv[7], "%ld", &label);
	sscanf(argv[8], "%s", &sys_type);
	sscanf(argv[9], "%i", &i_start);

	// printf("Total number of arguments is %ld \n", argc);
	// printf("Number of time steps is n = %ld \n", n);
	// printf("Initial time is t0 = %lg \n", t0);
	// printf("Mass ratio of the inspiral body = %lg \n", mu_body);
	// printf("Mass and spin of central BH are = %lg %lg \n", M, astar);
	// printf("Inital Js are: %lg %lg %lg\n", J_r_ini, J_theta_ini, J_phi_ini);
	// printf("System label and type are: %ld %s \n", label, sys_type);

	
	J_r_final = (double *)malloc(sizeof(double) * n);
	J_theta_final = (double *)malloc(sizeof(double) * n);
	J_phi_final = (double *)malloc(sizeof(double) * n);

	/* Define time step array, starts at inital time (input) */
    t = (double *)malloc(sizeof(double) * (n));
    // t[0] = t0;

	// printf("\n");
	rk4_J2Jdot_restart(t_start, t, i_start, n, J_r_ini, J_theta_ini, J_phi_ini, J_r_final, J_theta_final, J_phi_final, RESTART, mu_body, M, astar);
		 
 
	//printf("t \t J_r \t J_theta \t J_phi \n------------\n");
	/* for (i = 0; i < n; i++) {
		t = t0 + h * i;
		printf("%i \t %g \t %e \t %e \t %e\n", i, t, J_r_final[i], J_theta_final[i], J_phi_final[i]);
	} */
 
 	free((char*)J_r_final);
 	free((char*)J_theta_final);
	free((char*)J_phi_final);
	free((char*)t);
	return(0);

}
#endif

#if IS_DELTA_EQL
int main(int argc, char **argv)
{
	int i, j, k;
	double Delta_J_r, Delta_J_theta, Delta_J_phi, Delta_J[3], Delta_EQL[3];
	double J_r_res, J_theta_res, J_phi_res, J_res[3];
	double EQL_res[3], anc_res[3], Minv_res[9], M_res[9], Ident[9];
	double J_postres[3], EQL_postres[3], anc_postres[3];
	double M, astar; //Mass of body and BH parameters
	

	/* Inputs to be given in command line */
	/* Initial changes in action variables crossing tidal resonance action variables (Delta_J_i), 
	action variables at resonance (J_i_res), 
	mass and spin of central BH */
	sscanf(argv[1], "%lg", &Delta_J_r);
	sscanf(argv[2], "%lg", &Delta_J_theta);
	sscanf(argv[3], "%lg", &Delta_J_phi);
	sscanf(argv[4], "%lg", &J_r_res);
	sscanf(argv[5], "%lg", &J_theta_res);
	sscanf(argv[6], "%lg", &J_phi_res);
	sscanf(argv[7], "%lg", &M);
	sscanf(argv[8], "%lg", &astar);

	J_res[0] = J_r_res;
	J_res[1] = J_theta_res;
	J_res[2] = J_phi_res;

	CKerr_J2EQL(J_res, EQL_res, M, astar);
	CKerr_EQL2J(EQL_res, J_res, M, astar, anc_res);

	printf("At resonance: \n");
	printf("J_r = %lg, J_theta = %lg, J_phi = %lg \n", J_res[0], J_res[1], J_res[2]);
	printf("E = %lg, Q = %lg, L = %lg \n", EQL_res[0], EQL_res[1], EQL_res[2]);
	printf("apo = %lg, peri = %lg, inc = %lg \n", anc_res[2], anc_res[1], anc_res[0]);

	CKerr_Minverse(J_res, Minv_res, M, astar);

	invertMatrix(Minv_res, M_res);

	// for (int i = 0; i < 3; i++)
	// {
	// 	for (int j = 0; j < 3; j++)
	// 	{
	// 		for (int k = 0; k < 3; k++)
	// 		{
	// 			Ident[i + 3 * j] += Minv_res[i + 3 * j] * M_res[j + 3 * k]
	// 		}
			
	// 	}
		
	// }

	for (int i = 0; i < 3; i++) {
    	for (int j = 0; j < 3; j++) {
        	Ident[i * 3 + j] = 0.0;   // reset accumulator
        	for (int k = 0; k < 3; k++) {
            	Ident[i * 3 + j] += Minv_res[i * 3 + k] * M_res[k * 3 + j];
        	}
    	}
	}

	
	for (int i = 0; i < 9; i++)
	{
		printf("Ident[%d] = %lg \n", i, Ident[i]);
	}
	

	Delta_EQL[0] = M_res[0] * Delta_J_r + M_res[3] * Delta_J_theta + M_res[6] * Delta_J_phi;
	Delta_EQL[1] = M_res[1] * Delta_J_r + M_res[4] * Delta_J_theta + M_res[7] * Delta_J_phi;
	Delta_EQL[2] = M_res[2] * Delta_J_r + M_res[5] * Delta_J_theta + M_res[8] * Delta_J_phi;

	printf("Delta_EQL = %lg %lg %lg \n", Delta_EQL[0], Delta_EQL[1], Delta_EQL[2]);

	printf("After resonance: \n");
	for (int i = 0; i < 3; i++)
	{
		EQL_postres[i] = EQL_res[i] + Delta_EQL[i];
	}
	
	CKerr_EQL2J(EQL_postres, J_postres, M, astar, anc_postres);

	printf("J_r = %lg, J_theta = %lg, J_phi = %lg \n", J_postres[0], J_postres[1], J_postres[2]);
	printf("E = %lg, Q = %lg, L = %lg \n", EQL_postres[0], EQL_postres[1], EQL_postres[2]);
	printf("apo = %lg, peri = %lg, inc = %lg \n", anc_postres[2], anc_postres[1], anc_postres[0]);

	return(0);

}
#endif