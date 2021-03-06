#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "CKerr.h"

#ifdef IS_GAMMA
int main(){
	double jstep[3];
	int n, k, m, nl;
	double apo_res, peri, incline, mass, spin, radius_outer, guess1, guess2;
	double Gamma[2];

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

	
	//omega_dot(nl, n, k, m, apo_res, peri, incline, spin, mass, radius_outer, 1e-4, Gamma);
	printf("Gamma for inner and outer body are: %lg %lg \n", Gamma[0], Gamma[1]);
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

	guess1 = peri + 2.;
	guess2 = radius_outer - 2.;

	apo_res = find_resonance_apo_OuterCirc(n, k, m, radius_outer, guess1, guess2, peri, incline, spin, mass);
	printf("Apocenter for this resonance is = %lg \n", apo_res);
	#if 0
	for (j = 0; j <= 500; j++){
		step = (guess2 - guess1)/500;
		next_step = guess1 + step * j;
		Omega_condition = ra_rp_I2Omega(n, k, m, radius_outer, next_step, peri, incline, spin, mass);
		printf("%i %lg %lg \n", j, next_step, Omega_condition);
	}
	
	printf("Delta_J_r_tidal Re & Im = %lg %lg \n", Delta_J_r_tidal[0], Delta_J_r_tidal[1]);
	printf("Delta_J_theta_tidal = %lg %lg \n", Delta_J_theta_tidal[0], Delta_J_theta_tidal[1]);
	printf("Delta_J_phi_tidal = %lg %lg \n", Delta_J_phi_tidal[0], Delta_J_phi_tidal[1]);
	#endif
	//for (j=0;j<nl*nmax*kmax*mmax;j++){printf("%2d \t\t %19.12lE \n", j, J_dot_r[j]);}
	return(0);
}
#endif

#ifdef IS_J_DOT
int main(){
	int j;
	double J_dot_r, J_dot_theta, J_dot_phi;
	double J_dot_r_tidal, J_dot_theta_tidal, J_dot_phi_tidal;
	double J_dot_sf[3], J_dot_td[3];
	double Delta_J_r_tidal[2], Delta_J_theta_tidal[2], Delta_J_phi_tidal[2];
	double angle_space, angle_step, L_dot[50], Kepler_torque[50];
	int nl, nmax, kmax, mmax, N_res;
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

	#if 0
	printf("Number of modes and max n,k,m: ");
	scanf("%i %i %i %i", &nl, &nmax, &kmax, &mmax);
	#endif

	printf("Number of modes and mode vector for inner orbit: ");
	scanf("%i %i %i %i", &nl, &n_res_inner, &k_res_inner, &m_res_inner);

	printf("Number of resonance terms: ");
	scanf("%i", &N_res);

	printf("Mode vector for outer orbit: ");
	scanf("%i %i %i", &n_res_outer, &k_res_outer, &m_res_outer);
	
	
	#if 0
	printf("About to compute tidal J_dots \n");
	J_dot_tidal(nl, N_res, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, apo_res, peri, radius_outer, incline, mass, spin, -M_PI/2., &J_dot_r_tidal, &J_dot_theta_tidal, &J_dot_phi_tidal);
	printf("J_dot_r_tidal = %lg \n", J_dot_r_tidal);
	printf("J_dot_theta_tidal = %lg \n", J_dot_theta_tidal);
	printf("J_dot_phi_tidal = %lg \n", J_dot_phi_tidal);
	printf("Keplerian J_dot_phi at resonance = %lg \n", J_dot_phi_Kepler(1.0, radius_outer, apo_res, peri, incline));
	#endif

	
	#if 0
	J_dot_selfforce(nl, nmax, kmax, mmax, apo_res, peri, radius_outer, incline, mass, spin, J_dot_sf);
	printf("J_dot_r_sf = %lg \n", J_dot_sf[0]);
	printf("J_dot_theta_sf = %lg \n", J_dot_sf[1]);
	printf("J_dot_phi_sf = %lg \n", J_dot_sf[2]);
	#endif

	
	printf("About to start L_dot \n");
	for (j = 0; j < 50; j++){
		angle_space = 2 * M_PI / 50;
		angle_step = j * angle_space;
		J_dot_tidal(nl, N_res, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, apo_res, peri, radius_outer, incline, mass, spin, angle_step, J_dot_td);
		L_dot[j] = J_dot_td[2];
		Kepler_torque[j] = J_dot_phi_Kepler(1.0, radius_outer, apo_res, peri, incline, angle_step);
		printf("%lg %lg %lg \n", angle_step, L_dot[j], Kepler_torque[j]);
	}
	//for (j=0;j<nl*nmax*kmax*mmax;j++){printf("%2d \t\t %19.12lE \n", j, J_dot_r[j]);}
	return(0);
}
#endif

#ifdef IS_DELTA_J
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

#ifdef IS_J2J_DOT
int main(){
	int j;
	double J_dot_sf[3], J_inital[3];
	double Delta_J_r_tidal[2], Delta_J_theta_tidal[2], Delta_J_phi_tidal[2];
	double angle_space, angle_step, L_dot[50], Kepler_torque[50];
	int nl, nmax, kmax, mmax, N_res;
	int n_res_inner, k_res_inner, m_res_inner;
	int n_res_outer, k_res_outer, m_res_outer;
	double apo_res, peri, incline, mass, spin, radius_outer, guess1, guess2;

	#if 0
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
	#endif
	
	printf("Enter central mass: ");
	scanf("%lf", &mass);

	printf("Enter spin parameter of BH: ");
	scanf("%lf", &spin);

	printf("Number of modes and max n,k,m: ");
	scanf("%i %i %i %i", &nl, &nmax, &kmax, &mmax);

	printf("Initial J components: ");
	scanf("%lf %lf %lf", &J_inital[0], &J_inital[1], &J_inital[2]);

	#if 0
	printf("Number of modes and mode vector for inner orbit: ");
	scanf("%i %i %i %i", &nl, &n_res_inner, &k_res_inner, &m_res_inner);

	printf("Number of resonance terms: ");
	scanf("%i", &N_res);

	printf("Mode vector for outer orbit: ");
	scanf("%i %i %i", &n_res_outer, &k_res_outer, &m_res_outer);
	#endif
	
	#if 0
	printf("About to compute tidal J_dots \n");
	J_dot_tidal(nl, N_res, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, apo_res, peri, radius_outer, incline, mass, spin, -M_PI/2., &J_dot_r_tidal, &J_dot_theta_tidal, &J_dot_phi_tidal);
	printf("J_dot_r_tidal = %lg \n", J_dot_r_tidal);
	printf("J_dot_theta_tidal = %lg \n", J_dot_theta_tidal);
	printf("J_dot_phi_tidal = %lg \n", J_dot_phi_tidal);
	printf("Keplerian J_dot_phi at resonance = %lg \n", J_dot_phi_Kepler(1.0, radius_outer, apo_res, peri, incline));
	#endif

	
	#if 0
	J_dot_selfforce(nl, nmax, kmax, mmax, apo_res, peri, radius_outer, incline, mass, spin, J_dot_sf);
	printf("J_dot_r_sf = %lg \n", J_dot_sf[0]);
	printf("J_dot_theta_sf = %lg \n", J_dot_sf[1]);
	printf("J_dot_phi_sf = %lg \n", J_dot_sf[2]);
	#endif

	J2Jdot(nl, nmax, kmax, mmax, J_inital, J_dot_sf, mass, spin);
	printf("J_dot components: %lg %lg %lg \n", J_dot_sf[0], J_dot_sf[1], J_dot_sf[2]);

	#if 0
	printf("About to start L_dot \n");
	for (j = 0; j < 50; j++){
		angle_space = 2 * M_PI / 50;
		angle_step = j * angle_space;
		J_dot_tidal(nl, N_res, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, apo_res, peri, radius_outer, incline, mass, spin, angle_step, J_dot_td);
		L_dot[j] = J_dot_td[2];
		Kepler_torque[j] = J_dot_phi_Kepler(1.0, radius_outer, apo_res, peri, incline, angle_step);
		printf("%lg %lg %lg \n", angle_step, L_dot[j], Kepler_torque[j]);
	}
	//for (j=0;j<nl*nmax*kmax*mmax;j++){printf("%2d \t\t %19.12lE \n", j, J_dot_r[j]);}
	#endif
	return(0);
}
#endif

#if IS_RK4_J_DOT
int main()
{
	double h=0.1, t, t0, t1, J_r_ini, J_theta_ini, J_phi_ini, *J_r_final, *J_theta_final, *J_phi_final;
	int i, n;
	//t0, t1, J_r_ini, J_theta_ini, J_phi_ini
	printf("Starting and ending time t0 t1, and inital Js: ");
	scanf("%lf %lf %lf %lf %lf", &t0, &t1, &J_r_ini, &J_theta_ini, &J_phi_ini);
	n = 1 + (t1 - t0)/h;
	printf("n=%i \n", n);
	J_r_final = (double *)malloc(sizeof(double) * n);
	J_theta_final = (double *)malloc(sizeof(double) * n);
	J_phi_final = (double *)malloc(sizeof(double) * n);
	printf("Inital Js are: %lg %lg %lg\n", J_r_ini, J_theta_ini, J_phi_ini);

	rk4_J2Jdot(h, t0, n, J_r_ini, J_theta_ini, J_phi_ini, J_r_final, J_theta_final, J_phi_final);
	/*for (y[0] = 1, i = 1; i < n; i++){y[i] = rk4(h, x0 + h * (i - 1), y[i-1]);}*/
		
 
	printf("t \t J_r \t J_theta \t J_phi \n------------\n");
	for (i = 0; i < n; i++) {
		t = t0 + h * i;
		printf("%g\t%g\t%g\t%g\n", t, J_r_final, J_theta_final, J_phi_final);
	}
 
 	free((char*)J_r_final);
 	free((char*)J_theta_final);
	free((char*)J_phi_final);
	return 0;

}
#endif