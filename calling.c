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

#ifdef IS_J_DOT
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
	// sscanf(argv[4], "%lg", &J_outer[0]);
	// sscanf(argv[5], "%lg", &J_outer[1]);
	// sscanf(argv[6], "%lg", &J_outer[2]);
	// sscanf(argv[7], "%d", &n_res_inner);
	// sscanf(argv[8], "%d", &k_res_inner);
	// sscanf(argv[9], "%d", &m_res_inner);
	// sscanf(argv[10], "%d", &n_res_outer);
	// sscanf(argv[11], "%d", &k_res_outer);
	// sscanf(argv[12], "%d", &m_res_outer);
	// sscanf(argv[13], "%lg", &angle_torus);

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

	// CKerr_J2EQL(J_outer, EQL_outer, mass, spin);
	// CKerr_EQL2J(EQL_outer, J_outer, mass, spin, anc_outer);

	I_inner = anc_inner[0];
	rp_inner = anc_inner[1];
	ra_inner = anc_inner[2];

	// I_outer = anc_outer[0];
	// rp_outer = anc_outer[1];
	// ra_outer = anc_outer[2];


	// printf("Location on torus (in radians): ");
	// scanf("%lf", &angle_torus);



	#if 0
	printf("Number of modes and max n,k,m: ");
	scanf("%i %i %i %i", &nl, &nmax, &kmax, &mmax);
	#endif

	// printf("Mode vector for inner orbit: ");
	// scanf("%i %i %i", &n_res_inner, &k_res_inner, &m_res_inner);

	// printf("Number of resonance terms: ");
	// scanf("%i", &N_res);

	// printf("Mode vector for outer orbit: ");
	// scanf("%i %i %i", &n_res_outer, &k_res_outer, &m_res_outer);

	double e_inner = (ra_inner - rp_inner) / (ra_inner + rp_inner);

	printf("I_inner, rp_inner, ra_inner, eccentricity: %lg %lg %lg %lg \n", I_inner, rp_inner, ra_inner, e_inner);
	// printf("I_outer, rp_outer, ra_outer: %lg %lg %lg \n", I_outer, rp_outer, ra_outer);
	// printf("Inner and Outer EQL: %lg %lg %lg, %lg %lg %lg \n", EQL_inner[0], EQL_inner[1], EQL_inner[2], EQL_outer[0], EQL_outer[1], EQL_outer[2]);
	printf("EQL: %lg %lg %lg \n", EQL_inner[0], EQL_inner[1], EQL_inner[2]);
	
	#if 0
	printf("About to compute tidal J_dots \n");
	J_dot_tidal(nl, N_res, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, ra_inner, rp_inner, radius_outer, I_inner, ra_outer, rp_outer, I_outer, mass, spin, angle_torus, mu_outer, J_dot_td);
	printf("J_dot_r_tidal = %lg \n", J_dot_td[0]);
	printf("J_dot_theta_tidal = %lg \n", J_dot_td[1]);
	printf("J_dot_phi_tidal = %lg \n", J_dot_td[2]);
	// printf("Keplerian J_dot_phi at resonance = %lg \n", J_dot_phi_Kepler(1.0, radius_outer, apo_res, peri, incline));
	#endif
	

	
	#if 0
	printf("About to compute self-force J_dots \n");
	J_dot_selfforce(nl_self, nmax, kmax, mmax, ra_inner, rp_inner, radius_outer, I_inner, mass, spin, J_dot_sf);
	printf("J_dot_r_sf = %lg \n", J_dot_sf[0]);
	printf("J_dot_theta_sf = %lg \n", J_dot_sf[1]);
	printf("J_dot_phi_sf = %lg \n", J_dot_sf[2]);
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

#ifdef IS_DELTA_J2

typedef struct 
{
  // members for the student's type, name, age and average
  int  data_col0, data_col1, data_col8, data_col9, data_col10, data_col11, data_col12;
  double data_col2, data_col3, data_col4, data_col5, data_col6, data_col7, data_col13, data_col14, data_col15;
  double data_col16, data_col17, data_col18, data_col19, data_col20, data_col21, data_col22, data_col23, data_col24;

} data_vals;

//Filename wants 24 columns delimted by spaces
int readtxt(char FILENAME[], data_vals *data, long *number_row){
	// file pointer variable for accessing the file
  FILE *file;
  
  // attempt to open file.txt in read mode to read the file contents
  file = fopen(FILENAME, "r"); 
  
  // if the file failed to open, exit with an error message and status
  if (file == NULL)
  {
    printf("Error opening file.\n");
    return 1;
  }
  
  // array of structs for storing the data from the file
  //data_vals data_cols[N_LINES_MAX];
  data_vals *data_cols = data;
  
  // read will be used to ensure each line/record is read correctly
  int read = 0;
  
  // records will keep track of the number of Student records read from the file
  int records = 0;

  // read all records from the file and store them into the students array
  do 
  {
    // Read a line/record from the file with the above format, notice in 
    // particular how we read in the student's name with %49[^,] which matches
    // up to 49 characters NOT including the comma (so it will stop matching 
    // at the next comma).  The name member can store 50 characters, so 
    // factoring in the NULL terminator this is the maximum amount of characters
    // we can read in for a number.  fscanf() will return the number of values 
    // it was able to read successfully which we expect to be 4, and we store 
    // that into read.
    //
	// TODO: Rename the data_cols[].(USEFUL NAME)
    read = fscanf(file,
                  "%d %d %lg %lg %lg %lg %lg %lg %d %d %d %d %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
                  &data_cols[records].data_col0,
				  &data_cols[records].data_col1,
				  &data_cols[records].data_col2,
				  &data_cols[records].data_col3,
				  &data_cols[records].data_col4,
				  &data_cols[records].data_col5,
				  &data_cols[records].data_col6,
				  &data_cols[records].data_col7,
				  &data_cols[records].data_col8,
				  &data_cols[records].data_col9,
				  &data_cols[records].data_col10,
				  &data_cols[records].data_col11,
				  &data_cols[records].data_col12,
				  &data_cols[records].data_col13,
				  &data_cols[records].data_col14,
				  &data_cols[records].data_col15,
				  &data_cols[records].data_col16,
				  &data_cols[records].data_col17,
				  &data_cols[records].data_col18,
				  &data_cols[records].data_col19,
				  &data_cols[records].data_col20,
				  &data_cols[records].data_col21,
				  &data_cols[records].data_col22,
				  &data_cols[records].data_col23,
				  &data_cols[records].data_col24); 
    
    // if fscanf read 25 values from the file then we've successfully read 
    // in another record
	// Number of columns stored in "read"
    if (read == 25) records++;
    
    // The only time that fscanf should NOT read 4 values from the file is 
    // when we've reached the end of the file, so if fscanf did not read in 
    // exactly 4 values and we're not at the end of the file, there has been
    // an error (likely due to an incorrect file format) and so we exit with 
    // an error message and status.
    if (read != 25 && !feof(file) || records == N_LINES_MAX)
    {
      printf("File format incorrect.\n");
      return 1;
    }
    
    // if there was an error reading from the file exit with an error message 
    // and status
    if (ferror(file))
    {
      printf("Error reading file.\n");
      return 1;
    }

  } while (!feof(file));

  // close the file as we are done working with it
  fclose(file);
  
  // print out the number of records (number of rows in file ) read
  printf("\n%d records read.\n\n", records);
  *number_row = records;

  return(0);
}

int main(){
	int j;
	long i;
	double Delta_J_r_tidal, Delta_J_theta_tidal, Delta_J_phi_tidal;
	double Delta_J_tidal_value[3];
	int nl, nmax, kmax, mmax;
	int n_res_inner, k_res_inner, m_res_inner;
	int n_res_outer, k_res_outer, m_res_outer;
	double apo_res, peri, incline, mass, spin, radius_outer, guess1, guess2;
	double apo_inner, I_inner, peri_inner, theta_res_F;
	double apo_outer, I_outer, peri_outer;
	data_vals *data; //data_vals is defined as a data type, so we call the pointer "data" that data type
	long MAX_ROW; //Define the actual number of rows in the specific data file
	double ang_accel, Gamma, mu = 1e-5; //Value of angular acceleration (including body mass)


	printf("Enter central BH mass: ");
	scanf("%lf", &mass);
	printf("Enter spin parameter of BH: ");
	scanf("%lf", &spin);

	// printf("Data for inner body: \n");
	// printf("Enter pericenter: ");
	// scanf("%lf", &peri_inner);
	// printf("Enter inlincation angle (radians): ");
	// scanf("%lf", &I_inner);
	// printf("Enter corresponding apocenter at resonance: ");
	// scanf("%lf", &apo_inner);

	// printf("Data for outer body: \n");
	// printf("Enter pericenter: ");
	// scanf("%lf", &peri_outer);
	// printf("Enter inlincation angle (radians): ");
	// scanf("%lf", &I_outer);
	// printf("Enter corresponding apocenter at resonance: ");
	// scanf("%lf", &apo_outer);
	// printf("Enter outer radius: ");
	// scanf("%lf", &radius_outer);
	printf("Enter corresponding resonance angle: ");
	scanf("%lf", &theta_res_F);
	printf("Number of modes: ");
	scanf("%i", &nl);

	// //printf("Number of modes and max n,k,m: ");
	// //scanf("%i %i %i %i", &nl, &nmax, &kmax, &mmax);

	// printf("Number of modes and mode vector for inner orbit: ");
	// scanf("%i %i %i %i", &nl, &n_res_inner, &k_res_inner, &m_res_inner);

	// printf("Mode vector for outer orbit: ");
	// scanf("%i %i %i", &n_res_outer, &k_res_outer, &m_res_outer);

	data = (data_vals*)malloc((size_t)(N_LINES_MAX*sizeof(data_vals)));
	readtxt("gamma_vals_proper.txt", data, &MAX_ROW);
	// i = 0;
	// while (i < 10){
    // 	printf("%d %lg %d %lg\n", 
    //        	data[i].data_col1, 
    //        	data[i].data_col13,
    //        	data[i].data_col12,
    //        	data[i].data_col24);
 	// 	 printf("\n");
	// 	i++;
	// 	}
	// free((char*)data);
	//return(0);
	//TODO: Replace the appropriate arguments in Delta_J_tidal2 with the data[i].(whatever) that is appropriate for that value in the text file
	printf("i \t system label \t n_inner \t k_inner \t n_outer \t k_outer \t m \t Gamma \t Delta_J_r \t Delta_J_theta \t Delta_J_phi \n------------\n");
	for(i = 0; i < MAX_ROW; i++){

		//ang_accel = data[i].data_col5 * data[i].data_col7; 

		Delta_J_tidal2(nl, data[i].data_col9, data[i].data_col11, data[i].data_col10, data[i].data_col12, data[i].data_col8, data[i].data_col8, data[i].data_col21, data[i].data_col20, 0, data[i].data_col19, data[i].data_col24, data[i].data_col23, data[i].data_col22, mass, spin, theta_res_F, data[i].data_col7, data[i].data_col6, Delta_J_tidal_value);
		printf("%i \t %i \t %i \t %i \t %i \t %i \t %i \t %lg \t %lg \t %lg \t %lg \n", i, data[i].data_col1, data[i].data_col9, data[i].data_col10, data[i].data_col11, data[i].data_col12, data[i].data_col8, data[i].data_col7, Delta_J_tidal_value[0], Delta_J_tidal_value[1], Delta_J_tidal_value[2]);

	}
	// Delta_J_tidal2(nl, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, apo_inner, peri_inner, radius_outer, I_inner, apo_outer, peri_outer, I_outer, mass, spin, theta_res_F, ang_accel, Delta_J_tidal_value);
	// printf("Delta_J_r_tidal = %lg \n", Delta_J_tidal_value[0]);
	// printf("Delta_J_theta_tidal = %lg \n", Delta_J_tidal_value[1]);
	// printf("Delta_J_phi_tidal = %lg \n", Delta_J_tidal_value[2]);
	free((char*)data);
	return(0);
}
#endif

#ifdef IS_READ_TXT
typedef struct 
{
  // members for the student's type, name, age and average
  int data_col1, data_col8, data_col9, data_col10, data_col11, data_col12;
  double data_col2, data_col3, data_col4, data_col5, data_col6, data_col7, data_col13, data_col14, data_col15;
  double data_col16, data_col17, data_col18, data_col19, data_col20, data_col21, data_col22, data_col23, data_col24;

} data_vals;

//Filename wants 24 columns delimted by spaces
int readtxt(char FILENAME, struct (*data)){
	// file pointer variable for accessing the file
  FILE *file;
  
  // attempt to open file.txt in read mode to read the file contents
  file = fopen(FILENAME, "r"); 
  
  // if the file failed to open, exit with an error message and status
  if (file == NULL)
  {
    printf("Error opening file.\n");
    return 1;
  }
  
  // array of structs for storing the Student data from the file
  data_vals data_cols[N_LINES_MAX];
  
  // read will be used to ensure each line/record is read correctly
  int read = 0;
  
  // records will keep track of the number of Student records read from the file
  int records = 0;

  // read all records from the file and store them into the students array
  do 
  {
    // Read a line/record from the file with the above format, notice in 
    // particular how we read in the student's name with %49[^,] which matches
    // up to 49 characters NOT including the comma (so it will stop matching 
    // at the next comma).  The name member can store 50 characters, so 
    // factoring in the NULL terminator this is the maximum amount of characters
    // we can read in for a number.  fscanf() will return the number of values 
    // it was able to read successfully which we expect to be 4, and we store 
    // that into read.
    //
	// TODO: Rename the data_cols[].(USEFUL NAME)
    read = fscanf(file,
                  "%d %lg %lg %lg %lg %lg %lg %d %d %d %d %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
                  &data_cols[records].data_col1,
				  &data_cols[records].data_col2,
				  &data_cols[records].data_col3,
				  &data_cols[records].data_col4,
				  &data_cols[records].data_col5,
				  &data_cols[records].data_col6,
				  &data_cols[records].data_col7,
				  &data_cols[records].data_col8,
				  &data_cols[records].data_col9,
				  &data_cols[records].data_col10,
				  &data_cols[records].data_col11,
				  &data_cols[records].data_col12,
				  &data_cols[records].data_col13,
				  &data_cols[records].data_col14,
				  &data_cols[records].data_col15,
				  &data_cols[records].data_col16,
				  &data_cols[records].data_col17,
				  &data_cols[records].data_col18,
				  &data_cols[records].data_col19,
				  &data_cols[records].data_col20,
				  &data_cols[records].data_col21,
				  &data_cols[records].data_col22,
				  &data_cols[records].data_col23,
				  &data_cols[records].data_col24); 
    
    // if fscanf read 24 values from the file then we've successfully read 
    // in another record
	// TODO: Will add another column, then change to 25
	//printf("%d %d %d\n", read, feof(file), records);
    if (read == 24) records++;
    
    // The only time that fscanf should NOT read 4 values from the file is 
    // when we've reached the end of the file, so if fscanf did not read in 
    // exactly 4 values and we're not at the end of the file, there has been
    // an error (likely due to an incorrect file format) and so we exit with 
    // an error message and status.
    if (read != 24 && !feof(file))
    {
      printf("File format incorrect.\n");
      return 1;
    }
    
    // if there was an error reading from the file exit with an error message 
    // and status
    if (ferror(file))
    {
      printf("Error reading file.\n");
      return 1;
    }

  } while (!feof(file));

  // close the file as we are done working with it
  fclose(file);
  
  // print out the number of records read
  printf("\n%d records read.\n\n", records);

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

#if IS_RK4_J_DOT_norestart
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


	rk4_J2Jdot_restartII(t0, t, i_start, n, J_r_ini, J_theta_ini, J_phi_ini, J_r_final, J_theta_final, J_phi_final, RESTART, mu_body, M, astar);
	
		 
 
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
	rk4_J2Jdot_restartII(t_start, t, i_start, n, J_r_ini, J_theta_ini, J_phi_ini, J_r_final, J_theta_final, J_phi_final, RESTART, mu_body, M, astar);
		 
 
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

