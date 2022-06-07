#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CKerr.h"

/*#include "kerrtraj.c"
#include "kerrmode.c"
#include "kerrgwem.c"*/

/* Find resonance location (the apocenter) given the pericenter and inclination 
* using the constraint that 
* (-2\Omega^{\phi} + 2\Omega^{\theta} + \Omega^{\r})_{inner} = -2\Omega^{\phi}_{outer} 
*/


/* Make function that converts the input data from 
* (apocenter (ra), pericenter (rp), inclination (I), central BH mass (M), spin of BH (a)) --> (E,Q,L)
* EQL[0] = E, EQL[1] = Q, EQL[2] = L
*/

double ra_rp_I2EQL(double ra, double *EQL, double rp, double I, double astar, double M) {
	int i;
	double z;
	double a;
	/*double EQL[3];*/
	double Coeff_rp[4];
	double Coeff_ra[4];
	double Delta_rp;
	double Delta_ra;
	double squared_prefactor_ratio_LE;
	double linear_prefactor_ratio_LE;
	double constant_prefactor_ratio_LE;
	double discriminant;
	double relevant_ratio;
	double denom_plus;
	double denom_minus;
	double which_E[2];
	double which_ratio[2];


	/* Defining the coefficients that will appear in the conic section expression for EL relation */
	a = astar*M;
	z = sin(I);
	Delta_rp = rp*rp - 2*M*rp + a*a;
	Delta_ra = ra*ra - 2*M*ra + a*a;
	Coeff_rp[0] = (rp*rp + a*a)*(rp*rp + a*a) - Delta_rp*a*a*(1-z*z); /* E^2 coeff. */
	Coeff_rp[1] = -2*(a*(rp*rp + a*a) - Delta_rp*a); /* E*L coeff. */
	Coeff_rp[2] = a*a - Delta_rp/(1 - z*z); /* L^2 coeff. */
	Coeff_rp[3] = Delta_rp*(rp*rp + a*a*z*z); /* Constant term */
	Coeff_ra[0] = (ra*ra + a*a)*(ra*ra + a*a) - Delta_ra*a*a*(1-z*z); /* E^2 coeff. */
	Coeff_ra[1] = -2*(a*(ra*ra + a*a) - Delta_ra*a); /* E*L coeff. */
	Coeff_ra[2] = a*a - Delta_ra/(1 - z*z); /* L^2 coeff. */
	Coeff_ra[3] = Delta_ra*(ra*ra + a*a*z*z); /* Constant term */

	/*printf("For pericenter stuff: \nA = %E \nB = %E \nC = %E \nD = %E \n", Coeff_rp[0], Coeff_rp[1], Coeff_rp[2], Coeff_rp[3]);
	printf("For apocenter stuff: \nA = %E \nB = %E \nC = %E \nD = %E \n", Coeff_ra[0], Coeff_ra[1], Coeff_ra[2], Coeff_ra[3]);*/

	/* Defining the quadratic equation for the ratio of L/E */
	squared_prefactor_ratio_LE = Coeff_rp[2]*Coeff_ra[3] - Coeff_ra[2]*Coeff_rp[3];
	linear_prefactor_ratio_LE = (Coeff_rp[1]*Coeff_ra[3] - Coeff_ra[1]*Coeff_rp[3]);
	constant_prefactor_ratio_LE = Coeff_rp[0]*Coeff_ra[3] - Coeff_ra[0]*Coeff_rp[3];
	discriminant = linear_prefactor_ratio_LE*linear_prefactor_ratio_LE - 4*squared_prefactor_ratio_LE*constant_prefactor_ratio_LE;

	/*printf("squared coeff.: %E, linear coeff.: %E, constant coeff.: %E, discriminant: %E \n", squared_prefactor_ratio_LE, linear_prefactor_ratio_LE, constant_prefactor_ratio_LE, discriminant);*/

	/* The two solutions from the ratio quadratic */
	which_ratio[0] = (-linear_prefactor_ratio_LE + sqrt(discriminant))/(2*squared_prefactor_ratio_LE);
	which_ratio[1] = (-linear_prefactor_ratio_LE - sqrt(discriminant))/(2*squared_prefactor_ratio_LE);

	/*printf("Value of the ratio L/E are: %lf and %lf \n", which_ratio[0], which_ratio[1]);*/

	/* Get the four values of E from the two values of ratio */
	denom_plus = Coeff_rp[2]*which_ratio[0]*which_ratio[0] + Coeff_rp[1]*which_ratio[0] + Coeff_rp[0];
	denom_minus = Coeff_rp[2]*which_ratio[1]*which_ratio[1] + Coeff_rp[1]*which_ratio[1] + Coeff_rp[0];

	which_E[0] = sqrt(Coeff_rp[3]/denom_plus);
	/*E1_plus = sqrt(Coeff_rp[3]/denom_plus);*/
	/*E1_minus = -sqrt(Coeff_rp[3]/denom_plus);*/
	which_E[1] = sqrt(Coeff_rp[3]/denom_minus);
	/*E2_plus = sqrt(Coeff_rp[3]/denom_minus);*/
	/*E2_minus = -sqrt(Coeff_rp[3]/denom_minus);*/

	/*printf("values of denom_plus = %E \nvalues of denom_minus = %E \n", denom_plus, denom_minus);*/
	/*printf("E1: %lf, E2: %lf, E3: %lf, E4: %lf \n", E1_plus, E1_minus, E2_plus, E2_minus);*/
	/*printf("E1: %f, E2: %f \n", which_E[0], which_E[1]);*/

	/* Now we make conditions to choose which value of E is the true E */
	/* -- FILL IN LOGIC STATEMENT TO PICK WHICH E IS THE TRUE E -- */
	for (i=0; i<2; i++){
		if(which_E[i]<1){
			EQL[0] = which_E[i];
			relevant_ratio = which_ratio[i];
			/*printf("The relevant ratio of L/E is: %lf \n The relevant E is: %lf \n", relevant_ratio, EQL[0]);*/
		}
		/*else{
		*	EQL[0] = E2_plus;
		*	relevant_ratio = ratio_LE_minus;
		*	printf("The relevant ratio of L/E is: %lf \nThe relevant E is: %lf \n", relevant_ratio, EQL[0]);
		}*/
	}

	/* Now we get L from the chosen E */
	EQL[2] = relevant_ratio*EQL[0];

	/* Finally, we get Q from the computed E and L */
	EQL[1] = z*z/(1 - z*z)*EQL[2]*EQL[2] - a*a*z*z*EQL[0]*EQL[0] + a*a*z*z;

	/*printf("For a given apo, peri, incline, spin parameter, and central mass, we get the values \nE = %lf \nQ = %lf \nL = %lf \n", EQL[0], EQL[1], EQL[2]);*/
	/*return(EQL);*/
}

/* Defining the outer angular frequency using the expression derived in Chandrasekhar Ch. 7 for direct orbits */
double Omega_outer_direct(double radius, double M, double spin){
	double u;
	double omegaouter;

	u = 1./radius;
	omegaouter = sqrt(M*u*u*u)/(1. + spin * M * sqrt(M*u*u*u));

	return(omegaouter);

}

/* Define the resonance function */
/*double resonance_function(double *Omega_inner, double radius, double M, double spin){
	double function;
	double omega_out;
	int n,k,m;


	printf("Enter resonance coefficients: ");
	scanf("%i %i %i", &n, &k, &m);

	omega_out = Omega_outer_direct(radius, M, spin);
	function = n*Omega_inner[0] + k*Omega_inner[1] + m*Omega_inner[2] - m*omega_out;

	return(function);

}*/

/* Converts apocenter, pericenter, and inclination to Omega resonance condition */

/* For outer body on a circular orbit, use this function. Only needs to be fed a mode vector for the inner orbit and the outer orbit radius. */
double ra_rp_I2Omega_OuterCirc(int n, int k, int m, double radius, double ra, double rp, double I, double astar, double M){
	double EQL[3];
	double J[3];
	double Minv[9];
	double Omega_inner[3];
	double omega_out;
	double function;

	ra_rp_I2EQL(ra, EQL, rp, I, astar, M);
	CKerr_EQL2J(EQL, J, M, astar, NULL);
	CKerr_Minverse(J, Minv, M, astar);
	CKerr_Minv2Omega(Minv, Omega_inner);
	omega_out = Omega_outer_direct(radius, M, astar);
	function = - n*Omega_inner[0] - k*Omega_inner[1] - m*Omega_inner[2] + m*omega_out;
	return(function);
}

/* If you want a both orbits to be generic, then give apocenter, pericenter, and inclination for both */
double ra_rp_I2Omega_generic(int n_inner, int k_inner, int m_inner, int n_outer, int k_outer, double ra_inner, double rp_inner, double I_inner, double ra_outer, double rp_outer, double I_outer, double astar, double M){
	double EQL_inner[3], EQL_outer[3];
	double J_inner[3], J_outer[3];
	double Minv_inner[9], Minv_outer[9];
	double Omega_inner[3];
	double Omega_outer[3];
	//double omega_out;
	double function;

	/* Inner Data */
	if(rp_inner == 0 && I_inner == 0){
		Omega_inner[0] = 0;
		Omega_inner[1] = 0;
		Omega_inner[3] = Omega_outer_direct(ra_inner, M, astar);
	}
	else{
		ra_rp_I2EQL(ra_inner, EQL_inner, rp_inner, I_inner, astar, M);
		CKerr_EQL2J(EQL_inner, J_inner, M, astar, NULL);
		CKerr_Minverse(J_inner, Minv_inner, M, astar);
		CKerr_Minv2Omega(Minv_inner, Omega_inner);
	}
	
	/* Outer Data */
	if(rp_outer == 0 && I_outer == 0){
		Omega_outer[0] = 0;
		Omega_outer[1] = 0;
		Omega_outer[3] = Omega_outer_direct(ra_outer, M, astar);
	}
	else{
		ra_rp_I2EQL(ra_outer, EQL_outer, rp_outer, I_outer, astar, M);
		CKerr_EQL2J(EQL_outer, J_outer, M, astar, NULL);
		CKerr_Minverse(J_outer, Minv_outer, M, astar);
		CKerr_Minv2Omega(Minv_outer, Omega_outer);
	}
	

	//omega_out = Omega_outer_direct(radius, M, astar);
	//function = n*Omega_inner[0] + k*Omega_inner[1] + m*Omega_inner[2] - m*omega_out;
	function = (n_outer*Omega_outer[0] - n_inner*Omega_inner[0]) + (k_outer*Omega_outer[1] - k_inner*Omega_inner[1]) + m_inner * (Omega_outer[2] - Omega_inner[2]);
	//function = n*Omega_inner[0] + k*Omega_inner[1] + m*Omega_inner[2];
	return(function);
	//function = resonance_function(Omega_inner, radius, M, astar);
}


/* Bisection method for finding the apocenter of the resonance condition given the outer body is on a circular equatorial orbit */
double find_resonance_apo_OuterCirc(int n, int k, int m, double radius, double guess1, double guess2, double rp, double I, double astar, double M){
	double mid_apo=0;
	double error=0;
	int i=0;
	double mold=0;


	if(ra_rp_I2Omega_OuterCirc(n, k, m, radius, guess1, rp, I, astar, M)*ra_rp_I2Omega_OuterCirc(n, k, m, radius, guess2, rp, I, astar, M)>0){
        	printf("Value at first and second guess: %lf, %lf \n", ra_rp_I2Omega_OuterCirc(n, k, m, radius, guess1, rp, I, astar, M), ra_rp_I2Omega_OuterCirc(n, k, m, radius, guess2, rp, I, astar, M));
        	printf("Invalid Interval Exit! \n");       //to test whether search interval
        	exit(1);                                                        //is okay or not
	}
	else if(ra_rp_I2Omega_OuterCirc(n, k, m, radius, guess1, rp, I, astar, M)==0 || ra_rp_I2Omega_OuterCirc(n, k, m, radius, guess2, rp, I, astar, M)==0){
        	printf("Value at first and second guess: %lf, %lf \n", ra_rp_I2Omega_OuterCirc(n, k, m, radius, guess1, rp, I, astar, M), ra_rp_I2Omega_OuterCirc(n, k, m, radius, guess2, rp, I, astar, M));
        	printf("Root is one of interval bounds. Root is %f\n",ra_rp_I2Omega_OuterCirc(n, k, m, radius, guess1, rp, I, astar, M)==0?guess1:guess2);
        	exit(0);
	}
	printf("Ite\ta\t\tb\t\tmid_apo\t\tf(mid_apo)\t\terror\n");
	do{
        	mold = mid_apo;
        	mid_apo = (guess1 + guess2)/2;
        	printf("%2d\t%4.6f\t%4.6f\t%4.6f\t%lg\t",i++,guess1,guess2,mid_apo,ra_rp_I2Omega_OuterCirc(n, k, m, radius, mid_apo, rp, I, astar, M));
        	if(ra_rp_I2Omega_OuterCirc(n, k, m, radius, mid_apo, rp, I, astar, M)==0){
                	printf("Root is %4.6f\n",mid_apo);
        	}else if ((ra_rp_I2Omega_OuterCirc(n, k, m, radius, guess1, rp, I, astar, M)*ra_rp_I2Omega_OuterCirc(n, k, m, radius, mid_apo, rp, I, astar, M))<0){
                	guess2 = mid_apo;
        	}else guess1 = mid_apo;
        	error=fabs((mid_apo - mold)/mid_apo);
        	if(i==1){
                	printf("----\n");
        	}else printf("%4.6f\n",error);
	}while(error>1e-4);
	printf("Approximate Root is %4.6f \n",mid_apo);
	return(mid_apo);
	}
	