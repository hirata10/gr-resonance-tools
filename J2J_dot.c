#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CKerr.h"
#include "globalpars_c.h"

/* This will take a set of Js (pointer) for an orbit and return the time derivative, J_dot (pointer) */
/* J_initial[0] = J_r_inital, J_initial[1] = J_theta_inital, J_initial[2] = J_phi_initial */
int J2Jdot(int nl, int nmax, int kmax, int mmax, double *J_initial, double *J_dot_sf, double M, double astar){
    double EQL[3], J[3], ancillary[3];
    double radius_outer = 1.;

    CKerr_J2EQL(J_initial, EQL, M, astar);
    CKerr_EQL2J(EQL, J, M, astar, ancillary);
    J_dot_selfforce(nl, nmax, kmax, mmax, ancillary[2], ancillary[1], radius_outer, ancillary[0], M, astar, J_dot_sf);
    return(0);
}

/* Component version of J2Jdot function: takes in individual components to J_inital and returns individual components to J_dot */
int J2Jdot_component(int nl, int nmax, int kmax, int mmax, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_dot_r, double *J_dot_theta, double *J_dot_phi, double M, double astar){
    double EQL[3], J[3], ancillary[3];
    double J_inital[3], J_dot_sf[3];
    double radius_outer = 1.;
    int flag;

    J_inital[0] = J_r_ini;
    J_inital[1] = J_theta_ini;
    J_inital[2] = J_phi_ini;
    flag = CKerr_J2EQL(J_inital, EQL, M, astar);
    if (flag == 0)
    {
        return(1);
    }
    
    CKerr_EQL2J(EQL, J, M, astar, ancillary);
    J_dot_selfforce(nl, nmax, kmax, mmax, ancillary[2], ancillary[1], radius_outer, ancillary[0], M, astar, J_dot_sf);
    *J_dot_r = J_dot_sf[0];
    *J_dot_theta = J_dot_sf[1];
    *J_dot_phi = J_dot_sf[2];
    return(0);
}

/* This takes a step size (dt), an inital J_inital components, and return arrays of size n for each J components at each time step */
void rk4_J2Jdot(double *t, int n, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_r_final, double *J_theta_final, double *J_phi_final, FILE *fptr, double mu_body, double M, double astar){
    int i, j;
    int nl = GLOBALPAR_nl_self, nmax = GLOBALPAR_nmax, kmax = GLOBALPAR_kmax, mmax = GLOBALPAR_mmax;
    double dt, inspiral_scale = 0.01;
    double a[3];
    double J_Minv[3];
    double Minv[9], Omega[3];
    double k1r, k2r, k3r, k4r;
    double k1theta, k2theta, k3theta, k4theta;
    double k1phi, k2phi, k3phi, k4phi;
    double J_dot_r, J_dot_theta, J_dot_phi;
    double J_r_final_temp, J_theta_final_temp, J_phi_final_temp;
    int Totflag;
    //double time_factor;
    // double *k1r, *k2r, *k3r, *k4r;
    // double *k1theta, *k2theta, *k3theta, *k4theta;
    // double *k1phi, *k2phi, *k3phi, *k4phi;
    // double *J_dot_r0, *J_dot_theta0, *J_dot_phi0;
    // double *t;
    double t0;

    //J_Minv = (double *)malloc(sizeof(double) * (n));

     /* Define time step array, starts at inital time (input) */
    // t = (double *)malloc(sizeof(double) * (n));
    // t[0] = t0;

    t0 = t[0];
    

    /* Define the array for the final solution, they should start at the inital conditions given */
    J_r_final[0] = J_r_ini;
    J_theta_final[0] = J_theta_ini;
    J_phi_final[0] = J_phi_ini;

    /* Array to convert initial J to initial Omega */
    J_Minv[0] = J_r_ini;
    J_Minv[1] = J_theta_ini;
    J_Minv[2] = J_phi_ini;


    CKerr_Minverse(J_Minv, Minv, M, astar);
    CKerr_Minv2Omega(Minv, Omega);


    /*
    J_dot_r = (double *)malloc(sizeof(double) * n);
    J_dot_theta = (double *)malloc(sizeof(double) * n);
    J_dot_phi = (double *)malloc(sizeof(double) * n);
    */

    /* Create slopes for each direction r,theta,phi */
    #if 0
    k1r = (double *)malloc(sizeof(double) * n);
    k2r = (double *)malloc(sizeof(double) * n);
    k3r = (double *)malloc(sizeof(double) * n);
    k4r = (double *)malloc(sizeof(double) * n);

    k1theta = (double *)malloc(sizeof(double) * n);
    k2theta = (double *)malloc(sizeof(double) * n);
    k3theta = (double *)malloc(sizeof(double) * n);
    k4theta = (double *)malloc(sizeof(double) * n);

    k1phi = (double *)malloc(sizeof(double) * n);
    k2phi = (double *)malloc(sizeof(double) * n);
    k3phi = (double *)malloc(sizeof(double) * n);
    k4phi = (double *)malloc(sizeof(double) * n);
    #endif
    printf("i \t t \t J_r \t J_theta \t J_phi \t Omega_r \t Omega_theta \t Omega_phi \t dt \n------------\n");
    fprintf(fptr, "i \t t \t J_r \t J_theta \t J_phi \t Omega_r \t Omega_theta \t Omega_phi \t dt \n------------\n");
    printf("%i \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \n", 0, t[0], J_r_final[0], J_theta_final[0], J_phi_final[0], Omega[0], Omega[1], Omega[2], 0.0);
    fprintf(fptr, "%i \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \n", 0, t[0], J_r_final[0], J_theta_final[0], J_phi_final[0], Omega[0], Omega[1], Omega[2], 00.0);
    fflush(fptr);

    for (i = 1; i < (n); i++){
        
        /* k1_{r,theta,phi} will come from the same call of J2Jdot_component since no changes start yet */
        Totflag += J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        a[0] = inspiral_scale * fabs(J_r_final[i-1]/J_dot_r);
        a[1] = inspiral_scale * fabs(J_theta_final[i-1]/J_dot_theta);
        a[2] = inspiral_scale * fabs(J_phi_final[i-1]/J_dot_phi);
        dt = a[0]; //variable time step based on J/\dot{J}
        for (j = 1; j < 3; j++){dt = dt < a[j] ? dt:a[j];}
        //dt = dt * mu_body; //time step scaled with mass of orbiting body    
        k1r = dt * (J_dot_r);
        k1theta = dt * (J_dot_theta);
        k1phi = dt * (J_dot_phi);

        #if 1
        printf("k1r, k1theta, k1phi, dt, Totflag: %lg %lg %lg %lg %i \n", k1r, k1theta, k1phi, dt, Totflag);
        fprintf(fptr, "k1r, k1theta, k1phi, dt, Totflag: %lg %lg %lg %lg %i \n", k1r, k1theta, k1phi, dt, Totflag);
        fflush(fptr);
        #endif

        Totflag += J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k1r/2., J_theta_final[i-1] + k1theta/2., J_phi_final[i-1] + k1phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2r= dt * (J_dot_r);
       //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k1theta/2., J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2theta = dt * (J_dot_theta);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k1phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2phi = dt * (J_dot_phi);

        #if 1
        printf("k2r, k2theta, k2phi, Totflag: %lg %lg %lg %i \n", k2r, k2theta, k2phi, Totflag);
        fprintf(fptr, "k2r, k2theta, k2phi, Totflag: %lg %lg %lg %i \n", k2r, k2theta, k2phi, Totflag);
        fflush(fptr);
        #endif

        Totflag += J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k2r/2., J_theta_final[i-1] + k2theta/2., J_phi_final[i-1] + k2phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3r = dt * (J_dot_r);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k2theta/2., J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3theta = dt * (J_dot_theta);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k2phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3phi = dt * (J_dot_phi);

        #if 1
        printf("k3r, k3theta, k3phi, Totflag: %lg %lg %lg %i \n", k3r, k3theta, k3phi, Totflag);
        fprintf(fptr, "k3r, k3theta, k3phi, Totflag: %lg %lg %lg %i \n", k3r, k3theta, k3phi, Totflag);
        fflush(fptr);
        #endif

        Totflag += J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k3r, J_theta_final[i-1] + k3theta, J_phi_final[i-1] + k4phi, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4r = dt * (J_dot_r);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k3theta, J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4theta = dt * (J_dot_theta);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k3phi, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4phi = dt * (J_dot_phi);

        #if 1
        printf("k4r, k4theta, k4phi, Totflag: %lg %lg %lg %i \n", k4r, k4theta, k4phi, Totflag);
        fprintf(fptr, "k4r, k4theta, k4phi, Totflag: %lg %lg %lg %i \n", k4r, k4theta, k4phi, Totflag);
        fflush(fptr);
        #endif

        t[i] = t[i-1] + dt;
        // J_r_final[i] = J_r_final[i-1] + (k1r + 2. * k2r + 2. * k3r + k4r) / 6.;
        // J_theta_final[i] = J_theta_final[i-1] + (k1theta + 2. * k2theta + 2. * k3theta + k4theta) / 6.;
        // J_phi_final[i] = J_phi_final[i-1] + (k1phi + 2. * k2phi + 2. * k3phi + k4phi) / 6.;

        J_r_final_temp = J_r_final[i-1] + (k1r + 2. * k2r + 2. * k3r + k4r) / 6.;
        J_theta_final_temp = J_theta_final[i-1] + (k1theta + 2. * k2theta + 2. * k3theta + k4theta) / 6.;
        J_phi_final_temp = J_phi_final[i-1] + (k1phi + 2. * k2phi + 2. * k3phi + k4phi) / 6.;

        #if 1
        printf("J_r_temp, J_theta_temp, J_phi_temp: %lg %lg %lg \n", J_r_final_temp, J_theta_final_temp, J_phi_final_temp);
        fprintf(fptr, "J_r_temp, J_theta_temp, J_phi_temp: %lg %lg %lg \n", J_r_final_temp, J_theta_final_temp, J_phi_final_temp);
        fflush(fptr);
        #endif

        /* Condition for plunge i.e. gives NaN --- Update J_i pointers */
        if(isnan(J_r_final_temp) || isnan(J_theta_final_temp) || isnan(J_phi_final_temp) || Totflag > 0){
            printf("It plunged.\n");
            break;
        }
        else{
            J_r_final[i] = J_r_final_temp;
            J_theta_final[i] = J_theta_final_temp;
            J_phi_final[i] = J_phi_final_temp;

            /* Array to compute Omegas */
            J_Minv[0] = J_r_final_temp;
            J_Minv[1] = J_theta_final_temp;
            J_Minv[2] = J_phi_final_temp;
        }

        /* Convert J_i_final[i] to the Omega_i */
        CKerr_Minverse(J_Minv, Minv, M, astar);
        CKerr_Minv2Omega(Minv, Omega);

        //printf("%i\t %g\t \n", i, x[i]);
		printf("%i \t%f \t%e \t%e \t%e \t%e \t%e \t%e \t%e \n", i, t[i], J_r_final[i], J_theta_final[i], J_phi_final[i], Omega[0], Omega[1], Omega[2], dt);
        fprintf(fptr, "%i \t%f \t%e \t%e \t%e \t%e \t%e \t%e \t%e \n", i, t[i], J_r_final[i], J_theta_final[i], J_phi_final[i], Omega[0], Omega[1], Omega[2], dt);
        fflush(fptr);
    }
    // free((char*)t);
    #if 0
    for (i = 1; i < (n+1); i++){
        /* k1_{r,theta,phi} will come from the same call of J2Jdot_component since no changes start yet */
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k1r[i-1] = dt * (J_dot_r);
        k1theta[i-1] = dt * (J_dot_theta);
        k1phi[i-1] = dt * (J_dot_phi);

        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k1r[i-1]/2, J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2r[i-1] = dt * (J_dot_r);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k1theta[i-1]/2, J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2theta[i-1] = dt * (J_dot_theta);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k1phi[i-1]/2, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2phi[i-1] = dt * (J_dot_phi);
           
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k2r[i-1]/2, J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3r[i-1] = dt * (J_dot_r);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k2theta[i-1]/2, J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3theta[i-1] = dt * (J_dot_theta);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k2phi[i-1]/2, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3phi[i-1] = dt * (J_dot_phi);

        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k3r[i-1], J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4r[i-1] = dt * (J_dot_r);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k3theta[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4theta[i-1] = dt * (J_dot_theta);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k3phi[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4phi[i-1] = dt * (J_dot_phi);

        t[i] = t[i-1] + dt;
        J_r_final[i] = J_r_final[i-1] + (k1r[i-1] + 2 * k2r[i-1] + 2 * k3r[i-1] + k4r[i-1]) / 6;
        J_theta_final[i] = J_theta_final[i-1] + (k1theta[i-1] + 2 * k2theta[i-1] + 2 * k3theta[i-1] + k4theta[i-1]) / 6;
        J_phi_final[i] = J_phi_final[i-1] + (k1phi[i-1] + 2 * k2phi[i-1] + 2 * k3phi[i-1] + k4phi[i-1]) / 6;
        //printf("%i\t %g\t \n", i, x[i]);
		printf("%i \t%lf \t%lf \t%lf \t%lf \t%lf \t %lf \n", i-1, t[i-1], J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1], k1r[i-1], k1theta[i-1]);
    }
    free((char*)t);
    free((char*)k1r);
    free((char*)k1theta);
    free((char*)k1phi);
    free((char*)k2r);
    free((char*)k2theta);
    free((char*)k2phi);
    free((char*)k3r);
    free((char*)k3theta);
    free((char*)k3phi);
    free((char*)k4r);
    free((char*)k4theta);
    free((char*)k4phi);
    #endif
}

/* This is for restarting the RK4, an inital J_inital components, and return arrays of size n for each J components at each time step */
// TODO: Modify to take the stating index as an argument
void rk4_J2Jdot_restart(double *t, int i_start, int n, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_r_final, double *J_theta_final, double *J_phi_final, int RESTART, double mu_body, double M, double astar){
    int i, j;
    int nl = GLOBALPAR_nl_self, nmax = GLOBALPAR_nmax, kmax = GLOBALPAR_kmax, mmax = GLOBALPAR_mmax;
    double dt, inspiral_scale = 0.01;
    double a[3];
    double J_Minv[3];
    double Minv[9], Omega[3];
    double k1r, k2r, k3r, k4r;
    double k1theta, k2theta, k3theta, k4theta;
    double k1phi, k2phi, k3phi, k4phi;
    double J_dot_r, J_dot_theta, J_dot_phi;
    double J_r_final_temp, J_theta_final_temp, J_phi_final_temp;
    int Totflag;
    //double time_factor;
    // double *k1r, *k2r, *k3r, *k4r;
    // double *k1theta, *k2theta, *k3theta, *k4theta;
    // double *k1phi, *k2phi, *k3phi, *k4phi;
    // double *J_dot_r0, *J_dot_theta0, *J_dot_phi0;
    // double *t;
    double t0;

    //J_Minv = (double *)malloc(sizeof(double) * (n));

     /* Define time step array, starts at inital time (input) */
    // t = (double *)malloc(sizeof(double) * (n));
    // t[0] = t0;

    t0 = t[0];
    

    /* Define the array for the final solution, they should start at the inital conditions given */
    J_r_final[0] = J_r_ini;
    J_theta_final[0] = J_theta_ini;
    J_phi_final[0] = J_phi_ini;

    /* Array to convert initial J to initial Omega */
    J_Minv[0] = J_r_ini;
    J_Minv[1] = J_theta_ini;
    J_Minv[2] = J_phi_ini;


    CKerr_Minverse(J_Minv, Minv, M, astar);
    CKerr_Minv2Omega(Minv, Omega);


    /*
    J_dot_r = (double *)malloc(sizeof(double) * n);
    J_dot_theta = (double *)malloc(sizeof(double) * n);
    J_dot_phi = (double *)malloc(sizeof(double) * n);
    */

    /* Create slopes for each direction r,theta,phi */
    #if 0
    k1r = (double *)malloc(sizeof(double) * n);
    k2r = (double *)malloc(sizeof(double) * n);
    k3r = (double *)malloc(sizeof(double) * n);
    k4r = (double *)malloc(sizeof(double) * n);

    k1theta = (double *)malloc(sizeof(double) * n);
    k2theta = (double *)malloc(sizeof(double) * n);
    k3theta = (double *)malloc(sizeof(double) * n);
    k4theta = (double *)malloc(sizeof(double) * n);

    k1phi = (double *)malloc(sizeof(double) * n);
    k2phi = (double *)malloc(sizeof(double) * n);
    k3phi = (double *)malloc(sizeof(double) * n);
    k4phi = (double *)malloc(sizeof(double) * n);
    #endif
    if (RESTART == 1){
        printf("i \t t \t J_r \t J_theta \t J_phi \t Omega_r \t Omega_theta \t Omega_phi \t dt \n------------\n");
        printf("%i \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \n", 0, t[0], J_r_final[0], J_theta_final[0], J_phi_final[0], Omega[0], Omega[1], Omega[2], 0.0);
        fflush(stdout);
    }

    for (i = i_start; i < (n); i++){
        
        /* k1_{r,theta,phi} will come from the same call of J2Jdot_component since no changes start yet */
        Totflag += J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        a[0] = inspiral_scale * fabs(J_r_final[i-1]/J_dot_r);
        a[1] = inspiral_scale * fabs(J_theta_final[i-1]/J_dot_theta);
        a[2] = inspiral_scale * fabs(J_phi_final[i-1]/J_dot_phi);
        dt = a[0]; //variable time step based on J/\dot{J}
        for (j = 1; j < 3; j++){dt = dt < a[j] ? dt:a[j];}
        //dt = dt * mu_body; //time step scaled with mass of orbiting body    
        k1r = dt * (J_dot_r);
        k1theta = dt * (J_dot_theta);
        k1phi = dt * (J_dot_phi);

        #if 1
        printf("k1r, k1theta, k1phi, dt, Totflag: %lg %lg %lg %lg %i \n", k1r, k1theta, k1phi, dt, Totflag);
        fflush(stdout);
        #endif

        Totflag += J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k1r/2., J_theta_final[i-1] + k1theta/2., J_phi_final[i-1] + k1phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2r= dt * (J_dot_r);
       //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k1theta/2., J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2theta = dt * (J_dot_theta);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k1phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2phi = dt * (J_dot_phi);

        #if 1
        printf("k2r, k2theta, k2phi, Totflag: %lg %lg %lg %i \n", k2r, k2theta, k2phi, Totflag);
        fflush(stdout);
        #endif

        Totflag += J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k2r/2., J_theta_final[i-1] + k2theta/2., J_phi_final[i-1] + k2phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3r = dt * (J_dot_r);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k2theta/2., J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3theta = dt * (J_dot_theta);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k2phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3phi = dt * (J_dot_phi);

        #if 1
        printf("k3r, k3theta, k3phi, Totflag: %lg %lg %lg %i \n", k3r, k3theta, k3phi, Totflag);
        fflush(stdout);
        #endif

        Totflag += J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k3r, J_theta_final[i-1] + k3theta, J_phi_final[i-1] + k4phi, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4r = dt * (J_dot_r);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k3theta, J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4theta = dt * (J_dot_theta);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k3phi, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4phi = dt * (J_dot_phi);

        #if 1
        printf("k4r, k4theta, k4phi, Totflag: %lg %lg %lg %i \n", k4r, k4theta, k4phi, Totflag);
        fflush(stdout);
        #endif

        t[i] = t[i-1] + dt;

        // J_r_final[i] = J_r_final[i-1] + (k1r + 2. * k2r + 2. * k3r + k4r) / 6.;
        // J_theta_final[i] = J_theta_final[i-1] + (k1theta + 2. * k2theta + 2. * k3theta + k4theta) / 6.;
        // J_phi_final[i] = J_phi_final[i-1] + (k1phi + 2. * k2phi + 2. * k3phi + k4phi) / 6.;

        J_r_final_temp = J_r_final[i-1] + (k1r + 2. * k2r + 2. * k3r + k4r) / 6.;
        J_theta_final_temp = J_theta_final[i-1] + (k1theta + 2. * k2theta + 2. * k3theta + k4theta) / 6.;
        J_phi_final_temp = J_phi_final[i-1] + (k1phi + 2. * k2phi + 2. * k3phi + k4phi) / 6.;

        #if 1
        printf("t[%i], J_r_temp, J_theta_temp, J_phi_temp: %lg %lg %lg %lg \n", i, t[i], J_r_final_temp, J_theta_final_temp, J_phi_final_temp);
        fflush(stdout);
        #endif

        /* Condition for plunge i.e. gives NaN --- Update J_i pointers */
        if(isnan(J_r_final_temp) || isnan(J_theta_final_temp) || isnan(J_phi_final_temp) || Totflag > 0){
            printf("It plunged.\n");
            break;
        }
        else{
            J_r_final[i] = J_r_final_temp;
            J_theta_final[i] = J_theta_final_temp;
            J_phi_final[i] = J_phi_final_temp;

            /* Array to compute Omegas */
            J_Minv[0] = J_r_final_temp;
            J_Minv[1] = J_theta_final_temp;
            J_Minv[2] = J_phi_final_temp;
        }

        /* Convert J_i_final[i] to the Omega_i */
        CKerr_Minverse(J_Minv, Minv, M, astar);
        CKerr_Minv2Omega(Minv, Omega);

        //printf("%i\t %g\t \n", i, x[i]);
		printf("%i \t%f \t%e \t%e \t%e \t%e \t%e \t%e \t%e \n", i, t[i], J_r_final[i], J_theta_final[i], J_phi_final[i], Omega[0], Omega[1], Omega[2], dt);
        fflush(stdout);
    }
    // free((char*)t);
    #if 0
    for (i = 1; i < (n+1); i++){
        /* k1_{r,theta,phi} will come from the same call of J2Jdot_component since no changes start yet */
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k1r[i-1] = dt * (J_dot_r);
        k1theta[i-1] = dt * (J_dot_theta);
        k1phi[i-1] = dt * (J_dot_phi);

        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k1r[i-1]/2, J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2r[i-1] = dt * (J_dot_r);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k1theta[i-1]/2, J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2theta[i-1] = dt * (J_dot_theta);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k1phi[i-1]/2, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2phi[i-1] = dt * (J_dot_phi);
           
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k2r[i-1]/2, J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3r[i-1] = dt * (J_dot_r);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k2theta[i-1]/2, J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3theta[i-1] = dt * (J_dot_theta);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k2phi[i-1]/2, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3phi[i-1] = dt * (J_dot_phi);

        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k3r[i-1], J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4r[i-1] = dt * (J_dot_r);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k3theta[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4theta[i-1] = dt * (J_dot_theta);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k3phi[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4phi[i-1] = dt * (J_dot_phi);

        t[i] = t[i-1] + dt;
        J_r_final[i] = J_r_final[i-1] + (k1r[i-1] + 2 * k2r[i-1] + 2 * k3r[i-1] + k4r[i-1]) / 6;
        J_theta_final[i] = J_theta_final[i-1] + (k1theta[i-1] + 2 * k2theta[i-1] + 2 * k3theta[i-1] + k4theta[i-1]) / 6;
        J_phi_final[i] = J_phi_final[i-1] + (k1phi[i-1] + 2 * k2phi[i-1] + 2 * k3phi[i-1] + k4phi[i-1]) / 6;
        //printf("%i\t %g\t \n", i, x[i]);
		printf("%i \t%lf \t%lf \t%lf \t%lf \t%lf \t %lf \n", i-1, t[i-1], J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1], k1r[i-1], k1theta[i-1]);
    }
    free((char*)t);
    free((char*)k1r);
    free((char*)k1theta);
    free((char*)k1phi);
    free((char*)k2r);
    free((char*)k2theta);
    free((char*)k2phi);
    free((char*)k3r);
    free((char*)k3theta);
    free((char*)k3phi);
    free((char*)k4r);
    free((char*)k4theta);
    free((char*)k4phi);
    #endif
}

/* This is for restarting the RK4, an inital J_inital components, and return arrays of size n for each J components at each time step */
void rk4_J2Jdot_restartII(double t_start, double *t, int i_start, int n, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_r_final, double *J_theta_final, double *J_phi_final, int RESTART, double mu_body, double M, double astar){
    int i, j;
    int nl = GLOBALPAR_nl_self, nmax = GLOBALPAR_nmax, kmax = GLOBALPAR_kmax, mmax = GLOBALPAR_mmax;
    double dt, inspiral_scale = 0.01;
    double a[3];
    double J_Minv[3];
    double Minv[9], Omega[3];
    double k1r, k2r, k3r, k4r;
    double k1theta, k2theta, k3theta, k4theta;
    double k1phi, k2phi, k3phi, k4phi;
    double J_dot_r, J_dot_theta, J_dot_phi;
    double J_r_final_temp, J_theta_final_temp, J_phi_final_temp;
    int Totflag = 0, flag;
    double J_check[3], EQL_check[3];
    //double time_factor;
    // double *k1r, *k2r, *k3r, *k4r;
    // double *k1theta, *k2theta, *k3theta, *k4theta;
    // double *k1phi, *k2phi, *k3phi, *k4phi;
    // double *J_dot_r0, *J_dot_theta0, *J_dot_phi0;
    // double *t;
    // double t0;

    //J_Minv = (double *)malloc(sizeof(double) * (n));

     /* Define time step array, starts at inital time (input) */
    // t = (double *)malloc(sizeof(double) * (n));
    // t[0] = t0;


    /*
    J_dot_r = (double *)malloc(sizeof(double) * n);
    J_dot_theta = (double *)malloc(sizeof(double) * n);
    J_dot_phi = (double *)malloc(sizeof(double) * n);
    */

    /* Create slopes for each direction r,theta,phi */
    #if 0
    k1r = (double *)malloc(sizeof(double) * n);
    k2r = (double *)malloc(sizeof(double) * n);
    k3r = (double *)malloc(sizeof(double) * n);
    k4r = (double *)malloc(sizeof(double) * n);

    k1theta = (double *)malloc(sizeof(double) * n);
    k2theta = (double *)malloc(sizeof(double) * n);
    k3theta = (double *)malloc(sizeof(double) * n);
    k4theta = (double *)malloc(sizeof(double) * n);

    k1phi = (double *)malloc(sizeof(double) * n);
    k2phi = (double *)malloc(sizeof(double) * n);
    k3phi = (double *)malloc(sizeof(double) * n);
    k4phi = (double *)malloc(sizeof(double) * n);
    #endif
    if (RESTART == 1){
        t[0] = t_start;
    

        /* Define the array for the final solution, they should start at the inital conditions given */
        J_r_final[0] = J_r_ini;
        J_theta_final[0] = J_theta_ini;
        J_phi_final[0] = J_phi_ini;

        /* Array to convert initial J to initial Omega */
        J_Minv[0] = J_r_ini;
        J_Minv[1] = J_theta_ini;
        J_Minv[2] = J_phi_ini;


        CKerr_Minverse(J_Minv, Minv, M, astar);
        CKerr_Minv2Omega(Minv, Omega);
        printf("i \t t \t J_r \t J_theta \t J_phi \t Omega_r \t Omega_theta \t Omega_phi \t dt \n------------\n");
        printf("%i \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \n", 0, t[0], J_r_final[0], J_theta_final[0], J_phi_final[0], Omega[0], Omega[1], Omega[2], 0.0);
        fflush(stdout);
    }

    for (i = i_start; i < (n); i++){
        
        t[i-1] = t_start;
        J_r_final[i-1] = J_r_ini;
        J_theta_final[i-1] = J_theta_ini;
        J_phi_final[i-1] = J_phi_ini;

        #if 0
        printf("Starting time and J: %lg %lg %lg %lg \n", t[i-1], J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1]);
        fflush(stdout);
        #endif
        
        /* k1_{r,theta,phi} will come from the same call of J2Jdot_component since no changes start yet */
        Totflag += J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        a[0] = inspiral_scale * fabs(J_r_final[i-1]/J_dot_r);
        a[1] = inspiral_scale * fabs(J_theta_final[i-1]/J_dot_theta);
        a[2] = inspiral_scale * fabs(J_phi_final[i-1]/J_dot_phi);
        dt = a[0]; //variable time step based on J/\dot{J}
        for (j = 1; j < 3; j++){dt = dt < a[j] ? dt:a[j];}
        //dt = dt * mu_body; //time step scaled with mass of orbiting body    
        k1r = dt * (J_dot_r);
        k1theta = dt * (J_dot_theta);
        k1phi = dt * (J_dot_phi);

        t[i] = t[i-1] + dt; // Time stepped after computing appropriate jump
        #if 0
        printf("Next time step: %lf \n", t[i]);
        fflush(stdout);
        #endif

        #if 0
        printf("k1r, k1theta, k1phi, dt, Totflag: %lg %lg %lg %lf %i \n", k1r, k1theta, k1phi, dt, Totflag);
        fflush(stdout);
        #endif

        Totflag += J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k1r/2., J_theta_final[i-1] + k1theta/2., J_phi_final[i-1] + k1phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2r= dt * (J_dot_r);
       //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k1theta/2., J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2theta = dt * (J_dot_theta);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k1phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2phi = dt * (J_dot_phi);

        #if 0
        printf("k2r, k2theta, k2phi, Totflag: %lg %lg %lg %i \n", k2r, k2theta, k2phi, Totflag);
        fflush(stdout);
        #endif

        Totflag += J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k2r/2., J_theta_final[i-1] + k2theta/2., J_phi_final[i-1] + k2phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3r = dt * (J_dot_r);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k2theta/2., J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3theta = dt * (J_dot_theta);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k2phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3phi = dt * (J_dot_phi);

        #if 0
        printf("k3r, k3theta, k3phi, Totflag: %lg %lg %lg %i \n", k3r, k3theta, k3phi, Totflag);
        fflush(stdout);
        #endif

        Totflag += J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k3r, J_theta_final[i-1] + k3theta, J_phi_final[i-1] + k4phi, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4r = dt * (J_dot_r);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k3theta, J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4theta = dt * (J_dot_theta);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k3phi, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4phi = dt * (J_dot_phi);

        #if 0
        printf("k4r, k4theta, k4phi, Totflag: %lg %lg %lg %i \n", k4r, k4theta, k4phi, Totflag);
        fflush(stdout);
        #endif

        // J_r_final[i] = J_r_final[i-1] + (k1r + 2. * k2r + 2. * k3r + k4r) / 6.;
        // J_theta_final[i] = J_theta_final[i-1] + (k1theta + 2. * k2theta + 2. * k3theta + k4theta) / 6.;
        // J_phi_final[i] = J_phi_final[i-1] + (k1phi + 2. * k2phi + 2. * k3phi + k4phi) / 6.;

        J_r_final_temp = J_r_final[i-1] + (k1r + 2. * k2r + 2. * k3r + k4r) / 6.;
        J_theta_final_temp = J_theta_final[i-1] + (k1theta + 2. * k2theta + 2. * k3theta + k4theta) / 6.;
        J_phi_final_temp = J_phi_final[i-1] + (k1phi + 2. * k2phi + 2. * k3phi + k4phi) / 6.;

        J_check[0] = J_r_final_temp;
        J_check[1] = J_theta_final_temp;
        J_check[2] = J_phi_final_temp;


        flag = CKerr_J2EQL(J_check, EQL_check, M, astar); // If flag == 0, then the orbit is unbound after the RK4 procedure

        #if 0
        printf("t[%i], J_r_temp, J_theta_temp, J_phi_temp: %f %lg %lg %lg \n", i, t[i], J_r_final_temp, J_theta_final_temp, J_phi_final_temp);
        fflush(stdout);
        #endif

        /* Condition for plunge i.e. gives NaN --- Update J_i pointers */
        if(isnan(J_r_final_temp) || isnan(J_theta_final_temp) || isnan(J_phi_final_temp) || Totflag > 0 || flag == 0){
            printf("It plunged.\n");
            break;
        }
        else{
            J_r_final[i] = J_r_final_temp;
            J_theta_final[i] = J_theta_final_temp;
            J_phi_final[i] = J_phi_final_temp;

            /* Array to compute Omegas */
            J_Minv[0] = J_r_final_temp;
            J_Minv[1] = J_theta_final_temp;
            J_Minv[2] = J_phi_final_temp;
        }

        /* Convert J_i_final[i] to the Omega_i */
        CKerr_Minverse(J_Minv, Minv, M, astar);
        CKerr_Minv2Omega(Minv, Omega);

        //printf("%i\t %g\t \n", i, x[i]);
		printf("%i \t%f \t%e \t%e \t%e \t%e \t%e \t%e \t%e \n", i, t[i], J_r_final[i], J_theta_final[i], J_phi_final[i], Omega[0], Omega[1], Omega[2], dt);
        fflush(stdout);

        // Reset values at the beginning to start process for the next iteration
        t_start = t[i];
        J_r_ini = J_r_final[i];
        J_theta_ini = J_theta_final[i];
        J_phi_ini = J_phi_final[i];

    }
    // free((char*)t);
    #if 0
    for (i = 1; i < (n+1); i++){
        /* k1_{r,theta,phi} will come from the same call of J2Jdot_component since no changes start yet */
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k1r[i-1] = dt * (J_dot_r);
        k1theta[i-1] = dt * (J_dot_theta);
        k1phi[i-1] = dt * (J_dot_phi);

        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k1r[i-1]/2, J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2r[i-1] = dt * (J_dot_r);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k1theta[i-1]/2, J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2theta[i-1] = dt * (J_dot_theta);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k1phi[i-1]/2, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2phi[i-1] = dt * (J_dot_phi);
           
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k2r[i-1]/2, J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3r[i-1] = dt * (J_dot_r);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k2theta[i-1]/2, J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3theta[i-1] = dt * (J_dot_theta);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k2phi[i-1]/2, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3phi[i-1] = dt * (J_dot_phi);

        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k3r[i-1], J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4r[i-1] = dt * (J_dot_r);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k3theta[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4theta[i-1] = dt * (J_dot_theta);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k3phi[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4phi[i-1] = dt * (J_dot_phi);

        t[i] = t[i-1] + dt;
        J_r_final[i] = J_r_final[i-1] + (k1r[i-1] + 2 * k2r[i-1] + 2 * k3r[i-1] + k4r[i-1]) / 6;
        J_theta_final[i] = J_theta_final[i-1] + (k1theta[i-1] + 2 * k2theta[i-1] + 2 * k3theta[i-1] + k4theta[i-1]) / 6;
        J_phi_final[i] = J_phi_final[i-1] + (k1phi[i-1] + 2 * k2phi[i-1] + 2 * k3phi[i-1] + k4phi[i-1]) / 6;
        //printf("%i\t %g\t \n", i, x[i]);
		printf("%i \t%lf \t%lf \t%lf \t%lf \t%lf \t %lf \n", i-1, t[i-1], J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1], k1r[i-1], k1theta[i-1]);
    }
    free((char*)t);
    free((char*)k1r);
    free((char*)k1theta);
    free((char*)k1phi);
    free((char*)k2r);
    free((char*)k2theta);
    free((char*)k2phi);
    free((char*)k3r);
    free((char*)k3theta);
    free((char*)k3phi);
    free((char*)k4r);
    free((char*)k4theta);
    free((char*)k4phi);
    #endif
}