#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CKerr.h"
#include "globalpars_c.h"

/* LEGACY CODE (outdated; kept for reference): J2Jdot, J2Jdot_component, rk4_J2Jdot and rk4_J2Jdot_restart.
 * They are built on the serial J_dot_selfforce (J_dot.c), which has no mass ratio. The two RK4 routines take mu_body
 * but do not use it, so their time variable is mu_body * t, not physical time.
 * Neither the pipeline nor the Makefile's test tools use them: step 4 uses J2Jdot_component_openmp and
 * rk4_J2Jdot_restart_openmp (both in this file), which take mu_body and advance physical time. */

/* LEGACY (see the note above) */
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

/* LEGACY (see the note at the top of this file) */
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

/* Component version of J2Jdot function: takes in individual components to J_inital and returns individual components to J_dot_sf using openmp functions */
/* mu_body = (body mass)/M: J_dot is d(J/mu_body)/dt in physical time, so the RK4 that calls this advances physical time */
int J2Jdot_component_openmp(int nl, int nmax, int kmax, int mmax, int Nr, int Nt, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_dot_r, double *J_dot_theta, double *J_dot_phi, double mu_body, double M, double astar){
    double EQL[3], J[3], ancillary[3];
    double J_inital[3], J_dot_sf[3];
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
    J_dot_selfforce_openmp(nl, nmax, kmax, mmax, Nr, Nt, ancillary[2], ancillary[1], ancillary[0], mu_body, M, astar, J_dot_sf);
    *J_dot_r = J_dot_sf[0];
    *J_dot_theta = J_dot_sf[1];
    *J_dot_phi = J_dot_sf[2];
    return(0);
}

/* LEGACY (see the note at the top of this file): serial RK4, its time is mu_body * t */
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

        Totflag += J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k3r, J_theta_final[i-1] + k3theta, J_phi_final[i-1] + k3phi, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
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

/* LEGACY (see the note at the top of this file): serial RK4 with restart, its time is mu_body * t */
/* This is for restarting the RK4, an inital J_inital components, and return arrays of size n for each J components at each time step */
void rk4_J2Jdot_restart(double t_start, double *t, int i_start, int n, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_r_final, double *J_theta_final, double *J_phi_final, int RESTART, double mu_body, double M, double astar){
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

        Totflag += J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k3r, J_theta_final[i-1] + k3theta, J_phi_final[i-1] + k3phi, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
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
		printf("%i \t%f \t%.15g \t%.15g \t%.15g \t%.15g \t%.15g \t%.15g \t%g \n", i, t[i], J_r_final[i], J_theta_final[i], J_phi_final[i], Omega[0], Omega[1], Omega[2], dt);
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

/* This is for restarting the RK4, an inital J_inital components, and return arrays of size n for each J components at each time step */
void rk4_J2Jdot_restart_openmp(double t_start, double *t, int i_start, int n, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_r_final, double *J_theta_final, double *J_phi_final, int RESTART, double mu_body, int Nr, int Nt, double M, double astar){
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
        printf("%i \t%f \t%.15g \t%.15g \t%.15g \t%.15g \t%.15g \t%.15g \t%g \n", 0, t[0], J_r_final[0], J_theta_final[0], J_phi_final[0], Omega[0], Omega[1], Omega[2], 0.0);
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
        Totflag += J2Jdot_component_openmp(nl, nmax, kmax, mmax, Nr, Nt, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, mu_body, M, astar);
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

        Totflag += J2Jdot_component_openmp(nl, nmax, kmax, mmax, Nr, Nt, J_r_final[i-1] + k1r/2., J_theta_final[i-1] + k1theta/2., J_phi_final[i-1] + k1phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, mu_body, M, astar);
        k2r= dt * (J_dot_r);
       //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k1theta/2., J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2theta = dt * (J_dot_theta);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k1phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2phi = dt * (J_dot_phi);

        #if 0
        printf("k2r, k2theta, k2phi, Totflag: %lg %lg %lg %i \n", k2r, k2theta, k2phi, Totflag);
        fflush(stdout);
        #endif

        Totflag += J2Jdot_component_openmp(nl, nmax, kmax, mmax, Nr, Nt, J_r_final[i-1] + k2r/2., J_theta_final[i-1] + k2theta/2., J_phi_final[i-1] + k2phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, mu_body, M, astar);
        k3r = dt * (J_dot_r);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k2theta/2., J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3theta = dt * (J_dot_theta);
        //J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k2phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3phi = dt * (J_dot_phi);

        #if 0
        printf("k3r, k3theta, k3phi, Totflag: %lg %lg %lg %i \n", k3r, k3theta, k3phi, Totflag);
        fflush(stdout);
        #endif

        Totflag += J2Jdot_component_openmp(nl, nmax, kmax, mmax, Nr, Nt, J_r_final[i-1] + k3r, J_theta_final[i-1] + k3theta, J_phi_final[i-1] + k3phi, &J_dot_r, &J_dot_theta, &J_dot_phi, mu_body, M, astar);
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
		printf("%i \t%f \t%.15g \t%.15g \t%.15g \t%.15g \t%.15g \t%.15g \t%g \n", i, t[i], J_r_final[i], J_theta_final[i], J_phi_final[i], Omega[0], Omega[1], Omega[2], dt);
        fflush(stdout);

        // Reset values at the beginning to start process for the next iteration
        t_start = t[i];
        J_r_ini = J_r_final[i];
        J_theta_ini = J_theta_final[i];
        J_phi_ini = J_phi_final[i];

    }
}

/*
------------------------------------------------------------
Actions and angles of both bodies: one RK4 step (resonance-aware evolution, src/python/J_evolve_resonant.py)
------------------------------------------------------------
Equations, for each body, in physical (Boyer-Lindquist) time:
    d(J_i)/dt   = J_dot_sf,i(J)     self-force rate, J2Jdot_component_openmp (the right-hand side of step 4)
    d(psi_i)/dt = Omega_i(J)        fundamental frequencies, CKerr_Minverse + CKerr_Minv2Omega
The angles psi = (psi_r, psi_theta, psi_phi) are measured from the torus origin of CKerr_TorusOrigin (kerrtraj.c):
psi = 0 is the point of the orbit at pericentre (u_r = 0), crossing the equator with u_theta > 0, at phi = 0.
The resonant angle of a mode pair is Theta = N_outer.psi_outer - N_inner.psi_inner (Eq. (5) of arXiv:2207.07733);
its value at a resonance crossing is the theta_res_F of Delta_J_from_amplitudes (Delta_J.c).
Only the frequencies drive the angles: the order-mu corrections of each body's own self-force to d(psi)/dt are not included.
*/

/* Right-hand side for one body at actions J[3]: J_dot[3] = d(J/mu_body)/dt from the self-force, Omega[3] = d(psi)/dt. */
/* Self-force sum and grid as in rk4_J2Jdot_restart_openmp. Returns 0, or 1 if J is not a bound orbit ("plunge"). */
int J_psi_rhs_openmp(const double *J, double mu_body, int Nr, int Nt, double M, double astar, double *J_dot, double *Omega){
    double J_copy[3], Minv[9];
    int j;

    if (J2Jdot_component_openmp(GLOBALPAR_nl_self, GLOBALPAR_nmax, GLOBALPAR_kmax, GLOBALPAR_mmax, Nr, Nt, J[0], J[1], J[2], &J_dot[0], &J_dot[1], &J_dot[2], mu_body, M, astar) != 0) return(1);
    for (j = 0; j < 3; j++){J_copy[j] = J[j];}
    if (CKerr_Minverse(J_copy, Minv, M, astar) == 0) return(1);
    CKerr_Minv2Omega(Minv, Omega);
    return(0);
}

/* One classical RK4 step of y[12] = (J_inner[3], psi_inner[3], J_outer[3], psi_outer[3]). */
/* Time step: dt = inspiral_scale * min |J_j / J_dot_j| over the actions of both bodies, with J_dot at the start of the step */
/* (the rule and the inspiral_scale of rk4_J2Jdot_restart_openmp), but not longer than dt_max when dt_max > 0. An action with */
/* J_j = 0 or J_dot_j = 0 is left out of the minimum. */
/* The stages and the sum of the J part are written as in rk4_J2Jdot_restart_openmp (k = dt * J_dot), so with dt_max <= 0 the */
/* inner body's J follows step 4; the angles do not feed back into J. */
/* Outputs: y_new[12]; *dt = the step taken; Omega_new[6] = frequencies of the inner, then the outer body at y_new. */
/* Returns 0, or 1 if a stage or the new state is not a bound orbit for either body ("plunge"); the outputs are then not valid. */
int rk4_J_psi_step_openmp(const double *y, double dt_max, double mu_inner, double mu_outer, int Nr, int Nt, double M, double astar, double *y_new, double *dt, double *Omega_new){
    double inspiral_scale = 0.01;
    double rate[12], k[4][12], y_stage[12], J_copy[3], Minv[9];
    double stage_weight[4] = {0., 0.5, 0.5, 1.}; /* stage s is evaluated at y + stage_weight[s] * k[s-1] */
    double mu[2], h = 0., a;
    int s, b, j;

    mu[0] = mu_inner; mu[1] = mu_outer;

    for (s = 0; s < 4; s++){
        for (j = 0; j < 12; j++){y_stage[j] = s == 0 ? y[j] : y[j] + stage_weight[s] * k[s-1][j];}

        /* d(y)/dt at this stage: (J_dot, Omega) of the inner body, then of the outer body */
        for (b = 0; b < 2; b++){
            if (J_psi_rhs_openmp(y_stage + 6*b, mu[b], Nr, Nt, M, astar, rate + 6*b, rate + 6*b + 3) != 0) return(1);
        }

        if (s == 0){
            /* variable time step based on J/\dot{J} at the start of the step */
            for (b = 0; b < 2; b++){
                for (j = 0; j < 3; j++){
                    if (y[6*b + j] != 0. && rate[6*b + j] != 0.){
                        a = inspiral_scale * fabs(y[6*b + j] / rate[6*b + j]);
                        if (h == 0. || a < h) h = a;
                    }
                }
            }
            if (dt_max > 0. && dt_max < h) h = dt_max;
            if (!(h > 0.) || !isfinite(h)) return(1);
        }

        for (j = 0; j < 12; j++){k[s][j] = h * rate[j];}
    }

    for (j = 0; j < 12; j++){y_new[j] = y[j] + (k[0][j] + 2. * k[1][j] + 2. * k[2][j] + k[3][j]) / 6.;}

    /* Frequencies at the new state. CKerr_Minverse returns 0 if the new J is not a bound orbit (the plunge test of step 4) */
    for (b = 0; b < 2; b++){
        for (j = 0; j < 3; j++){J_copy[j] = y_new[6*b + j];}
        if (!isfinite(J_copy[0]) || !isfinite(J_copy[1]) || !isfinite(J_copy[2])) return(1);
        if (CKerr_Minverse(J_copy, Minv, M, astar) == 0) return(1);
        CKerr_Minv2Omega(Minv, Omega_new + 3*b);
    }

    *dt = h;
    return(0);
}
