#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CKerr.h"

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

    J_inital[0] = J_r_ini;
    J_inital[1] = J_theta_ini;
    J_inital[2] = J_phi_ini;
    CKerr_J2EQL(J_inital, EQL, M, astar);
    CKerr_EQL2J(EQL, J, M, astar, ancillary);
    J_dot_selfforce(nl, nmax, kmax, mmax, ancillary[2], ancillary[1], radius_outer, ancillary[0], M, astar, J_dot_sf);
    *J_dot_r = J_dot_sf[0];
    *J_dot_theta = J_dot_sf[1];
    *J_dot_phi = J_dot_sf[2];
    return(0);
}

/* This takes a step size (dt), an inital J_inital components, and return arrays of size n for each J components at each time step */
void rk4_J2Jdot(double dt, double t0, int n, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_r_final, double *J_theta_final, double *J_phi_final){
    int i;
    int nl = 1, nmax = 1, kmax = 1, mmax = 1;
    double M = 1, astar = 0.5;
    double k1r, k2r, k3r, k4r;
    double k1theta, k2theta, k3theta, k4theta;
    double k1phi, k2phi, k3phi, k4phi;
    double J_dot_r, J_dot_theta, J_dot_phi;
    // double *k1r, *k2r, *k3r, *k4r;
    // double *k1theta, *k2theta, *k3theta, *k4theta;
    // double *k1phi, *k2phi, *k3phi, *k4phi;
    // double *J_dot_r, *J_dot_theta, *J_dot_phi;
    // double *J_dot_r, *J_dot_theta, *J_dot_phi;
    // double *J_dot_r, *J_dot_theta, *J_dot_phi;
    // double *J_dot_r, *J_dot_theta, *J_dot_phi;
    // double *J_dot_r, *J_dot_theta, *J_dot_phi;
    // double *J_dot_r, *J_dot_theta, *J_dot_phi;
    // double *J_dot_r, *J_dot_theta, *J_dot_phi;
    // double *J_dot_r, *J_dot_theta, *J_dot_phi;
    // double *J_dot_r0, *J_dot_theta0, *J_dot_phi0;
    double *t;

    /* Define the array for the final solution, they should start at the inital conditions given */
    J_r_final[0] = J_r_ini;
    J_theta_final[0] = J_theta_ini;
    J_phi_final[0] = J_phi_ini;

    /* Define time step array, starts at inital time (input) */
    t = (double *)malloc(sizeof(double) * n);
    t[0] = t0;

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
    printf("i \t t \t J_r \t J_theta \t J_phi \n------------\n");
    //printf("%i \t%lf \t%lf \t%lf \t%lf \t%lf \t %lf \n", 0, t[0], J_r_final[0], J_theta_final[0], J_phi_final[0], k1r[0], k1theta[0]);
    
    for (i = 1; i < (n + 1); i++){
        /* k1_{r,theta,phi} will come from the same call of J2Jdot_component since no changes start yet */
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k1r = dt * (J_dot_r);
        k1theta = dt * (J_dot_theta);
        k1phi = dt * (J_dot_phi);

        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k1r/2., J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2r= dt * (J_dot_r);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k1theta/2., J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2theta = dt * (J_dot_theta);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k1phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k2phi = dt * (J_dot_phi);
           
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k2r/2., J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3r = dt * (J_dot_r);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k2theta/2., J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3theta = dt * (J_dot_theta);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k2phi/2., &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k3phi = dt * (J_dot_phi);

        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1] + k3r, J_theta_final[i-1], J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4r = dt * (J_dot_r);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1] + k3theta, J_phi_final[i-1], &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4theta = dt * (J_dot_theta);
        J2Jdot_component(nl, nmax, kmax, mmax, J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1] + k3phi, &J_dot_r, &J_dot_theta, &J_dot_phi, M, astar);
        k4phi = dt * (J_dot_phi);

        t[i] = t[i-1] + dt;
        J_r_final[i] = J_r_final[i-1] + (k1r + 2. * k2r + 2. * k3r + k4r) / 6.;
        J_theta_final[i] = J_theta_final[i-1] + (k1theta + 2. * k2theta + 2. * k3theta + k4theta) / 6.;
        J_phi_final[i] = J_phi_final[i-1] + (k1phi + 2. * k2phi + 2. * k3phi + k4phi) / 6.;
        //printf("%i\t %g\t \n", i, x[i]);
		printf("%i \t%f \t%e \t%e \t%e \t%e \t %e \n", i-1, t[i-1], J_r_final[i-1], J_theta_final[i-1], J_phi_final[i-1], k1r, k1theta);
    }
    free((char*)t);
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