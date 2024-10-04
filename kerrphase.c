#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "CKerr.h"
#include "globalpars_c.h"

// REMEMBER WRAPPER!!!

/* SECTION 1: CONTOUR INTEGRAL */



// Function to calculate chi //
double calculate_chi(double a, double E, double L, double Q) {
    return (a * a * (1 - E * E)) / (L * L + Q * Q);
}

// Function to calculate r // 
double calculate_r(double L, double Q) {
    return (L * L - Q) / (L * L + Q);
}

// Function to calculate z_+ and z_- //
void calculate_z_roots(double chi, double r, complex double *z_minus, complex double *z_plus) {
    double discriminant = sqrt(1 + chi * chi + 2 * r * chi);
    double z_minus_squared = (1 + chi - discriminant) / (2 * chi);
    double z_plus_squared = (1 + chi + discriminant) / (2 * chi);

    *z_minus = sqrt(z_minus_squared);
    *z_plus = sqrt(z_plus_squared);
}

// Polynomial input to the integrals (equation A10 in the paper.) //
complex double P(complex double z, double a, double E, double L, double Q) {
    return (Q * (1 - z * z)) - (a * a * (1 - E * E) * z * z * (1 - z * z) * (1 - z * z)) - (L * L * z * z);
}

// Square root of polynomial function above //
complex double sqrt_P(complex double z, complex double z_minus, complex double z_plus, double a, double E, double L, double Q) {
    return I * a * sqrt(1 - E * E) * z * csqrt(1 - ((z_minus * z_minus)/(z * z))) * csqrt((z_plus * z_plus) - (z * z));
}

// uz value present in integrals //
complex double uz(complex double z, double a, double E, double L, double Q){
    return csqrt(P(z, a, E, L, Q)) / (1 - z * z);
}

// Function to evaluate equation A7 in the paper //
complex double d2JdE2(complex double z, complex double z_minus, complex double z_plus, double a, double E, double L, double Q) {
    complex double numerator = (P(z, a, E, L, Q) - (a * a * z * z * E * E * (1 - z * z))) * (a * a * z * z);
    complex double denominator = 2 * M_PI * sqrt_P(z, z_minus, z_plus, a, E, L, Q) * P(z, a, E, L, Q);

    return numerator / denominator;
}

// Function to evaluate equation A8 in the paper //
complex double d2JdL2(complex double z, complex double z_minus, complex double z_plus, double a, double E, double L, double Q) {
    complex double numerator = - (P(z, a, E, L, Q) + (z * z * L * L)) * (z * z);
    complex double denominator = 2 * M_PI * sqrt_P(z, z_minus, z_plus, a, E, L, Q) * P(z, a, E, L, Q) * (1 - z * z);

    return numerator / denominator;
}

// Function to evaluate equation A8 in the paper //
complex double d2JdEdL(complex double z, complex double z_minus, complex double z_plus, double a, double E, double L, double Q) {
    complex double numerator = a * a * z * z * z * z * E * L;
    complex double denominator = 2 * M_PI * sqrt_P(z, z_minus, z_plus, a, E, L, Q) * P(z, a, E, L, Q);

    return numerator / denominator;
}

// Function to perform the contour integration //
complex double contour_integral(complex double (*f)(complex double, complex double, complex double, double, double, double, double),
                                            complex double z_start, double radius, int num_points,
                                            double a, double E, double L, double Q,
                                            complex double z_minus, complex double z_plus) {
    complex double integral = 0.0 + 0.0 * I;
    double theta;
    complex double z;
    double dtheta = 2 * M_PI / num_points;  // angle step

    for (int i = 0; i < num_points; i++) {
        theta = i * dtheta;
        z = radius * (cos(theta) + I * sin(theta));  // point on the contour
        // Now passing a, E, L, Q into f
        integral += f(z, z_minus, z_plus, a, E, L, Q) * radius * dtheta * (-sin(theta) + I * cos(theta));  // trapezoidal rule
    }
    return integral;
}

// Contour integral calculation main function //
int main(int argc, char **argv) {
    // double a = 0.5;  // spin
    // double E = 0.9;  // energy
    // double L = 4.0;  // angular momentum
    // double Q = 2.0;  // Carter constant

    double a;  // spin
    double E;  // energy
    double L;  // angular momentum
    double Q;  // Carter constant

    sscanf(argv[1], "%lg", &a);
    sscanf(argv[2], "%lg", &E);
    sscanf(argv[3], "%lg", &L);
    sscanf(argv[4], "%lg", &Q);


    double chi = calculate_chi(a, E, L, Q);
    double r = calculate_r(L, Q);
    
    complex double z_minus, z_plus;
    calculate_z_roots(chi, r, &z_minus, &z_plus);
    
    double radius = 2 * ((1 / z_minus) + (1 / z_plus));  // (harmonic mean) radius of the contour
    complex double z_start = radius + 0.0 * I;  // starting point on the contour
    int num_points = 1000;  // number of points on the contour

    // Perform the contour integration
    complex double integral_d2JdE2 = contour_integral(d2JdE2, z_start, radius, num_points, a, E, L, Q, z_minus, z_plus); 
    complex double integral_d2JdEdL = contour_integral(d2JdEdL, z_start, radius, num_points, a, E, L, Q, z_minus, z_plus);
    complex double integral_d2JdL2 = contour_integral(d2JdL2, z_start, radius, num_points, a, E, L, Q, z_minus, z_plus);

    // Print the result
    printf("d2JdE2 contour integral: %.10f + %.10f I\n", creal(integral_d2JdE2), cimag(integral_d2JdE2));
    printf("d2JdEdL contour integral: %.10f + %.10f I\n", creal(integral_d2JdEdL), cimag(integral_d2JdEdL));
    printf("d2JdL2 contour integral: %.10f + %.10f I\n", creal(integral_d2JdL2), cimag(integral_d2JdL2));

    return 0;
}



/* SECTION 2: dM/dJ CALCULATION*/
// dM/dJ calculation (explain that the double derivative is just Minverse w/o having to explicitly calculate)
// inverse of Minverse





/* SECTION 3: dM/dJ -> dOmega/dJ CALCULATION */






/* SECTION 4: PHASE CHANGE CALCULATION */




