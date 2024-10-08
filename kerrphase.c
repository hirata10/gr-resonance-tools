#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
//#include "CKerr.h"
#include "globalpars_c.h"

#define CKERR_NBISECT_ITER 56
#define CKERR_NBISECT_ITER2 32
#define CKERR_RESOLUTION_J_INTEG 64
#define CKERR_J_TOL 1e-9
#define CKERR_ACTION_MAX 1e49
#define CKERR_NPOINT_DERIV 3
#define CKERR_DACTION_DERIV 5e-4

// REMEMBER WRAPPER!!!

/* SECTION 1: Polar CONTOUR INTEGRAL */


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

// Polynomial input to the integrals (equation A13 in the paper.) //
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
complex double d2J_zdE2(complex double z, complex double z_minus, complex double z_plus, double a, double E, double L, double Q) {
    complex double numerator = (P(z, a, E, L, Q) - (a * a * z * z * E * E * (1 - z * z))) * (a * a * z * z);
    complex double denominator = 2 * M_PI * sqrt_P(z, z_minus, z_plus, a, E, L, Q) * P(z, a, E, L, Q);

    return numerator / denominator;
}

// Function to evaluate equation A8 in the paper //
complex double d2J_zdL2(complex double z, complex double z_minus, complex double z_plus, double a, double E, double L, double Q) {
    complex double numerator = - (P(z, a, E, L, Q) + (z * z * L * L)) * (z * z);
    complex double denominator = 2 * M_PI * sqrt_P(z, z_minus, z_plus, a, E, L, Q) * P(z, a, E, L, Q) * (1 - z * z);

    return numerator / denominator;
}

// Function to evaluate equation A9 in the paper //
complex double d2J_zdEdL(complex double z, complex double z_minus, complex double z_plus, double a, double E, double L, double Q) {
    complex double numerator = a * a * z * z * z * z * E * L;
    complex double denominator = 2 * M_PI * sqrt_P(z, z_minus, z_plus, a, E, L, Q) * P(z, a, E, L, Q);

    return numerator / denominator;
}

// Function to evaluate equation A10 in the paper //
complex double d2J_zdQ2(complex double z, complex double z_minus, complex double z_plus, double a, double E, double L, double Q) {
    complex double numerator = - (1 - z*z);
    complex double denominator = 8 * M_PI * sqrt_P(z, z_minus, z_plus, a, E, L, Q) * P(z, a, E, L, Q);

    return numerator / denominator;
}

// Function to evaluate equation A11 in the paper //
complex double d2J_zdEQ(complex double z, complex double z_minus, complex double z_plus, double a, double E, double L, double Q) {
    complex double numerator = - (1 - z*z) * z*z * a*a * E;
    complex double denominator = 4 * M_PI * sqrt_P(z, z_minus, z_plus, a, E, L, Q) * P(z, a, E, L, Q);

    return numerator / denominator;
}

// Function to evaluate equation A12 in the paper //
complex double d2J_zdLQ(complex double z, complex double z_minus, complex double z_plus, double a, double E, double L, double Q) {
    complex double numerator = L * z;
    complex double denominator = 4 * M_PI * sqrt_P(z, z_minus, z_plus, a, E, L, Q) * P(z, a, E, L, Q);

    return numerator / denominator;
}

// Function to perform the contour integration //
complex double contour_integral_z(complex double (*f)(complex double, complex double, complex double, double, double, double, double),
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

/* Radial Contour Integration */

int radial_roots(double *EQL, double *J, double M, double astar, double *ancillary, double *r_roots) {
  long i;
  double zminus2, zplusinv2, dz, value, a;
  double zminus, zplusinv, zmzp;
  double Jt_integ, y, cosy, siny, arg, dy, secharg;
  double K;
  double c[5], disc, rinfl, rinfl2, rpeak, rdip, r, rh;
  double rturn[4];
  double rmid, drmax, Jr, ur;
  double denom;
  double r_h[2];
  int inflsign = 0;
  int returnflag = 1;

  a = M*astar;

  /* Reject unbound orbits */
  if (EQL[0]<=-1) return(0);
  if (EQL[0]>=1) return(2);

  /* Azimuthal action is easy */
  J[2] = EQL[2];
  /* Radial action -- first get roots of r */
  K = EQL[1] + (EQL[2]-a*EQL[0])*(EQL[2]-a*EQL[0]);  
  /* We know that for |E|<1 the polynomial u_r^2/Delta^2 is 4th order, is positive at outer horizon,
   * and negative at infinity.  There is a stable orbit or not depending on whether the polynomial
   * has 3 roots at r>rh (stable) or 1 root (unstable).  The coefficients are:
   */
  c[0] = -a*a*EQL[1];                                             /* negative or 0 */
  c[1] = 2*M*(EQL[1] + (EQL[2]-a*EQL[0])*(EQL[2]-a*EQL[0]));      /* positive or 0 */
  c[2] = -a*a*(1.-EQL[0]*EQL[0]) - EQL[2]*EQL[2] - EQL[1];        /* negative      */
  c[3] = 2.*M;                                                    /* positive      */
  c[4] = EQL[0]*EQL[0]-1.;                                        /* negative      */

  /* Find inflection point */
  disc = 36.*c[3]*c[3] - 96.*c[4]*c[2];
  if (disc<=0) {disc=0.; returnflag=0;} /* no inflection point */
  rinfl = (-6.*c[3]-sqrt(disc))/(24.*c[4]);
  rinfl2 = (-6.*c[3]+sqrt(disc))/(24.*c[4]);
  if (rinfl <= M + sqrt(M*M-a*a)) inflsign = 1;

  /* Find the greatest local maximum & local minimum of the polynomial if it is > rinfl */
  rpeak = rinfl;
  if (((4.*c[4]*rpeak + 3.*c[3])*rpeak + 2.*c[2])*rpeak + c[1] < 0) {
    inflsign = 1;
    returnflag=0; /* no peak after 2nd inflection point; no stable orbits */
  }
  while (((4.*c[4]*rpeak + 3.*c[3])*rpeak + 2.*c[2])*rpeak + c[1] >0) {
    rpeak *= 2;
    if (rpeak>1e20*M) {
#if 1
      fprintf(stderr, "Warning: peak out of range, breaking\n");
#endif
      break;
    }
  }
  dz = 0.5;  
  for(i=0; i<CKERR_NBISECT_ITER;i++) {
    rpeak *= ((4.*c[4]*rpeak + 3.*c[3])*rpeak + 2.*c[2])*rpeak + c[1] > 0? pow(2,dz): pow(2,-dz);    
    dz/=2.;
  }
  rdip = 0.5*(rinfl+rinfl2);
  dz = 0.25*(rinfl-rinfl2);
  for(i=0; i<CKERR_NBISECT_ITER;i++) {
    rdip += ((4.*c[4]*rdip + 3.*c[3])*rdip + 2.*c[2])*rdip + c[1] < 0? dz: -dz;
    dz/=2.;
  }

  /* Are these peaks successively below and above zero? */
  r = rpeak;
  if ( (((c[4]*r+c[3])*r+c[2])*r+c[1])*r+c[0] < 0 ) returnflag=0;
  r = rdip;
  if ( (((c[4]*r+c[3])*r+c[2])*r+c[1])*r+c[0] > 0 ) {
    returnflag=0;
    if (inflsign==0) returnflag = 2;
  }
#if 0
fprintf(stderr, "%12.5le,%12.5le,%12.5le\n", r, (((c[4]*r+c[3])*r+c[2])*r+c[1])*r+c[0], ((4.*c[4]*r + 3.*c[3])*r + 2.*c[2])*r + c[1]);
#endif

  /* Find turning points */
  r = 2*rpeak;
  while ( (((c[4]*r+c[3])*r+c[2])*r+c[1])*r+c[0] >= 0 ) {
    r*=2;
  }
  dz = 0.5;
  for(i=0; i<CKERR_NBISECT_ITER;i++) {
    r *= (((c[4]*r+c[3])*r+c[2])*r+c[1])*r+c[0] >= 0? pow(2,dz): pow(2,-dz);
    dz/=2.;
  }
  rturn[3] = r;
  r = 0.5*(rpeak+rdip);
  dz = 0.25*(rpeak-rdip);
  for(i=0; i<CKERR_NBISECT_ITER;i++) {
    r += (((c[4]*r+c[3])*r+c[2])*r+c[1])*r+c[0] >= 0? -dz: dz;
    dz/=2.;
  }
  rturn[2] = r;
  rh = M+sqrt(M*M-a*a);
  r = 0.5*(rdip+rh);
  dz = 0.25*(rdip-rh);
  for(i=0; i<CKERR_NBISECT_ITER;i++) {
    r += (((c[4]*r+c[3])*r+c[2])*r+c[1])*r+c[0] >= 0? dz: -dz;
    dz/=2.;
  }
  rturn[1] = r;
  r = 0.5*rh;
  dz = 0.25*rh;
  for(i=0; i<CKERR_NBISECT_ITER;i++) {
    r += (((c[4]*r+c[3])*r+c[2])*r+c[1])*r+c[0] >= 0? -dz: dz;
    dz/=2.;
  }
  rturn[0] = r;

  /* Defining inner [0] and outer [1] horizon */
  r_h[0] = M - sqrt(M*M - a*a);
  r_h[1] = M + sqrt(M*M - a*a);

  /* Output into r_roots array */
  r_roots[0] = rturn[0];
  r_roots[1] = rturn[1];
  r_roots[2] = rturn[2];
  r_roots[3] = rturn[3];

  //fprintf(stderr, "r_infl = %12.5le,%12.5le r_dip,peak = %12.5le,%12.5le inner and outer horizon = %12.5le %12.5le roots = %12.5le %12.5le %12.5le %12.5le\n", rinfl2, rinfl, rdip, rpeak, r_h[0], r_h[1], rturn[0], rturn[1], rturn[2], rturn[3]);
    
}

// Polynomial input to the radial integrals (equation A20 in the paper.) //
complex double P_r(complex double z, double a, double M, double E, double Q, double L) {
    return ((z*z + a*a) * E - a * L) * ((z*z + a*a) * E - a * L) - (z*z -2*M*z + a*a) * (z*z + Q + (L - a*E) * (L - a*E));
}

// Square root of polynomial function above //
complex double sqrt_P_r(complex double z, complex double r_root1, complex double r_root2, complex double r_root3, complex double r_root4, double a, double E) {
    return I * sqrt(1 - E * E) * z * csqrt(1 - r_root1/z) * csqrt(1 - r_root2/z) * csqrt(z - r_root3) * csqrt(z - r_root4);
}

// Function to evaluate equation A14 in the paper //
complex double d2J_rdE2(complex double z, complex double r_root1, complex double r_root2, complex double r_root3, complex double r_root4, double a, double M, double E, double Q, double L) {
    complex double numerator = ((z*z + a*a) * (z*z + a*a) - a*a * (z*z - 2*M*z +a*a)) * P_r(z, a, M, E, Q, L) - (((z*z + a*a) * E - a * L) * (z*z + a*a) + a * (z*z - 2*M*z + a*a) * (L - a * E)) * (((z*z + a*a) * E - a * L) * (z*z + a*a) + a * (z*z - 2*M*z + a*a) * (L - a * E));
    complex double denominator = 2 * M_PI * (z*z - 2*M*z + a*a) * sqrt_P_r(z, r_root1, r_root2, r_root3, r_root4, a, E) * P_r(z, a, M, E, L, Q);

    return numerator / denominator;
}

// Function to evaluate equation A15 in the paper //
complex double d2J_rdL2(complex double z, complex double r_root1, complex double r_root2, complex double r_root3, complex double r_root4, double a, double M, double E, double Q, double L) {
    complex double numerator = (a*a - (z*z - 2*M*z + a*a)) * P_r(z, a, M, E, Q, L) - (((z*z + a*a) * E - a * L) * a + (z*z - 2*M*z + a*a) * (L - a * E)) * (((z*z + a*a) * E - a * L) * a + (z*z - 2*M*z + a*a) * (L - a * E));
    complex double denominator = 2 * M_PI * (z*z - 2*M*z + a*a) * sqrt_P_r(z, r_root1, r_root2, r_root3, r_root4, a, E) * P_r(z, a, M, E, L, Q);

    return numerator / denominator;
}

// Function to evaluate equation A16 in the paper //
complex double d2J_rdQ2(complex double z, complex double r_root1, complex double r_root2, complex double r_root3, complex double r_root4, double a, double M, double E, double Q, double L) {
    complex double numerator = - (z*z - 2*M*z + a*a);
    complex double denominator = 8 * M_PI * sqrt_P_r(z, r_root1, r_root2, r_root3, r_root4, a, E) * P_r(z, a, M, E, L, Q);

    return numerator / denominator;
}

// Function to evaluate equation A17 in the paper //
complex double d2J_rdEL(complex double z, complex double r_root1, complex double r_root2, complex double r_root3, complex double r_root4, double a, double M, double E, double Q, double L) {
    complex double numerator = (a * (z*z - 2*M*z + a*a) - a * (z*z + a*a)) * P_r(z, a, M, E, Q, L) + (((z*z + a*a) * E - a * L) * (z*z + a*a) + a * (z*z - 2*M*z + a*a) * (L - a * E)) * (((z*z + a*a) * E - a * L) * a + (z*z - 2*M*z + a*a) * (L - a * E));
    complex double denominator = 2 * M_PI * (z*z - 2*M*z + a*a) * sqrt_P_r(z, r_root1, r_root2, r_root3, r_root4, a, E) * P_r(z, a, M, E, L, Q);

    return numerator / denominator;
}

// Function to evaluate equation A18 in the paper //
complex double d2J_rdEQ(complex double z, complex double r_root1, complex double r_root2, complex double r_root3, complex double r_root4, double a, double M, double E, double Q, double L) {
    complex double numerator = (((z*z + a*a) * E - a * L) * (z*z + a*a) + a * (z*z - 2*M*z + a*a) * (L - a * E));
    complex double denominator = 4 * M_PI * sqrt_P_r(z, r_root1, r_root2, r_root3, r_root4, a, E) * P_r(z, a, M, E, L, Q);

    return numerator / denominator;
}

// Function to evaluate equation A19 in the paper //
complex double d2J_rdLQ(complex double z, complex double r_root1, complex double r_root2, complex double r_root3, complex double r_root4, double a, double M, double E, double Q, double L) {
    complex double numerator =  - (((z*z + a*a) * E - a * L) * a + (z*z - 2*M*z + a*a) * (L - a * E));
    complex double denominator = 4 * M_PI * sqrt_P_r(z, r_root1, r_root2, r_root3, r_root4, a, E) * P_r(z, a, M, E, L, Q);

    return numerator / denominator;
}

// Function to perform the contour integration //
complex double contour_integral_r(complex double (*f)(complex double, complex double, complex double, complex double, complex double, double, double, double, double, double),
                                            complex double r_start, double radius, int num_points, double center,
                                            double a, double M, double E, double L, double Q,
                                            complex double r_root1, complex double r_root2, complex double r_root3, complex double r_root4) {
    complex double integral = 0.0 + 0.0 * I;
    double theta;
    complex double z;
    double dtheta = 2 * M_PI / num_points;  // angle step

    for (int i = 0; i < num_points; i++) {
        theta = i * dtheta;
        z =  center + radius * (cos(theta) + I * sin(theta));  // point on the contour
        // Now passing a, E, L, Q into f
        integral += f(z, r_root1, r_root2, r_root3, r_root4, a, M, E, Q, L) * radius * dtheta * (-sin(theta) + I * cos(theta));  // trapezoidal rule
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
    double M; // mass of black hole
    double EQL[3], ancillary[3], r_roots[4]; // array for EQL and ancillary data
    double J[3]; // array for action variables

    sscanf(argv[1], "%lg", &a);
    sscanf(argv[2], "%lg", &E);
    sscanf(argv[3], "%lg", &L);
    sscanf(argv[4], "%lg", &Q);
    sscanf(argv[5], "%lg", &M);

    EQL[0] = E;
    EQL[1] = Q;
    EQL[2] = L;

    #if 0
    /* For polar integral */
    double chi = calculate_chi(a, E, L, Q);
    double r = calculate_r(L, Q);
    
    complex double z_minus, z_plus;
    calculate_z_roots(chi, r, &z_minus, &z_plus);
    
    double radius = 2 * ((1 / z_minus) + (1 / z_plus));  // (harmonic mean) radius of the contour
    complex double z_start = radius + 0.0 * I;  // starting point on the contour
    int num_points = 1000;  // number of points on the contour

    // Perform the contour integration
    complex double integral_d2JdE2 = contour_integral_z(d2J_zdE2, z_start, radius, num_points, a, E, L, Q, z_minus, z_plus); 
    complex double integral_d2JdEdL = contour_integral_z(d2J_zdEdL, z_start, radius, num_points, a, E, L, Q, z_minus, z_plus);
    complex double integral_d2JdL2 = contour_integral_z(d2J_zdL2, z_start, radius, num_points, a, E, L, Q, z_minus, z_plus);
    complex double integral_d2JdEQ = contour_integral_z(d2J_zdEQ, z_start, radius, num_points, a, E, L, Q, z_minus, z_plus);
    complex double integral_d2JdQ2 = contour_integral_z(d2J_zdQ2, z_start, radius, num_points, a, E, L, Q, z_minus, z_plus);
    complex double integral_d2JdLQ = contour_integral_z(d2J_zdLQ, z_start, radius, num_points, a, E, L, Q, z_minus, z_plus);
     

    // Print the result
    printf("d2JdE2 contour integral: %.10f + %.10f I\n", creal(integral_d2JdE2), cimag(integral_d2JdE2));
    printf("d2JdEdL contour integral: %.10f + %.10f I\n", creal(integral_d2JdEdL), cimag(integral_d2JdEdL));
    printf("d2JdL2 contour integral: %.10f + %.10f I\n", creal(integral_d2JdL2), cimag(integral_d2JdL2));
    printf("d2JdQ2 contour integral: %.10f + %.10f I\n", creal(integral_d2JdQ2), cimag(integral_d2JdQ2));
    printf("d2JdEQ contour integral: %.10f + %.10f I\n", creal(integral_d2JdEQ), cimag(integral_d2JdEQ));
    printf("d2JdLQ contour integral: %.10f + %.10f I\n", creal(integral_d2JdLQ), cimag(integral_d2JdLQ));
    #endif

    /* For radial integral */
    radial_roots(EQL, J, M, a, ancillary, r_roots); //compute radial roots

    double chi = calculate_chi(a, E, L, Q);
    double r = calculate_r(L, Q);
    
    complex double z_minus, z_plus;
    calculate_z_roots(chi, r, &z_minus, &z_plus);
    
    double radius = 2 * ((1 / z_minus) + (1 / z_plus));  // (harmonic mean) radius of the contour
    complex double z_start = radius + 0.0 * I;  // starting point on the contour

    int num_points = 1000;  // number of points on the contour
    double radius_r = (r_roots[2] - r_roots[1]); //radius of radial contour
    double center = (r_roots[2] + r_roots[3])/2.0; //center of circular contour
    complex double r_start = radius_r + 0.0 * I; //starting point on contour

    double r_H_outer = M + sqrt(M*M - a*a);
    double r_H_inner = M - sqrt(M*M - a*a);

    printf("sqrtP(r) real and imaginary %0.1f + %0.1f I \n", creal(sqrt_P_r(r_start, r_roots[0], r_roots[1], r_roots[2], r_roots[3], a, E)), cimag(sqrt_P_r(r_start, r_roots[0], r_roots[1], r_roots[2], r_roots[3], a, E)));

    printf("sqrtP(z) real and imaginary %0.1f + %0.1f I \n", creal(sqrt_P(z_start, z_minus, z_plus, a, E, L, Q)), cimag(sqrt_P_r(z_start, z_minus, z_plus, a, E, L, Q)));

    printf("inner and outer horizon: %.10f and %.10f \n", r_H_inner, r_H_outer);
    complex double integral_d2J_r_dE2 = contour_integral_r(d2J_rdE2, r_start, radius_r, num_points, center, a, M, E, Q, L, r_roots[0], r_roots[1], r_roots[2], r_roots[3]);

    printf("d2J_rdE2 contour integral: %.10f + %.10f I\n", creal(integral_d2J_r_dE2), cimag(integral_d2J_r_dE2));
    return 0;
}



/* SECTION 2: dM/dJ CALCULATION*/
// dM/dJ calculation (explain that the double derivative is just Minverse w/o having to explicitly calculate)
// inverse of Minverse





/* SECTION 3: dM/dJ -> dOmega/dJ CALCULATION */






/* SECTION 4: PHASE CHANGE CALCULATION */




