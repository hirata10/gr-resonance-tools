#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#define CKERR_NBISECT_ITER 56
#define CKERR_NBISECT_ITER2 32
#define CKERR_RESOLUTION_J_INTEG 64
#define CKERR_J_TOL 1e-9
#define CKERR_ACTION_MAX 1e49
#define CKERR_NPOINT_DERIV 3
#define CKERR_DACTION_DERIV 5e-4
#undef I // Undeclare I since we use it for inclination in CKerr.h
#define complex_I _Complex_I //define complex number I and complex_I

#define SIZE 3  // Size of the matrix (3x3)
#define TOTAL_ELEMENTS (SIZE * SIZE)  // Total elements in the matrix


#include "CKerr.h"
#include "globalpars_c.h"

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
    return complex_I * a * sqrt(1 - E * E) * z * csqrt(1 - ((z_minus * z_minus)/(z * z))) * csqrt((z_plus * z_plus) - (z * z));
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
    complex double integral = 0.0 + 0.0 * complex_I;
    double theta;
    complex double z;
    double dtheta = 2 * M_PI / num_points;  // angle step

    for (int i = 0; i < num_points; i++) {
        theta = i * dtheta;
        z = radius * (cos(theta) + complex_I * sin(theta));  // point on the contour
        // Now passing a, E, L, Q into f
        integral += f(z, z_minus, z_plus, a, E, L, Q) * radius * dtheta * (-sin(theta) + complex_I * cos(theta));  // trapezoidal rule
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
    return complex_I * sqrt(1 - E * E) * z * csqrt(1 - r_root1/z) * csqrt(1 - r_root2/z) * csqrt(z - r_root3) * csqrt(z - r_root4);
    //return complex_I * sqrt(1 - E * E) * csqrt(z - r_root1) * csqrt(z - r_root2) * csqrt(z - r_root3) * csqrt(z - r_root4);
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
                                            double a, double M, double E, double Q, double L,
                                            complex double r_root1, complex double r_root2, complex double r_root3, complex double r_root4) {
    complex double integral = 0.0 + 0.0 * complex_I;
    double theta;
    complex double z;
    double dtheta = 2 * M_PI / num_points;  // angle step

    for (int i = 0; i < num_points; i++) {
        theta = i * dtheta;
        z =  center + radius * (cos(theta) + complex_I * sin(theta));  // point on the contour
        // Now passing a, E, L, Q into f
        integral += f(z, r_root1, r_root2, r_root3, r_root4, a, M, E, Q, L) * radius * dtheta * (-sin(theta) + complex_I * cos(theta));  // trapezoidal rule
    }
    return integral;
}



/* SECTION 2: dM/dJ CALCULATION*/
// dM/dJ calculation (explain that the double derivative is just Minverse w/o having to explicitly calculate)

// inverse of Minverse

int invertMatrix(double *Minverse, double *M_matrix) {
    double determinant = Minverse[0] * (Minverse[4] * Minverse[8] - Minverse[5] * Minverse[7]) -
                         Minverse[1] * (Minverse[3] * Minverse[8] - Minverse[5] * Minverse[6]) +
                         Minverse[2] * (Minverse[3] * Minverse[7] - Minverse[4] * Minverse[6]);

    if (determinant == 0) {
        printf("Matrix is singular and cannot be inverted.\n");
        return 0; // Indicate failure
    }

    double invDet = 1.0 / determinant;

    M_matrix[0] =  (Minverse[4] * Minverse[8] - Minverse[5] * Minverse[7]) * invDet;
    M_matrix[1] = -(Minverse[1] * Minverse[8] - Minverse[2] * Minverse[7]) * invDet;
    M_matrix[2] =  (Minverse[1] * Minverse[5] - Minverse[2] * Minverse[4]) * invDet;

    M_matrix[3] = -(Minverse[3] * Minverse[8] - Minverse[5] * Minverse[6]) * invDet;
    M_matrix[4] =  (Minverse[0] * Minverse[8] - Minverse[2] * Minverse[6]) * invDet;
    M_matrix[5] = -(Minverse[0] * Minverse[5] - Minverse[2] * Minverse[3]) * invDet;

    M_matrix[6] =  (Minverse[3] * Minverse[7] - Minverse[4] * Minverse[6]) * invDet;
    M_matrix[7] = -(Minverse[0] * Minverse[7] - Minverse[1] * Minverse[6]) * invDet;
    M_matrix[8] =  (Minverse[0] * Minverse[4] - Minverse[1] * Minverse[3]) * invDet;

    return 1; // Indicate success
}

/* Input array is the array on contour integrals and then reorganized into the appropriate array for phase change calculations */
int partial_J_derivatives(double *input_array, double *partial_derivatives){
  /* input_array structure:  */
  /* 
  integral_d2J_r_dE2 - 0
  integral_d2J_r_dL2 - 1
  integral_d2J_r_dQ2 - 2
  integral_d2J_r_dEL - 3
  integral_d2J_r_dEQ - 4
  integral_d2J_r_dLQ - 5
  integral_d2JdE2 - 6
  integral_d2JdL2 - 7
  integral_d2JdQ2 - 8
  integral_d2JdEdL - 9
  integral_d2JdEQ - 10
  integral_d2JdLQ - 11
  */
  partial_derivatives[0] = input_array[0]; // jr_ee
  partial_derivatives[1] = input_array[4]; // jr_eq
  partial_derivatives[2] = input_array[3]; // jr_el
  partial_derivatives[3] = input_array[4]; // jr_qe
  partial_derivatives[4] = input_array[2]; // jr_qq
  partial_derivatives[5] = input_array[5]; // jr_ql
  partial_derivatives[6] = input_array[3]; //jr_le
  partial_derivatives[7] = input_array[5]; // jr_lq
  partial_derivatives[8] = input_array[1]; //jr_ll
  partial_derivatives[9] = input_array[6]; //jtheta_ee
  partial_derivatives[10] = input_array[10]; // jtheta_eq
  partial_derivatives[11] = input_array[9]; // jtheta_el
  partial_derivatives[12] = input_array[10]; // jtheta_qe
  partial_derivatives[13] = input_array[8]; // jtheta_qq
  partial_derivatives[14] = input_array[11]; // jtheta_ql
  partial_derivatives[15] = input_array[9]; // jtheta_le
  partial_derivatives[16] = input_array[11]; // jtheta_lq
  partial_derivatives[17] = input_array[7]; // jtheta_ll
  partial_derivatives[18] = 0; // jphi_{anything}, since u_phi is a constant, then all derivatives of J_phi is zero
  partial_derivatives[19] = 0;
  partial_derivatives[20] = 0;
  partial_derivatives[21] = 0;
  partial_derivatives[22] = 0;
  partial_derivatives[23] = 0;
  partial_derivatives[24] = 0;
  partial_derivatives[25] = 0;
  partial_derivatives[26] = 0;
  
  
  return 0;
}

/* SECTION 3: dM/dJ -> dOmega/dJ CALCULATION */

int dM_dJ(double M[TOTAL_ELEMENTS], double partial_derivatives[3 * TOTAL_ELEMENTS], double result[3 * TOTAL_ELEMENTS]) {

    for (int A = 0; A < 3; A++) {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 3; k++) {
                // Initialize the result for this (A, i, k) to zero
                result[A * 3 + i + 3 * 3 * k] = 0;

                for (int j = 0; j < 3; j++) { // Loop over j
                    for (int B = 0; B < 3; B++) {
                        for (int D = 0; D < 3; D++) {
                            // Access M elements based on flattened index
                            double M_Aj = M[A * 3 + j]; // M_Aj for specific A and j
                            double M_Di = M[D * 3 + i]; // M_Di for specific D and i
                            double M_Bk = M[B * 3 + k]; // M_Bk for specific B and k

                            // Update result tensor
                            result[A * 3 + i + 3 * 3 * k] += - M_Aj * M_Di * M_Bk * partial_derivatives[j * 3 * 3 + B * 3 + D];
                    }
                }
            }
        }
    }
}


    return 0;
}





/* SECTION 4: PHASE CHANGE CALCULATION */

// Compute \Delta Phi^i = sum_i J_i/\dot{J}_i,sf * dOmega^i\dJ_k \Delta J_k,td

// Input: EQL, Delta_J_td_i, dOmega_dJ
// Output: array with Delta_Phi for each spatial direction

int Delta_Phi(double E, double Q, double L, double Delta_J_td_r, double Delta_J_td_theta, double Delta_J_td_phi , double *dOmega_dJ, double *Delta_Phi_components, double astar, double M){
  double J[3], ancillary[3], J_dot_sf[3], EQL[3], Delta_J_td[3];
  double t_scale;

  EQL[0] = E;
  EQL[1] = Q;
  EQL[2] = L;

  Delta_J_td[0] = Delta_J_td_r;
  Delta_J_td[1] = Delta_J_td_theta;
  Delta_J_td[2] = Delta_J_td_phi;



  CKerr_EQL2J(EQL, J, M, astar, ancillary); // Convert EQL to J

  J_dot_selfforce(GLOBALPAR_nl_self, GLOBALPAR_nmax, GLOBALPAR_kmax, GLOBALPAR_mmax, ancillary[2], ancillary[1], 0, ancillary[0], M, astar, J_dot_sf); // Compute J_dot_selfforce for time scale

  for (int i = 0; i < 3; i++){
    printf("%d: %lg \n", i, J[i] / J_dot_sf[i]);
    t_scale += J[i] * J[i] / (J_dot_sf[i] * J_dot_sf[i]); // J_i/J_dot_i added in quadrature
    t_scale = sqrt(t_scale); // Time scale for resonance crossing
  }

  for (int i = 0; i < 3; i++){
    for (int k = 0; k < 3; k++){
      Delta_Phi_components[i] += t_scale * dOmega_dJ[i + 9 * k] * Delta_J_td[k]; // Computing the change in Phi for each component
    }
  }


return (0);
}



// Contour integral calculation main function //
int main(int argc, char **argv) {
    // double a = 0.5;  // spin
    // double E = 0.9;  // energy
    // double L = 4.0;  // angular momentum
    // double Q = 2.0;  // Carter constant
    int i, j;
    int n_inner, k_inner, m_inner;
    double a;  // spin
    double E;  // energy
    double L;  // angular momentum
    double Q;  // Carter constant
    double M; // mass of black hole
    double EQL[3], ancillary[3], r_roots[4]; // array for EQL and ancillary data
    double J[3], Minv[9], M_matrix[9], J_r, J_theta, J_phi; // array for action variables and Minv matrix
    double input_array[12], partial_derivatives[27], dM_dJ_result[27]; //Arrays for storing partial derivatives of J
    double Delta_J_td_r, Delta_J_td_theta, Delta_J_td_phi;
    double Delta_Phi_components[3], tot_Delta_Phi;

    sscanf(argv[1], "%lg", &a);
    // sscanf(argv[2], "%lg", &E);
    // sscanf(argv[3], "%lg", &Q);
    // sscanf(argv[4], "%lg", &L);
    sscanf(argv[2], "%lg", &J_r);
    sscanf(argv[3], "%lg", &J_theta);
    sscanf(argv[4], "%lg", &J_phi);
    sscanf(argv[5], "%lg", &M);
    sscanf(argv[6], "%lg", &Delta_J_td_r);
    sscanf(argv[7], "%lg", &Delta_J_td_theta);
    sscanf(argv[8], "%lg", &Delta_J_td_phi);
    sscanf(argv[9], "%i", &n_inner);
    sscanf(argv[10], "%i", &k_inner);
    sscanf(argv[11], "%i", &m_inner);

    
    J[0] = J_r;
    J[1] = J_theta;
    J[2] = J_phi;

    //CKerr_EQL2J(EQL, J, M, a, ancillary);
    CKerr_J2EQL(J, EQL, M, a);
    CKerr_Minverse(J, Minv, M, a);
    invertMatrix(Minv, M_matrix);

    E = EQL[0];
    Q = EQL[1];
    L = EQL[2];

    //printf("%lg %lg %lg \n", E, Q, L);
    
    printf("Minv components: \n");
    for (i = 0; i < 9; i++){
      printf("%lg", Minv[i]);
      printf("\n");
    }

    printf("M_matrix components: \n");
    for (i = 0; i < 9; i++){
      printf("%lg", M_matrix[i]);
      printf("\n");
    }
    /* For polar integral */
    double chi = calculate_chi(a, E, L, Q);
    double r = calculate_r(L, Q);
    
    complex double z_minus, z_plus;
    calculate_z_roots(chi, r, &z_minus, &z_plus);
    
    double radius = 2 * ((1 / z_minus) + (1 / z_plus));  // (harmonic mean) radius of the contour
    complex double z_start = radius + 0.0 * complex_I;  // starting point on the contour
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
    

    /* For radial integral */
    radial_roots(EQL, J, M, a, ancillary, r_roots); //compute radial roots

    double radius_r = (r_roots[3] - r_roots[1])/2.0; //radius of circular contour
    double center = (r_roots[2] + r_roots[3])/2.0; //center of circular contour
    complex double r_start = radius_r + 0.0 * complex_I; //starting point on contour

    double r_H_outer = M + sqrt(M*M - a*a);
    double r_H_inner = M - sqrt(M*M - a*a);

    // printf("radial roots: %.10f %.10f %.10f %.10f inner and outer horizon: %.10f %.10f \n", r_roots[0], r_roots[1], r_roots[2], r_roots[3], r_H_inner, r_H_outer);

    // printf("sqrtP(r) real and imaginary %.10f + %.10f I \n", creal(sqrt_P_r(r_start, r_roots[0], r_roots[1], r_roots[2], r_roots[3], a, E)), cimag(sqrt_P_r(r_start, r_roots[0], r_roots[1], r_roots[2], r_roots[3], a, E)));

    // printf("sqrtP(z) real and imaginary %.10f + %.10f I \n", creal(sqrt_P(z_start, z_minus, z_plus, a, E, L, Q)), cimag(sqrt_P_r(z_start, z_minus, z_plus, a, E, L, Q)));

    //printf("inner and outer horizon: %.10f and %.10f \n", r_H_inner, r_H_outer);
   
   /* Compute the radial contour integrals */
    complex double integral_d2J_r_dE2 = contour_integral_r(d2J_rdE2, r_start, radius_r, num_points, center, a, M, E, Q, L, r_roots[0], r_roots[1], r_roots[2], r_roots[3]);
    complex double integral_d2J_r_dL2 = contour_integral_r(d2J_rdL2, r_start, radius_r, num_points, center, a, M, E, Q, L, r_roots[0], r_roots[1], r_roots[2], r_roots[3]);
    complex double integral_d2J_r_dQ2 = contour_integral_r(d2J_rdQ2, r_start, radius_r, num_points, center, a, M, E, Q, L, r_roots[0], r_roots[1], r_roots[2], r_roots[3]);
    complex double integral_d2J_r_dEL = contour_integral_r(d2J_rdEL, r_start, radius_r, num_points, center, a, M, E, Q, L, r_roots[0], r_roots[1], r_roots[2], r_roots[3]);
    complex double integral_d2J_r_dEQ = contour_integral_r(d2J_rdEQ, r_start, radius_r, num_points, center, a, M, E, Q, L, r_roots[0], r_roots[1], r_roots[2], r_roots[3]);
    complex double integral_d2J_r_dLQ = contour_integral_r(d2J_rdLQ, r_start, radius_r, num_points, center, a, M, E, Q, L, r_roots[0], r_roots[1], r_roots[2], r_roots[3]);
   
    printf("d2J_rdE2 contour integral: %.10f + %.10f I\n", creal(integral_d2J_r_dE2), cimag(integral_d2J_r_dE2));
    printf("d2J_rdL2 contour integral: %.10f + %.10f I\n", creal(integral_d2J_r_dL2), cimag(integral_d2J_r_dL2));
    printf("d2J_rdQ2 contour integral: %.10f + %.10f I\n", creal(integral_d2J_r_dQ2), cimag(integral_d2J_r_dQ2));
    printf("d2J_rdEL contour integral: %.10f + %.10f I\n", creal(integral_d2J_r_dEL), cimag(integral_d2J_r_dEL));
    printf("d2J_rdEQ contour integral: %.10f + %.10f I\n", creal(integral_d2J_r_dEQ), cimag(integral_d2J_r_dEQ));
    printf("d2J_rdLQ contour integral: %.10f + %.10f I\n", creal(integral_d2J_r_dLQ), cimag(integral_d2J_r_dLQ));

    /* Define input array for partial derivatives of J array */
    input_array[0] = integral_d2J_r_dE2;
    input_array[1] = integral_d2J_r_dL2;
    input_array[2] = integral_d2J_r_dQ2;
    input_array[3] = integral_d2J_r_dEL;
    input_array[4] = integral_d2J_r_dEQ;
    input_array[5] = integral_d2J_r_dLQ;
    input_array[6] = integral_d2JdE2;
    input_array[7] = integral_d2JdL2;
    input_array[8] = integral_d2JdQ2;
    input_array[9] = integral_d2JdEdL;
    input_array[10] = integral_d2JdEQ;
    input_array[11] = integral_d2JdLQ;

    partial_J_derivatives(input_array, partial_derivatives); // Compute the array of partial derivatives of J

    printf("dJ_dBD components: \n");
    for (i = 0; i < 27; i++){
      printf("%lg", partial_derivatives[i]);
      printf("\n");
    }

    dM_dJ(M_matrix, partial_derivatives, dM_dJ_result); // Compute dM_dJ

    printf("dM_dJ components: \n");
    for (i = 0; i < 27; i++){
      printf("%lg", dM_dJ_result[i]);
      printf("\n");
    }

    printf("dOmega_dJ components: \n"); // Computes dOmega_dJ = dM_dJ[i + 0 * A + 9 * k]
    for (i = 0; i < 3; i++){
      for (j = 0; j < 3; j++){
        printf("%lg", dM_dJ_result[i + 9 * j]);
        printf("\n");
      }
    }

    Delta_Phi(E, Q, L, Delta_J_td_r, Delta_J_td_theta, Delta_J_td_phi, dM_dJ_result, Delta_Phi_components, a, M); // Compute Delta Phi for each direction

    printf("Change in phase components: \n");
    for (i = 0; i < 3; i++){
      printf("%lg \n", Delta_Phi_components[i]);
    }
    
    tot_Delta_Phi = n_inner * Delta_Phi_components[0] + k_inner * Delta_Phi_components[1] + m_inner * Delta_Phi_components[2]; // Total change in Phi over a resonance crossing with mode (n,k,m)_inner
    printf("Total change in phase over resonance mode (%i, %i, %i): %lg \n", n_inner, k_inner, m_inner, tot_Delta_Phi);
    return 0;
}

