#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#undef I // Undeclare I since we use it for inclination in CKerr.h
#define complex_I _Complex_I //define complex number I and complex_I

#include "globalpars_c.h"
#include "CKerr.h"


// Function to solve Kepler's equation E - e * sin(E) = M_kepler using Newton's method
double kepler_eqn(double M_kepler, double e) {
    double E = M_kepler;  // Initial guess for E
    double delta, f_E, f_prime_E;
    int max_iter = 100;
    double tolerance = 1e-10;

    for (int i = 0; i < max_iter; i++) {
        f_E = E - e * sin(E) - M_kepler;
        f_prime_E = 1 - e * cos(E);
        delta = f_E / f_prime_E;
        E -= delta;

        if (fabs(delta) < tolerance) {
            break;
        }
    }
    return E;
}

// Function to calculate the true anomaly from the eccentric anomaly
double true_anomaly(double E, double e) {
    if (e >= 1.0) {
        printf("Warning: Eccentricity too large for true anomaly calculation.\n");
        return NAN;
    }
    // return 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
    return atan2(sqrt(1 - e*e) * sin(E), cos(E - e)); // A more stable atan except for when e is close to 1
}

// Complex integrand f(nu, M_kepler) = exp(i M_kepler (n + k + m)) * exp(i (m + k) nu) * (1 + e * cos(nu))^(1 + ell)
// Integrating over mean anomaly (M_kepler) since we can express nu (true anomaly) as a function of M_kepler
double complex integrand_Z_down(double M_kepler, double e, double semi_major, int n, int k, int m, int ell, double mass_body, double omega_gw) {
    // Solve Kepler's equation to get eccentric anomaly E
    double E = kepler_eqn(M_kepler, e);
    
    // Compute true anomaly nu from eccentric anomaly E (although M = nu for circular orbits)
    // In practice, we can use M and E directly for elliptical orbits.

    if (isnan(E)) {
        printf("NaN detected in Kepler's equation solution.\n");
        return 0.0 + 0.0 * complex_I;
    }

    double nu = true_anomaly(E, e);
    if (isnan(nu)) {
        printf("NaN detected in true anomaly calculation.\n");
        return 0.0 + 0.0 * complex_I;
    }

    // Check if the arguments of complex exponentials are reasonable
    // printf("M: %lg, nu: %lg, E: %lg\n", M, nu, E);  // Debugging statement
    
    // Complex exponentials
    double complex term1 = cexp(complex_I * M_kepler * (n + k + m)); // exp(i M (n + k + m))
    double complex term2 = cexp(complex_I * (m + k) * nu);   // exp(i (m + k) nu)
    
    // Real term (1 + e * cos(nu))^(-1 + ell)
    double real_term = pow(1 + e * cos(nu), (1 + ell));

    double prefactor_without_k1angularparts = - mass_body / (2 * (2 * ell + 1) * 2 * M_PI) * sqrt(2 * M_PI * tgamma(ell + 3) / tgamma(ell - 1)) *  pow(1/(semi_major * (1-e*e)), (1 + ell));
    
    // Return the product of the complex exponentials and the real term
    // printf("%lg %lg \n", creal(prefactor_without_angularparts * term1 * term2 * real_term), cimag(prefactor_without_angularparts * term1 * term2 * real_term));
    return prefactor_without_k1angularparts * term1 * term2 * real_term;

    // If prefactor or other terms result in a NaN, we should print them
    if (isnan(prefactor_without_k1angularparts) || isnan(nu) || isnan(E)) {
        printf("NaN detected: prefactor=%lg, nu=%lg, E=%lg\n", prefactor_without_k1angularparts, nu, E);
        return 0.0 + 0.0 * complex_I;  // Return zero if NaN detected
    }
}

double complex integrand_Z_out(double M_kepler, double e, double semi_major, int n, int k, int m, int ell, double mass_body, double omega_gw) {
    // Solve Kepler's equation to get eccentric anomaly E
    double E = kepler_eqn(M_kepler, e);
    
    // Compute true anomaly nu from eccentric anomaly E (although M = nu for circular orbits)
    // In practice, we can use M and E directly for elliptical orbits.

    if (isnan(E)) {
        printf("NaN detected in Kepler's equation solution.\n");
        return 0.0 + 0.0 * complex_I;
    }

    double nu = true_anomaly(E, e);
    if (isnan(nu)) {
        printf("NaN detected in true anomaly calculation.\n");
        return 0.0 + 0.0 * complex_I;
    }

    // Check if the arguments of complex exponentials are reasonable
    // printf("M: %lg, nu: %lg, E: %lg\n", M, nu, E);  // Debugging statement
    
    // Complex exponentials
    double complex term1 = cexp(complex_I * M_kepler * (n + k + m)); // exp(i M (n + k + m))
    double complex term2 = cexp(complex_I * (m + k) * nu);   // exp(i (m + k) nu)
    
    // Real term (1 + e * cos(nu))^(-ell)
    double real_term = pow(1 + e * cos(nu), (-ell));
    // printf("integrand_Z_out real term: %lg \n", real_term);

    double complex prefactor_without_angularparts = - mass_body * cpow(complex_I,(ell - 2)) * pow(2 * omega_gw, (ell + 2)) * tgamma(ell - 1) / (2 * tgamma(2 * ell + 2) * 2 * M_PI) * sqrt(2 * M_PI * tgamma(ell + 3) / tgamma(ell - 1)) *  pow((semi_major * (1-e*e)), (ell));
    // printf("integrand_Z_out prefactor: Re: %lg Im: %lg \n", creal(prefactor_without_angularparts), cimag(prefactor_without_angularparts));

    // Return the product of the complex exponentials and the real term
    // printf("%lg %lg \n", creal(prefactor_without_angularparts * term1 * term2 * real_term), cimag(prefactor_without_angularparts * term1 * term2 * real_term));
    return prefactor_without_angularparts * term1 * term2 * real_term;

    // If prefactor or other terms result in a NaN, we should print them
    if (isnan(creal(prefactor_without_angularparts)) || isnan(cimag(prefactor_without_angularparts)) || isnan(nu) || isnan(E)) {
        printf("NaN detected: Re(prefactor)=%lg, Im(prefactor)=%lg, nu=%lg, E=%lg\n", creal(prefactor_without_angularparts), cimag(prefactor_without_angularparts), nu, E);
        return 0.0 + 0.0 * complex_I;  // Return zero if NaN detected
    }
}

// Trapezoidal rule for numerical integration over complex numbers
double complex trapezoidal_rule(double complex (*f)(double, double, double, int, int, int, int, double, double), 
                                double a, double b, int n, double e, double semi_major, int n_param, int k_param, int m_param, int ell_param, double mass_body, double omega_gw) {
    double h = (b - a) / n;
    double complex sum = 0.5 * f(a, e, semi_major, n_param, k_param, m_param, ell_param, mass_body, omega_gw) + f(b, e, semi_major, n_param, k_param, m_param, ell_param, mass_body, omega_gw);

    for (int i = 1; i < n; i++) {
        sum += f(a + i * h, e, semi_major, n_param, k_param, m_param, ell_param, mass_body, omega_gw);
    }

    return sum * h;
}

int main() {
    double e_inner, semi_major_inner, e_outer, semi_major_outer;  // eccentricity and semi-major axis
    double a = 0.0;  // Lower limit of the integral (M_kepler = 0)
    double b = 2 * M_PI;  // Upper limit of the integral (M_kepler = 2 * pi)
    double  mu_inner = GLOBALPAR_mu_inner, mu_outer = GLOBALPAR_mu_outer; // mass of body
    double complex Z_down_inner, Z_out_inner, Z_down_outer, Z_out_outer;
    double J_inner[3], EQL_inner[3], anc_inner[3], Minv_inner[9], Omega_inner[3], J_outer[3], EQL_outer[3], anc_outer[3], Minv_outer[9], Omega_outer[3];
    double mass = GLOBALPAR_M, spin = GLOBALPAR_astar, omega_nkm_inner, omega_nkm_outer;
    int n = 10000, nl = GLOBALPAR_nl_res, il, i;  // Number of steps for integration
    int ell_param = 2, N_res = GLOBALPAR_N_res, n_res_inner, k_res_inner, m_res_inner, n_res_outer, k_res_outer, m_res_outer;  // start at ell = 2 for gravitational modes

    printf("Enter J_inner components: ");
	scanf("%lg %lg %lg", &J_inner[0], &J_inner[1], &J_inner[2]);
	printf("Enter J_outer components: ");
	scanf("%lg %lg %lg", &J_outer[0], &J_outer[1], &J_outer[2]);

    printf("Mode vector for inner orbit: ");
	scanf("%i %i %i", &n_res_inner, &k_res_inner, &m_res_inner);

	printf("Mode vector for outer orbit: ");
	scanf("%i %i %i", &n_res_outer, &k_res_outer, &m_res_outer);

	CKerr_J2EQL(J_inner, EQL_inner, mass, spin);
	CKerr_EQL2J(EQL_inner, J_inner, mass, spin, anc_inner);
    CKerr_Minverse(J_inner, Minv_inner, mass, spin);
    CKerr_Minv2Omega(Minv_inner, Omega_inner);

	CKerr_J2EQL(J_outer, EQL_outer, mass, spin);
	CKerr_EQL2J(EQL_outer, J_outer, mass, spin, anc_outer);
    CKerr_Minverse(J_outer, Minv_outer, mass, spin);
    CKerr_Minv2Omega(Minv_outer, Omega_outer);

	double I_inner = anc_inner[0];
	double rp_inner = anc_inner[1];
	double ra_inner = anc_inner[2];

	double I_outer = anc_outer[0];
	double rp_outer = anc_outer[1];
	double ra_outer = anc_outer[2];

    e_inner = (ra_inner - rp_inner) / (ra_inner + rp_inner);
    semi_major_inner = (ra_inner + rp_inner) / 2.0;

    e_outer = (ra_outer - rp_outer) / (ra_outer + rp_outer);
    semi_major_outer = (ra_outer + rp_outer) / 2.0;

    printf("J_inner: %lg %lg %lg \n", J_inner[0], J_inner[1], J_inner[2]);
    printf("J_outer: %lg %lg %lg \n", J_outer[0], J_outer[1], J_outer[2]);
    printf("Omega_inner: %lg %lg %lg \n", Omega_inner[0], Omega_inner[1], Omega_inner[2]);
    printf("Omega_outer: %lg %lg %lg \n", Omega_outer[0], Omega_outer[1], Omega_outer[2]);
    printf("Eccentricity and semi-major axis for inner: %lg %lg \n", e_inner, semi_major_inner);
    printf("Eccentricity and semi-major axis for outer: %lg %lg \n", e_outer, semi_major_outer);

    // Compute the complex integral using the trapezoidal rule
    // double complex result = trapezoidal_rule(integrand, a, b, n, e, semi_major, n_param, k_param, m_param, ell_param);
    for (i = -N_res; i <= N_res; i++){
		int i_n_inner = i*n_res_inner;
		int i_k_inner = i*k_res_inner;
		int i_m_inner = i*m_res_inner;
        int i_n_outer = i*n_res_outer;
		int i_k_outer = i*k_res_outer;
		int i_m_outer = i*m_res_outer;
        if((i_n_outer == 0 && i_k_outer == 0 && i_m_outer == 0) || (i_n_inner == 0 && i_k_inner == 0 && i_m_inner == 0))
			continue;
        for (il = 0; (il) < nl; il++){
            int ell = ell_param + il;
            omega_nkm_inner = i_n_inner * Omega_inner[0] + i_k_inner * Omega_inner[1] + i_m_inner * Omega_inner[2];
            omega_nkm_outer = i_n_outer * Omega_outer[0] + i_k_outer * Omega_outer[1] + i_m_outer * Omega_outer[2];
            Z_down_inner = trapezoidal_rule(integrand_Z_down, a, b, n, e_inner, semi_major_inner, i_n_inner, i_k_inner, i_m_inner, ell, mu_inner, omega_nkm_inner);
            Z_out_inner = trapezoidal_rule(integrand_Z_out, a, b, n, e_inner, semi_major_inner, i_n_inner, i_k_inner, i_m_inner, ell, mu_inner, omega_nkm_inner);
            Z_down_outer = trapezoidal_rule(integrand_Z_down, a, b, n, e_outer, semi_major_outer, i_n_outer, i_k_outer, i_m_outer, ell, mu_outer, omega_nkm_outer);
            Z_out_outer = trapezoidal_rule(integrand_Z_out, a, b, n, e_outer, semi_major_outer, i_n_outer, i_k_outer, i_m_outer, ell, mu_outer, omega_nkm_outer);
            printf("(%i, %i, %i, %i) omega_gw_inner = %lg: \n", i_n_inner, i_k_inner, i_m_inner, ell, omega_nkm_inner);
            printf("Inner: Re(Z_down): %lg \t Im(Z_down): %lg \t Re(Z_out): %lg \t Im(Z_out): %lg \n", creal(Z_down_inner), cimag(Z_down_inner), creal(Z_out_inner), cimag(Z_out_inner));
            printf("(%i, %i, %i, %i) omega_gw_outer = %lg: \n", i_n_outer, i_k_outer, i_m_outer, ell, omega_nkm_outer);
            printf("Outer: Re(Z_down): %lg \t Im(Z_down): %lg \t Re(Z_out): %lg \t Im(Z_out): %lg \n", creal(Z_down_outer), cimag(Z_down_outer), creal(Z_out_outer), cimag(Z_out_outer));

        }
    }

    // Print the real and imaginary parts of the result
    // printf("Integral result for (%i, %i, %i, %i): \n", n_param, k_param, m_param, ell_param);
    // printf("Real part: %lg\n", creal(result));
    // printf("Imaginary part: %lg\n", cimag(result));
    

    // double Eccen_anomaly = kepler_eqn(1.5, e);
    // double True_anomaly = true_anomaly(Eccen_anomaly, e);
    // double complex Integrand = integrand(1.5, e, semi_major, n_param, k_param, m_param, ell_param);
    // printf("Eccentric and True anomaly are %lg %lg \n", Eccen_anomaly, True_anomaly);
    // printf("Real and Imaginary parts of the integrand are %lg %lg \n", creal(Integrand), cimag(Integrand));
    return 0;
}

