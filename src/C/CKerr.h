/* kerrtraj.c: Particle trajectories in Kerr */

int CKerr_EQL2J(double *EQL, double *J, double M, double astar, double *ancillary);
double CKerr_Emin(double Q, double L, double M, double astar);
double CKerr_Emax(double Q, double L, double M, double astar, double *Jrmax);
double CKerr_QLJr2E(double Q, double L, double Jr, double M, double astar, double *Jtheta);
int CKerr_J2EQL(double *J, double *EQL, double M, double astar);
int CKerr_dEQdJ(double *J, double *EQL, double *dEdJ, double *dQdJ, double M, double astar);
int CKerr_MinverseND(double *J, double *Minv, double M, double astar);
int CKerr_Minverse(double *J, double *Minv, double M, double astar);
int CKerr_TorusOrigin(double *J, double *xu, double M, double astar);
void CKerr_Minv2Omega(double *Minv, double *Omega);
void CKerr_FlowVectors(double *xu, double *xudot_H, double *xudot_Q, double M, double astar);
void CKerr_AngleFlowVectors(double *xu, double *xudot_r, double *xudot_theta, double *Minv, double M, double astar);
void CKerr_AngleStep(double *xu_in, double *xu_out, double *Minv, double dpsi_r, double dpsi_t, double M, double astar, long N);
double CKerr_getL_CircEq(double M, double astar, double r);
int CKerr_getData_CircEq(double M, double astar, double r, double *info);
double CKerr_FindLindbladResonance(double M, double astar, double r0, double mm);
int CKerr_FindEQL_IRR(double Ic, double rp, double ra, double *EQL, double M, double astar);
int CKerr_FindEQL_IRCirc(double Ic, double r, double *EQL, double M, double astar);
int CKerr_FindResCirc(int fratio, double Ic, double *J, double M, double astar);
double CKerr_FindISCO(double Ic, double M, double astar);

/* kerrmode.c: Angular & radial mode functions for the Teukolsky equations */

void CKerr_Yslm(double theta, long l, long s, double *Y, double *dYdtheta);
long CKerr_SbCoef(long s, long m, long nb, double *E, double *b, double chi, double errtol);
long CKerr_SbCoefN(long s, long m, long nu, double *E, double **b, double chi, double errtol);
void CKerr_SpheroidalYSLM(double theta, long s, long m, long nu, long nb, double *b, double *S, double *dSdtheta);
double CKerr_r2rstar(double r, double M, double astar);
double CKerr_rstar2r(double rstar, double M, double astar, double *Delta);
void CKerr_Vr(double rstar, double M, double astar, double omega, long mm, double eigenlm,
  double *ReV, double *ImV);
void CKerr_IngoingR1(double rf, double M, double astar, double omega, long mm, double eigenlm,
  double *R1);
void CKerr_AsympR3(double *r, double M, double astar, double omega, long mm,
  double eigenlm, double *R3);
void CKerr_Radial2nd(double *r, double M, double astar, double omega, long mm,
  double eigenlm, double *RR, double *d2Rdr2);
void CKerr_RadialRKStep(double *r, double *h, double M, double astar, double omega, long mm,
  double eigenlm, double *RR_old, double *RR_new);
void CKerr_RadialMidptStep(double *r, double *h, double M, double astar, double omega, long mm,
  double eigenlm, double *RR_old, double *RR_new, long N);
void CKerr_RadialMidptExtrStep(double *r, double *h, double M, double astar, double omega, long mm,
  double eigenlm, double *RR_old, double *RR_new);
void CKerr_OutgoingR3(double rf, double M, double astar, double omega, long mm, double eigenlm,
  double *R3);
void CKerr_RadialStepMany(double *r_old, double *r_new, double M, double astar, double omega, long mm,
  double eigenlm, double *RR_old, double *RR_new, long Nstep);
void CKerr_GWScatMatrix(double M, double astar, double omega, long mm, double eigenlm,
  double *cscat, double *aux);

/* kerrgwem.c: Gravitational wave emission from particles on geodesic orbits in Kerr */

void CKerr_xu2intT(double *xu, double *T, double M, double astar);
int CKerr_RadialFunc(double *Minv, double *xuorig, double M, double astar, long qr, long qt, long mm,
  long Nr, long Nt, long nl, double *Coef, double *omega, double *E);

/* Orbit cache for the self-force mode sum (radialfunc_cached.c): torus samples and Y^{-2}_{lm}, l = 2..lc,
 * computed once per J_dot; CKerr_RadialFunc_cached reuses them for every (n,k,m) and takes both radial
 * functions from one stored 4096-step sweep per direction. Agrees with CKerr_RadialFunc to ~1e-10 relative. */
typedef struct {
  long Nr, Nt, Ns, lc, B;
  double r_in, r_out;
  double *XU, *Y, *dY;
} CKerr_OrbitCache;
void CKerr_OrbitCacheBuild(double *Minv, double *xuorig, double M, double astar, long Nr, long Nt, long lc, CKerr_OrbitCache *o);
void CKerr_OrbitCacheFree(CKerr_OrbitCache *o);
int CKerr_RadialFunc_cached(double *Minv, const CKerr_OrbitCache *o, double M, double astar, long qr, long qt, long mm,
  long nl, double *Coef, double *omega, double *E);
/* Largest l the cached mode function reads from the Y table for this mode and nl (the lc a cache needs for it) */
long CKerr_RadialFunc_cached_lmax(double *Minv, double M, double astar, long qr, long qt, long mm, long nl);
double CKerr_HorizonFluxNorm(double M, double astar, double omega, long mm, double eigenlm);
void CKerr_QuantaEmitted(double M, double astar, long qr, long qt, long mm, long nl,
  double *E, double *Coef, double omega, double *QuantaInf, double *QuantaHor);
double CKerr_LindbladResonanceStrength(double M, double astar, double r0, long mm, long nl,
  double jr1, double *r1);

/* resonance_find.c: finding the apocenter that satisfies the resonance condition for */
/* circular, equitorial orbiting perturber */

void ra_rp_I2EQL(double ra, double *EQL, double rp, double I, double astar, double M);
double R_geodesic(double r, double E, double Q, double L, double a, double M);
double Omega_outer_direct(double radius, double M, double spin);
double ra_rp_I2Omega_OuterCirc(int n, int k, int m, double radius, double ra, double rp, double I, double astar, double M);
double ra_rp_I2Omega_generic(int n_inner, int k_inner, int m_inner, int n_outer, int k_outer, double ra_inner, double rp_inner, double I_inner, double ra_outer, double rp_outer, double I_outer, double astar, double M);
double find_resonance_apo_OuterCirc(int n, int k, int m, double radius, double guess1, double guess2, double rp, double I, double astar, double M);
double find_resonance_apo_OuterGeneric(int n_inner, int k_inner, int m_inner, int n_outer, int k_outer, double guess1, double guess2, double rp_inner, double I_inner, double ra_outer, double rp_outer, double I_outer, double astar, double M);


/* J_dot.c: finds the time derivative of the J components for the self-force and tidal field */

/* LEGACY (outdated; not used by the pipeline): serial self-force for a unit mass ratio on a fixed 28 x 14 grid. See J_dot.c */
void J_dot_selfforce(int nl, int nmax, int kmax, int mmax, double apo, double rp, double radius_outer, double I, double M, double astar, double *J_dot_sf);

/* Returns d(J/mu_body)/dt in physical time (units of M): the rate for unit mass ratio times mu_body = (body mass)/M,
 * like mu_inner and mu_outer in pars.txt. The serial J_dot_selfforce above has no mu_body and always returns the
 * rate per unit mass ratio. */
void J_dot_selfforce_openmp(int nl, int nmax, int kmax, int mmax, int Nr, int Nt, double apo, double rp, double I, double mu_body, double M, double astar, double *J_dot_sf);

// void J_dot_tidal(int nl, int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double apo, double rp, double radius_outer, double I, double M, double astar, double theta_res_F, double *J_dot_td);
void J_dot_tidal(int nl, int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double ra_inner, double rp_inner, double radius_outer, double I_inner, double ra_outer, double rp_outer, double I_outer, double M, double astar, double theta_res_F, double mu_outer, double *J_dot_td);

// OpenMP version for J_dot_tidal
void J_dot_tidal_openmp(
    int nl,
    int N_res,
    int n_res_inner, int n_res_outer,
    int k_res_inner, int k_res_outer,
    int m_res_inner, int m_res_outer,
    double ra_inner, double rp_inner,
    double radius_outer,
    double I_inner, double ra_outer, double rp_outer,
    double I_outer,
    double M, double astar,
    double theta_res_F,
    double mu_inner, double mu_outer,
    double *J_dot_td_inner, double *J_dot_td_outer);
    
/* Tidal resonance amplitudes: the part of the tidal J_dot and of Delta J that does not depend on the resonant angle
 * or on Gamma. Same mode sum as J_dot_tidal_openmp (harmonics s = 1..N_res, counted twice for -s), on an Nr x Nt
 * orbit grid. For each harmonic s, amp[4*(s-1) + 0..3] = R1_inner, R2_inner, R1_outer, R2_outer, where
 * R1 = sum_l termR1/omega_s^3 and R2 = sum_l termR2/omega_s^3 (per unit mass ratio of both bodies), so that
 *   J_dot_td_inner,j(Theta) = 2 sum_s (-s N_inner,j) mu_outer (R1_inner cos(s Theta) + R2_inner sin(s Theta)),
 *   J_dot_td_outer,j(Theta) = 2 sum_s (-s N_outer,j) mu_inner (R1_outer cos(s Theta) + R2_outer sin(s Theta)).
 * l range of harmonic s: l = l0 .. l0 + nl - 1 with l0 = max(2, s|m|) when l_extend = 0 (as J_dot_tidal_openmp and
 * Delta_J_tidal2); when l_extend != 0 it is carried on to l_th(s) + nl - 1, l_th(s) = max(2, s|m|, s|m + k_inner|,
 * s|m + k_outer|) being the first multipole the harmonic couples to (terms are added, never dropped).
 * Optional outputs (pass NULL to skip): amp_l[4*((s-1)*nl_max + il) + 0..3] = the same four numbers for each l
 * (l = l0 + il, nl_max = J_dot_tidal_nl_max(...)), omegagw[2*(s-1) + 0..1] = harmonic frequency s N.Omega of the inner
 * and outer body, *n_used = harmonics summed.
 * The mode functions come from the orbit cache (radialfunc_cached.c) unless built with -DTIDAL_AMP_ORIGINAL.
 * Returns 0, or 1 (bad arguments), 2 (no orbit for the given ra, rp, I), 3 (out of memory). Prints only to stderr. */
int J_dot_tidal_nl_harmonic(int nl, int l_extend, int s, int k_res_inner, int k_res_outer, int m_res);
int J_dot_tidal_nl_max(int nl, int l_extend, int N_res, int k_res_inner, int k_res_outer, int m_res);
int J_dot_tidal_amplitudes_openmp(int nl, int l_extend, int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double ra_inner, double rp_inner, double I_inner, double ra_outer, double rp_outer, double I_outer, double M, double astar, int Nr, int Nt, double *amp, double *amp_l, double *omegagw, int *n_used);

/* Tidal J_dot of both bodies at resonant angle theta_res_F from the amplitudes above: d(J/mu_body)/dt in physical time
 * (inner ~ mu_outer, outer ~ mu_inner), like J_dot_tidal_openmp */
void J_dot_tidal_from_amplitudes(int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, const double *amp, double theta_res_F, double mu_inner, double mu_outer, double *J_dot_td_inner, double *J_dot_td_outer);

// Compute Keplerian torque
double J_dot_phi_Kepler(double mu_outer, double r_outer, double apo, double peri, double incline, double Theta_res);

/* Gamma.c: Computes the change in omega for orbits */

/* LEGACY (outdated; not used by the pipeline or the test tools): unit mass ratios, mode sum fixed at nmax = kmax = mmax = 1. See Gamma.c */
int omega_dot(int nl, int n_res_inner, int k_res_inner, int m_res_inner, int n_res_outer, int k_res_outer, double ra_inner, double rp_inner, double I_inner, double ra_outer, double rp_outer, double I_outer, double astar, double M, double radius_outer, double delta_t, double *Gamma);

// Eq. (17) of arXiv:2207.07733 using dOmega/dJ from kerrphase.c, J_dots supplied
// The J_dots are d(J/mu_body)/dt in physical time with the mass ratios already in them (as J_dot_selfforce_openmp returns them); nothing is rescaled here
int Gamma_dOmegadJ(int n_inner, int k_inner, int n_outer, int k_outer, int m_res, double *J_inner, double *J_outer, double *J_dot_inner, double *J_dot_outer, double M, double astar, double *Gamma);

// Same as Gamma_dOmegadJ, with the J_dots computed by J_dot_selfforce_openmp for the mass ratios mu_inner, mu_outer
int Gamma_dOmegadJ_selfforce(int nl, int nmax, int kmax, int mmax, int Nr, int Nt, int n_inner, int k_inner, int n_outer, int k_outer, int m_res, double *J_inner, double *J_outer, double mu_inner, double mu_outer, double M, double astar, double *J_dot_inner, double *J_dot_outer, double *Gamma);

/* The one place where the Gamma used by the Delta J jump prescription (Eq. (38) of arXiv:2207.07733) is chosen.
 * source = GAMMA_SOURCE_SUPPLIED: Gamma[2] = Gamma_supplied (e.g. gamma of potential_resonances_*.txt); Gamma[0], Gamma[1]
 *          and the J_dots are set to NaN, and the self-force arguments are not used.
 * source = GAMMA_SOURCE_EQ17:     Gamma_dOmegadJ_selfforce, i.e. Eq. (17) with the self-force of both bodies.
 * Gamma is in physical time (1/M^2) in both cases. Returns 0, or 1 for an unknown source (Gamma = NaN). See Gamma.c */
#define GAMMA_SOURCE_SUPPLIED 0
#define GAMMA_SOURCE_EQ17 1
int Gamma_resonance(int source, double Gamma_supplied, int nl_self, int nmax, int kmax, int mmax, int Nr_self, int Nt_self, int n_inner, int k_inner, int n_outer, int k_outer, int m_res, double *J_inner, double *J_outer, double mu_inner, double mu_outer, double M, double astar, double *J_dot_inner, double *J_dot_outer, double *Gamma);

/* Delta_J.c computes the change in J due to external field (Eq. (12) in arXiv:1905.00030v2) */

/* LEGACY (outdated; nothing calls it): uses the legacy omega_dot with unit mass ratios. The pipeline uses Delta_J_tidal2. See Delta_J.c */
void Delta_J_tidal(int nl, int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double ra_inner, double rp_inner, double radius_outer, double I_inner, double ra_outer, double rp_outer, double I_outer, double M, double astar, double theta_res_F, double *Delta_J_tidal);

void Delta_J_tidal2(int nl, int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double ra_inner, double rp_inner, double radius_outer, double I_inner, double ra_outer, double rp_outer, double I_outer, double M, double astar, double theta_res_F, double ang_accel, double mu_outer, double *Delta_J);

/* Jump of the actions of both bodies across a resonance, Eq. (38) of arXiv:2207.07733, from the tidal resonance
 * amplitudes of J_dot_tidal_amplitudes_openmp (J_dot.c). theta_res_F = resonant angle at the crossing, ang_accel = Gamma
 * in physical time (1/M^2). Delta_J_inner = Delta(J_inner/mu_inner) ~ mu_outer, Delta_J_outer = Delta(J_outer/mu_outer) ~ mu_inner.
 * Returns 0, or 1 if Gamma is zero or not finite (both outputs are then NaN). */
int Delta_J_from_amplitudes(int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, const double *amp, double theta_res_F, double ang_accel, double mu_inner, double mu_outer, double *Delta_J_inner, double *Delta_J_outer);

/* Amplitudes + jump in one call: the OpenMP, two-body counterpart of Delta_J_tidal2 on an Nr x Nt orbit grid
 * (l_extend as in J_dot_tidal_amplitudes_openmp; 0 gives the l range of Delta_J_tidal2).
 * Returns the status of J_dot_tidal_amplitudes_openmp or of Delta_J_from_amplitudes (outputs are NaN unless it is 0). */
int Delta_J_tidal_openmp(int nl, int l_extend, int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double ra_inner, double rp_inner, double I_inner, double ra_outer, double rp_outer, double I_outer, double M, double astar, double theta_res_F, double ang_accel, double mu_inner, double mu_outer, int Nr, int Nt, double *Delta_J_inner, double *Delta_J_outer);

/* J2J_dot.c computes the evolution of the J_{i}'s for a specific orbiting body over a some time interval */
/* LEGACY (outdated; not used by the pipeline or the test tools): J2Jdot, J2Jdot_component, rk4_J2Jdot and rk4_J2Jdot_restart
 * are built on the serial J_dot_selfforce (no mass ratio; the RK4 time is mu_body * t). The pipeline uses the two
 * _openmp routines. See J2J_dot.c */

int J2Jdot(int nl, int nmax, int kmax, int mmax, double *J_initial, double *J_dot_sf, double M, double astar);

int J2Jdot_component(int nl, int nmax, int kmax, int mmax, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_dot_r, double *J_dot_theta, double *J_dot_phi, double M, double astar);

int J2Jdot_component_openmp(int nl, int nmax, int kmax, int mmax, int Nr, int Nt, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_dot_r, double *J_dot_theta, double *J_dot_phi, double mu_body, double M, double astar);

void rk4_J2Jdot(double *t, int n, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_r_final, double *J_theta_final, double *J_phi_final, FILE *fptr, double mu_body, double M, double astar);

void rk4_J2Jdot_restart(double t_start, double *t, int i_start, int n, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_r_final, double *J_theta_final, double *J_phi_final, int RESTART, double mu_body, double M, double astar);

void rk4_J2Jdot_restart_openmp(double t_start, double *t, int i_start, int n, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_r_final, double *J_theta_final, double *J_phi_final, int RESTART, double mu_body, int Nr, int Nt, double M, double astar);

/* Actions and angles of both bodies (resonance-aware evolution): dJ/dt from the self-force, d(psi)/dt = Omega(J). */
/* The angles are measured from the torus origin of CKerr_TorusOrigin: pericentre, equatorial crossing with u_theta > 0, phi = 0. */
// Right-hand side for one body: J_dot[3] = d(J/mu_body)/dt, Omega[3] = d(psi)/dt. Returns 1 if J is not a bound orbit
int J_psi_rhs_openmp(const double *J, double mu_body, int Nr, int Nt, double M, double astar, double *J_dot, double *Omega);
// One RK4 step of y[12] = (J_inner, psi_inner, J_outer, psi_outer); time step as in rk4_J2Jdot_restart_openmp, not longer than dt_max if dt_max > 0
// Outputs y_new[12], *dt, Omega_new[6] (inner, outer). Returns 1 if a body has no bound orbit
int rk4_J_psi_step_openmp(const double *y, double dt_max, double mu_inner, double mu_outer, int Nr, int Nt, double M, double astar, double *y_new, double *dt, double *Omega_new);

/* kerrphase.c computes changes in phase crossing a resonance */

int invertMatrix(double *Minverse, double *M_matrix);

// dOmega_dJ[3 * i + k] = dOmega_i/dJ_k from the contour integrals
int J2dOmega_dJ(double *J, double M, double astar, double *dOmega_dJ);
