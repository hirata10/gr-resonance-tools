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
double Ckerr_FindISCO(double Ic, double M, double astar);

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
double CKerr_HorizonFluxNorm(double M, double astar, double omega, long mm, double eigenlm);
void CKerr_QuantaEmitted(double M, double astar, long qr, long qt, long mm, long nl,
  double *E, double *Coef, double omega, double *QuantaInf, double *QuantaHor);
double CKerr_LindbladResonanceStrength(double M, double astar, double r0, long mm, long nl,
  double jr1, double *r1);

/* resonance_find.c: finding the apocenter that satisfies the resonance condition for */
/* circular, equitorial orbiting perturber */

void ra_rp_I2EQL(double ra, double *EQL, double rp, double I, double astar, double M);
double Omega_outer_direct(double radius, double M, double spin);
double ra_rp_I2Omega_OuterCirc(int n, int k, int m, double radius, double ra, double rp, double I, double astar, double M);
double ra_rp_I2Omega_generic(int n_inner, int k_inner, int m_inner, int n_outer, int k_outer, double ra_inner, double rp_inner, double I_inner, double ra_outer, double rp_outer, double I_outer, double astar, double M);
double find_resonance_apo_OuterCirc(int n, int k, int m, double radius, double guess1, double guess2, double rp, double I, double astar, double M);

/* J_dot.c: finds the time derivative of the J components for the self-force and tidal field */

void J_dot_selfforce(int nl, int nmax, int kmax, int mmax, double apo, double rp, double radius_outer, double I, double M, double astar, double *J_dot_sf);

void J_dot_tidal(int nl, int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double apo, double rp, double radius_outer, double I, double M, double astar, double theta_res_F, double *J_dot_td);

double J_dot_phi_Kepler(double mu_outer, double r_outer, double apo, double peri, double incline, double Theta_res);

/* Gamma.c: Computes the change in omega for orbits */

int omega_dot(int nl, int n_res_inner, int k_res_inner, int m_res_inner, int n_res_outer, int k_res_outer, double ra_inner, double rp_inner, double I_inner, double ra_outer, double rp_outer, double I_outer, double astar, double M, double radius_outer, double delta_t, double *Gamma);

/* Delta_J.c computes the change in J due to external field (Eq. (12) in arXiv:1905.00030v2) */

void Delta_J_tidal(int nl, int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double ra_inner, double rp_inner, double radius_outer, double I_inner, double ra_outer, double rp_outer, double I_outer, double M, double astar, double theta_res_F, double *Delta_J_tidal);

void Delta_J_tidal2(int nl, int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double ra_inner, double rp_inner, double radius_outer, double I_inner, double ra_outer, double rp_outer, double I_outer, double M, double astar, double theta_res_F, double ang_accel, double mu_outer, double *Delta_J);

/* J2J_dot.c computes the evolution of the J_{i}'s for a specific orbiting body over a some time interval */

int J2Jdot(int nl, int nmax, int kmax, int mmax, double *J_initial, double *J_dot_sf, double M, double astar);

int J2Jdot_component(int nl, int nmax, int kmax, int mmax, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_dot_r, double *J_dot_theta, double *J_dot_phi, double M, double astar);

void rk4_J2Jdot(double t0, int n, double J_r_ini, double J_theta_ini, double J_phi_ini, double *J_r_final, double *J_theta_final, double *J_phi_final, FILE *fptr, double mu_body, double M, double astar);

//int readtxt(char FILENAME, struct data_vals *data);

