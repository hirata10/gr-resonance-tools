/* Cached-orbit Teukolsky source integration for the self-force mode sum.
 * Moved into the repo 2026-11-09. The source arithmetic is copied from CKerr_RadialFunc (kerrgwem.c):
 *  - the torus samples (CKerr_AngleStep) and the spin-weighted harmonics Y^{-2}_{lm} are computed once per J_dot
 *    (CKerr_OrbitCacheBuild) instead of once per (n,k,m,l);
 *  - the ingoing (R1) and outgoing (R3) radial functions are taken from one stored 4096-step sweep per direction
 *    plus one partial extrapolated-midpoint step to each sample radius, instead of integrating from r_in/r_out
 *    to every sample.
 * Checked against CKerr_RadialFunc: the summed J_dot agrees to <= 6e-13 RELATIVE (nmax 48 run, local J_dot_sf cost study 2026-11-03;
 * production box 2026-11-09: 5e-14 relative, 1e-20 absolute on J_dot ~ 1e-7). Per-mode Z agree to <= 7e-10 relative except on
 * cancelled modes 1e12 below the strongest.
 *
 * HOW IT IS USED (J_dot.c, J_dot_selfforce_openmp):
 *   CKerr_OrbitCache cache;
 *   CKerr_OrbitCacheBuild(Minv, xuorig, M, astar, Nr, Nt, lc, &cache);     once per J_dot (per orbit)
 *   ... for every mode (n,k,m), in parallel, each thread with its own Coef/E arrays:
 *   CKerr_RadialFunc_cached(Minv, &cache, M, astar, n, k, m, nl, Coef, &omega, E);
 *   CKerr_OrbitCacheFree(&cache);
 * The cache is read-only inside CKerr_RadialFunc_cached, so one cache can be shared by all OpenMP threads.
 * Building with -DJDOT_SF_ORIGINAL (make JDOT_SF=original) makes J_dot.c call CKerr_RadialFunc instead; this file is
 * still compiled but unused. That build is the A/B reference for any change made here.
 *
 * WHAT A CALL RETURNS: for mode (qr,qt,mm) = (n,k,m) and each il = 0..nl-1 (l = max(2,|m|) + il),
 *   Coef[4*il+0], Coef[4*il+1] = Re, Im of the horizon ("down") amplitude Z^H_{lmkn},
 *   Coef[4*il+2], Coef[4*il+3] = Re, Im of the infinity ("up")  amplitude Z^inf_{lmkn},
 *   E[il] = spheroidal eigenvalue, *omega = n Omega_r + k Omega_theta + m Omega_phi,
 * exactly as CKerr_RadialFunc returns them. J_dot.c turns these into fluxes with CKerr_HorizonFluxNorm (alpha).
 *
 * FILE LAYOUT: yoff (index into the Y table) -> CKerr_OrbitCacheBuild / CKerr_OrbitCacheFree (the cache)
 *   -> sph_cached (spheroidal harmonics from the table) -> step_to (one partial radial step)
 *   -> CKerr_RadialFunc_cached (the mode: sweeps, source integral, Wronskian normalization). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CKerr.h"

/* These four mirror kerrgwem.c lines 5-8 and MUST be changed together with them: they fix the torus-sample
 * positions (CKERR_NSTEP_ANGLE, _ANGLE2: RK4 steps in CKerr_AngleStep), the spheroidal-expansion tolerance
 * (CKERR_ERRTOL_ANG, fed to CKerr_SbCoefN) and the radial step count (the sweep below uses 4*CKERR_NSTEP_RINTEG
 * steps, the same count CKerr_RadialFunc uses for its Wronskian integration). If the two files disagree, the
 * cached and original paths no longer compute the same numbers. NBMAX is the largest spheroidal basis size
 * (nb) that sph_cached serves from its stack arrays; larger nb falls back to CKerr_SpheroidalYSLM. */
#define CKERR_ERRTOL_ANG 1.0e-8
#define CKERR_NSTEP_ANGLE 4096
#define CKERR_NSTEP_ANGLE2 256
#define CKERR_NSTEP_RINTEG 1024
#define NBMAX 64


/* Layout of the Y table (o->Y and o->dY). For each torus sample s there is one block of B = (lc+1)^2 - 4 doubles
 * holding Y^{-2}_{lm}(theta_s) for every l = 2..lc and every m = -l..l, in the order CKerr_Yslm writes them
 * (a full row of 2l+1 values per l). The row for l starts at offset yoff(l) = sum_{l'=2}^{l-1} (2l'+1) = l^2 - 4,
 * so the entry for (l, m) is  Y[s*B + yoff(l) + l + m]. All m are stored because CKerr_Yslm computes a whole
 * row at once and because different (n,k,m) modes need different m from the same sample. This is what makes
 * the table quadratic in lc: raising lc by 6 (see sph_cached) adds 16*Ns*[(lc+7)^2 - (lc+1)^2] bytes. */
static long yoff(long l) { return l * l - 4; }

/* Builds the per-orbit cache: everything in CKerr_RadialFunc that depends on the orbit but not on the mode.
 * Inputs: Minv and xuorig from CKerr_Minverse and CKerr_TorusOrigin (the orbit's action-angle map and its
 * origin point), the torus grid Nr x Nt, and lc = the largest l for which Y^{-2}_{lm} is tabulated.
 * J_dot.c passes lc = max(2, mmax) + nl - 1, the largest l the mode sum itself uses.
 * Fills o with:
 *   o->XU[6*s .. 6*s+5] = (r, theta, phi, u^r, u^theta, u^phi) at torus sample s = ir*Nt + it, i.e. at the angles
 *                         (psi_r, psi_theta) = (2 pi ir/Nr, 2 pi it/Nt);
 *   o->r_in, o->r_out    = the radii the radial functions are integrated from (the same choice as
 *                         CKerr_RadialFunc: 1% inside pericenter toward the horizon, 1% outside apocenter);
 *   o->Y, o->dY          = Y^{-2}_{lm}(theta_s) and dY/dtheta, l = 2..lc, m = -l..l (layout: see yoff).
 * Memory: 48*Ns + 16*Ns*B bytes, Ns = Nr*Nt, B = (lc+1)^2 - 4. For Nr 28, Nt 14, lc 7: 0.4 MB.
 * Called once per J_dot. Must be followed by CKerr_OrbitCacheFree. */
void CKerr_OrbitCacheBuild(double *Minv, double *xuorig, double M, double astar, long Nr, long Nt, long lc, CKerr_OrbitCache *o) {
  double xpi[6];
  o->Nr = Nr; o->Nt = Nt; o->Ns = Nr * Nt; o->lc = lc; o->B = yoff(lc + 1);
  /* Apocenter point (psi_r = pi) sets r_out; the origin (pericenter, psi_r = 0) sets r_in. */
  CKerr_AngleStep(xuorig, xpi, Minv, M_PI, 0, M, astar, CKERR_NSTEP_ANGLE);
  double rhp = M + M * sqrt(1 - astar * astar);
  o->r_in = 0.99 * xuorig[0] + 0.01 * rhp;
  o->r_out = 1.01 * xpi[0];
  /* Torus samples, stepped exactly as in CKerr_RadialFunc so the positions are bitwise the same: for each ir
   * one long step from the origin to (psi_r, 0) with CKERR_NSTEP_ANGLE RK4 steps, then a chain of short steps
   * of 2 pi/Nt in psi_theta with CKERR_NSTEP_ANGLE2 steps each. The chains are independent across ir, so the
   * parallel loop is over ir; within one ir the it loop is sequential by construction. */
  o->XU = malloc(6 * o->Ns * sizeof(double));
#pragma omp parallel for schedule(dynamic, 1)
  for (long ir = 0; ir < Nr; ir++) {
    double x[6], xt[6];
    for (long it = 0; it < Nt; it++) {
      if (it == 0) {
        CKerr_AngleStep(xuorig, x, Minv, 2. * M_PI * ir / (double)Nr, 2. * M_PI * it / (double)Nt, M, astar, CKERR_NSTEP_ANGLE);
      } else {
        for (int j = 0; j < 6; j++) xt[j] = x[j];
        CKerr_AngleStep(xt, x, Minv, 0., 2. * M_PI / (double)Nt, M, astar, CKERR_NSTEP_ANGLE2);
      }
      for (int j = 0; j < 6; j++) o->XU[6 * (ir * Nt + it) + j] = x[j];
    }
  }
  /* Spin-weight -2 spherical harmonics at every sample's theta, one full m-row per l. These depend only on
   * theta_s, not on the mode, so they are computed once here instead of once per (n,k,m,l) and sample. */
  o->Y = malloc(o->Ns * o->B * sizeof(double));
  o->dY = malloc(o->Ns * o->B * sizeof(double));
#pragma omp parallel for schedule(dynamic, 8)
  for (long s = 0; s < o->Ns; s++)
    for (long l = 2; l <= lc; l++)
      CKerr_Yslm(o->XU[6 * s + 1], l, -2, o->Y + s * o->B + yoff(l), o->dY + s * o->B + yoff(l));
}

/* Releases the three arrays allocated by CKerr_OrbitCacheBuild. */
void CKerr_OrbitCacheFree(CKerr_OrbitCache *o) {
  free(o->XU); free(o->Y); free(o->dY);
  o->XU = o->Y = o->dY = NULL;
}

/* Spin-weighted spheroidal harmonics S_{lm}(theta_s) and dS/dtheta for the nu output modes l = lmin..lmin+nu-1
 * of azimuthal number m, at torus sample s. Same arithmetic and summation order as CKerr_SpheroidalYSLM
 * (kerrmode.c), which it replaces: a spheroidal harmonic is the linear combination
 *     S_i(theta) = sum_{j=0}^{nb-1} b[i*nb + j] * Y^{-2}_{lmin+j, m}(theta)
 * with the coefficients b and the basis size nb from CKerr_SbCoefN (called once per mode by the caller).
 * The only difference from CKerr_SpheroidalYSLM is where the Y's come from: the cache table instead of a fresh
 * CKerr_Yslm evaluation per call.
 *
 * Reason for fallback: the table holds l <= o->lc, which J_dot.c sets to the largest l in the MODE SUM. But the
 * expansion needs MORE spherical harmonics than it has output modes: CKerr_SbCoefN starts at nb = nu + 1 and
 * raises nb until the nu eigenvalues converge to CKERR_ERRTOL_ANG, and the required nb grows with the
 * spheroidicity chi = a*omega. So the expansion reaches up to l = lmin + nb - 1 >= lmin + nu. For |m| = mmax
 * that is at least max(2,mmax) + nl, one past the table (lc = max(2,mmax) + nl - 1), and for smaller |m| it
 * can also exceed lc when chi is large. In that case, or if nb would overflow the stack arrays, the call is
 * forwarded to CKerr_SpheroidalYSLM, which gives the identical result (verified to 15 digits) at the original
 * per-sample cost. */
static void sph_cached(const CKerr_OrbitCache *o, long s, long m, long nu, long nb, double *b, double *S, double *dS) {
  long lmin = labs(m) < 2 ? 2 : labs(m), lmax = lmin + nb - 1;   /* l range of the expansion, not of the sum */
  if (lmax > o->lc || nb > NBMAX) {
    CKerr_SpheroidalYSLM(o->XU[6 * s + 1], -2, m, nu, nb, b, S, dS);   /* exact, just slower */
    return;
  }
  double Ym[NBMAX], dYm[NBMAX];
  for (long l = lmin; l <= lmax; l++) {                          /* gather the m-column of the table */
    Ym[l - lmin] = o->Y[s * o->B + yoff(l) + l + m];
    dYm[l - lmin] = o->dY[s * o->B + yoff(l) + l + m];
  }
  for (long i = 0; i < nu; i++) {                                /* S_i = sum_j b_ij Y_j, as in the original */
    S[i] = dS[i] = 0.;
    for (long j = 0; j < nb; j++) {
      S[i] += b[i * nb + j] * Ym[j];
      dS[i] += b[i * nb + j] * dYm[j];
    }
  }
}

/* Carries a radial solution from a stored sweep node at r0 to a sample radius rs in ONE step of the usual
 * integrator (CKerr_RadialMidptExtrStep: four modified-midpoint passes of 4, 6, 8, 16 substeps, Richardson-
 * combined; kerrmode.c). RRnode and RRout are (Re R, Im R, Re dR/dr, Im dR/dr). The step length |rs - r0| is at
 * most one sweep spacing h, shorter than any step the original per-sample integration takes (it uses 1024 steps
 * of (rs - r_in)/1024, i.e. up to 4h), so this is not an extrapolation and loses no accuracy. r is real here;
 * the integrator itself allows complex r (hh[1] = 0 is the imaginary part of the step). */
static void step_to(double r0, double rs, double M, double astar, double omega, long mm, double eig, double *RRnode, double *RRout) {
  double rr[2] = {r0, 0.}, hh[2] = {rs - r0, 0.};
  if (hh[0] == 0.) { for (int k = 0; k < 4; k++) RRout[k] = RRnode[k]; return; }   /* sample sits on a node */
  CKerr_RadialMidptExtrStep(rr, hh, M, astar, omega, mm, eig, RRnode, RRout);
}

/* The Teukolsky amplitudes Z^H and Z^inf of one mode (qr, qt, mm) = (n, k, m) for the nl values l = max(2,|m|)..
 * max(2,|m|)+nl-1, computed from the cached orbit. Same interface and output as CKerr_RadialFunc (kerrgwem.c)
 * except that the orbit (xuorig, Nr, Nt) is replaced by the cache. Returns 1. Output layout: see the file header.
 *
 * Stages, with the one that differs from CKerr_RadialFunc marked:
 *   1. omega, chi = a*omega, and the spheroidal expansion (eigenvalues E, coefficients b, basis size nb);
 *   2. S and dS/dtheta at every torus sample from the Y table (sph_cached)             [cached: table lookup]
 *   3. per l: boundary values R1(r_in), R3(r_out) from CKerr_IngoingR1/OutgoingR3, then ONE sweep of each
 *      solution across [r_in, r_out] in NSW = 4*CKERR_NSTEP_RINTEG equal steps, every node stored   [cached]
 *      (the R1 sweep also yields R1(r_out), which gives the Wronskian aleph exactly as the original does);
 *   4. per sample: the source coefficients A0, A1, A2 (verbatim from CKerr_RadialFunc) and the radial
 *      functions at the sample radius from the nearest stored node plus one partial step (step_to) [cached]
 *      instead of a fresh 1024-step integration from the boundary;
 *   5. the sum over samples, divided by the Wronskian and by Nr*Nt (the torus average).
 * Memory per call: nl*Ns*16 bytes for S, dS and 2*4*(NSW+1)*8 = 262 kB for the two sweeps.
 * Thread safety: reads the cache, writes only its own heap arrays and the caller's Coef/E, so it can be called
 * concurrently from OpenMP threads sharing one cache (J_dot.c does this). */
int CKerr_RadialFunc_cached(double *Minv, const CKerr_OrbitCache *o, double M, double astar, long qr, long qt, long mm,
  long nl, double *Coef, double *omega, double *E) {
  double a, chi;
  double TOmega[3], xu[6], intT[6];
  double *b, Spheroidal[3];
  long nb, il, ir, it, j, s, Nr = o->Nr, Nt = o->Nt, Ns = o->Ns, NSW = CKERR_NSTEP_RINTEG * 4;
  double phase, cosphase, sinphase;
  double costheta, sintheta;
  double absrho, phaserho, K, Delta, Sigma, rho_iKD[2];
  double temp[2], A0[2], A1[2], A2[2], tempc[2];
  double deriv_KD;
  double r_in = o->r_in, r_out = o->r_out;
  double R1_r_in[4], R3_r_out[4], R1_r_out[4], aleph[2];
  double R1loc[6], R3loc[6];
  double xi_norm, beta_norm[2];

  /* 1. Mode frequency and the spheroidal expansion. CKerr_SbCoefN allocates b (nb x nb, freed below) and
   *    returns the eigenvalues E[0..nl-1]; nb is the basis size it needed for tolerance CKERR_ERRTOL_ANG. */
  a = astar * M;
  CKerr_Minv2Omega(Minv, TOmega);
  *omega = qr * TOmega[0] + qt * TOmega[1] + mm * TOmega[2];
  chi = a * (*omega);
  nb = CKerr_SbCoefN(-2, mm, nl, E, &b, chi, CKERR_ERRTOL_ANG);

  /* 2. Spheroidal harmonics S_l(theta_s) and dS/dtheta for all nl values of l at every sample, from the cached
   *    Y table (or the fallback inside sph_cached). SS[nl*s + il] is sample s, l-index il. */
  double *SS = malloc(nl * Ns * sizeof(double)), *dSS = malloc(nl * Ns * sizeof(double));
  for (s = 0; s < Ns; s++) sph_cached(o, s, mm, nl, nb, b, SS + nl * s, dSS + nl * s);

  /* 3. Storage for the two radial sweeps: R1n[4*i .. 4*i+3] = (Re R1, Im R1, Re R1', Im R1') at r_in + i*h1,
   *    R3n likewise at r_out + i*h3, i = 0..NSW. h1 > 0 runs outward from r_in, h3 = -h1 runs inward from r_out. */
  double *R1n = malloc(4 * (NSW + 1) * sizeof(double)), *R3n = malloc(4 * (NSW + 1) * sizeof(double));
  double h1 = (r_out - r_in) / (double)NSW, h3 = (r_in - r_out) / (double)NSW;

  for (il = 0; il < nl; il++) {
    /* Boundary values: R1 = ingoing at the horizon, started at r_in; R3 = outgoing at infinity, started at
     * r_out (both integrate in from their asymptotic regions inside CKerr_IngoingR1/OutgoingR3). */
    CKerr_IngoingR1(r_in, M, astar, *omega, mm, E[il], R1_r_in);
    CKerr_OutgoingR3(r_out, M, astar, *omega, mm, E[il], R3_r_out);
    for (j = 0; j < 4; j++) { R1n[j] = R1_r_in[j]; R3n[j] = R3_r_out[j]; }
    /* The sweeps: one pass per solution across the whole radial range of the orbit, every node kept.
     * The R1 sweep is the same integration CKerr_RadialFunc does once to get R1(r_out) for the Wronskian;
     * here it also serves every sample. */
    for (long i = 0; i < NSW; i++) {
      double rr[2] = {r_in + i * h1, 0.}, hh[2] = {h1, 0.};
      CKerr_RadialMidptExtrStep(rr, hh, M, astar, *omega, mm, E[il], R1n + 4 * i, R1n + 4 * (i + 1));
      double rr3[2] = {r_out + i * h3, 0.}, hh3[2] = {h3, 0.};
      CKerr_RadialMidptExtrStep(rr3, hh3, M, astar, *omega, mm, E[il], R3n + 4 * i, R3n + 4 * (i + 1));
    }
    /* Wronskian aleph = (R3 R1' - R1 R3') / Delta at r_out, used to normalize the Green's-function integral.
     * Same expression as CKerr_RadialFunc. */
    for (j = 0; j < 4; j++) R1_r_out[j] = R1n[4 * NSW + j];
    aleph[0] = R3_r_out[0]*R1_r_out[2] - R3_r_out[1]*R1_r_out[3] - R1_r_out[0]*R3_r_out[2] + R1_r_out[1]*R3_r_out[3];
    aleph[1] = R3_r_out[0]*R1_r_out[3] + R3_r_out[1]*R1_r_out[2] - R1_r_out[0]*R3_r_out[3] - R1_r_out[1]*R3_r_out[2];
    aleph[0] /= r_out*(r_out-2*M)+a*a;
    aleph[1] /= r_out*(r_out-2*M)+a*a;

    Coef[4*il] = Coef[4*il+1] = Coef[4*il+2] = Coef[4*il+3] = 0.;

    /* 4. Source integral over the torus. Each sample contributes source x radial function, where the
     *    Z^H integrand uses R3 and the Z^inf integrand uses R1 (Green's function construction). */
    for (ir = 0; ir < Nr; ir++) {
      for (it = 0; it < Nt; it++) {
        s = ir * Nt + it;
        for (j = 0; j < 6; j++) xu[j] = o->XU[6 * s + j];

        /* Fourier phase of the mode on the torus, e^{i(n psi_r + k psi_theta)} e^{-i m phi(s)} */
        phase = 2.*M_PI*(qr*ir/(double)Nr+qt*it/(double)Nt) - mm*xu[2];
        cosphase = cos(phase); sinphase = sin(phase);

        /* Stress-energy projections (n n, n mbar, mbar mbar) at the sample, rotated by the phase.
         * From here to the end of the A0 block the code is CKerr_RadialFunc's, line for line. */
        CKerr_xu2intT(xu,intT,M,astar);
        for(j=0;j<6;j+=2) {
          temp[0] = cosphase*intT[j] - sinphase*intT[j+1];
          temp[1] = sinphase*intT[j] + cosphase*intT[j+1];
          intT[j] = temp[0];
          intT[j+1] = temp[1];
        }

        /* Spheroidal = {S, L2 S, L1 L2 S} (the spin-lowering operators applied to S), built from the cached
         * S and dS/dtheta instead of a fresh CKerr_SpheroidalYSLM call. */
        costheta = cos(xu[1]); sintheta = sin(xu[1]);
        Spheroidal[0] = SS[nl*s+il];
        Spheroidal[1] = dSS[nl*s+il] + (-mm/sintheta + chi*sintheta + 2*costheta/sintheta)*SS[nl*s+il];
        Spheroidal[2] = 2.*(-mm/sintheta + chi*sintheta + costheta/sintheta)*dSS[nl*s+il]
                        + (-chi*chi*(costheta*costheta-sintheta*sintheta) - 2*mm*chi + 2*mm*mm/sintheta/sintheta
                           - 6*mm*costheta/sintheta/sintheta - 2. + 4/sintheta/sintheta - E[il])*SS[nl*s+il];

        /* Source coefficients A2, A1, A0 (the coefficients of R, R', R'' in the Teukolsky source; the formalism is
         * arXiv:1010.0759). rho = -1/(r - i a cos theta) is kept as |rho| and its phase. */
        absrho = 1./sqrt(xu[0]*xu[0]+a*a*costheta*costheta);
        phaserho = M_PI - atan2(-a*costheta,xu[0]);
        temp[0] = -Spheroidal[0]/absrho/absrho*cos(4*phaserho);
        temp[1] =  Spheroidal[0]/absrho/absrho*sin(4*phaserho);
        A2[0] = temp[0]*intT[4] - temp[1]*intT[5];
        A2[1] = temp[0]*intT[5] + temp[1]*intT[4];

        K = *omega*(xu[0]*xu[0]+a*a) - a*mm;
        Delta = xu[0]*(xu[0]-2.*M)+a*a;
        Sigma = xu[0]*xu[0]+a*a*costheta*costheta;
        rho_iKD[0] = absrho*cos(phaserho);
        rho_iKD[1] = absrho*sin(phaserho) - K/Delta;
        A1[0] = -2*( A2[0]*rho_iKD[0] - A2[1]*rho_iKD[1] );
        A1[1] = -2*( A2[1]*rho_iKD[0] + A2[0]*rho_iKD[1] );

        temp[0] =  Spheroidal[1]/absrho/absrho/absrho*cos(3*phaserho) + 2.*a*a*sintheta*costheta*Spheroidal[0]/Sigma/absrho/absrho/absrho*cos(3*phaserho);
        temp[1] = -Spheroidal[1]/absrho/absrho/absrho*sin(3*phaserho) - 2.*a*a*sintheta*costheta*Spheroidal[0]/Sigma/absrho/absrho/absrho*sin(3*phaserho);
        temp[0] *= -sqrt(8.)/Delta;
        temp[1] *= -sqrt(8.)/Delta;
        A1[0] += temp[0]*intT[2] - temp[1]*intT[3];
        A1[1] += temp[0]*intT[3] + temp[1]*intT[2];

        deriv_KD = (Delta*2.*(*omega)*xu[0] - K*2.*(xu[0]-M))/Delta/Delta;
        tempc[0] = K*K/Delta/Delta - 2.*absrho*sin(phaserho)*K/Delta;
        tempc[1] = 2*absrho*cos(phaserho)*K/Delta + deriv_KD;
        temp[0] = ( tempc[0]*cos(4*phaserho) + tempc[1]*sin(4*phaserho))/absrho/absrho * Spheroidal[0];
        temp[1] = (-tempc[0]*sin(4*phaserho) + tempc[1]*cos(4*phaserho))/absrho/absrho * Spheroidal[0];
        A0[0] = temp[0]*intT[4] - temp[1]*intT[5];
        A0[1] = temp[0]*intT[5] + temp[1]*intT[4];

        tempc[0] = 2*xu[0]/Sigma * ( Spheroidal[1] - 2*a*a*sintheta*costheta/Sigma*Spheroidal[0] ) * sqrt(8.)/Delta/absrho/absrho/absrho;
        tempc[1] = K/Delta   * ( Spheroidal[1] + 2*a*a*sintheta*costheta/Sigma*Spheroidal[0] ) * sqrt(8.)/Delta/absrho/absrho/absrho;
        temp[0] = -(tempc[0]*cos(3*phaserho) + tempc[1]*sin(3*phaserho));
        temp[1] = -(tempc[1]*cos(3*phaserho) - tempc[0]*sin(3*phaserho));
        A0[0] += temp[0]*intT[2] - temp[1]*intT[3];
        A0[1] += temp[0]*intT[3] + temp[1]*intT[2];

        tempc[0] = -2*a*absrho*sintheta*sin(phaserho)*Spheroidal[1] + Spheroidal[2];
        tempc[1] =  2*a*absrho*sintheta*cos(phaserho)*Spheroidal[1];
        temp[0] = -2*Sigma/absrho/absrho/Delta/Delta*( tempc[0]*cos(2*phaserho) + tempc[1]*sin(2*phaserho) );
        temp[1] = -2*Sigma/absrho/absrho/Delta/Delta*( tempc[1]*cos(2*phaserho) - tempc[0]*sin(2*phaserho) );
        A0[0] += temp[0]*intT[0] - temp[1]*intT[1];
        A0[1] += temp[0]*intT[1] + temp[1]*intT[0];

        /* Radial functions at the sample radius rs (THE cached step). j1 is the sweep node just below rs on the
         * R1 sweep (counted from r_in), j3 the node just above rs on the R3 sweep (counted from r_out; h3 < 0,
         * so (rs - r_out)/h3 is positive). One partial step of length < h from that node reaches rs, and
         * CKerr_Radial2nd supplies R'' from the Teukolsky equation. R1loc/R3loc = (Re R, Im R, Re R', Im R',
         * Re R'', Im R''). CKerr_RadialFunc instead re-integrates from r_in and from r_out in 1024 steps here. */
        {
          double rs = xu[0], rr2[2] = {rs, 0.};
          long j1 = (long)floor((rs - r_in) / h1); if (j1 < 0) j1 = 0; if (j1 > NSW - 1) j1 = NSW - 1;
          step_to(r_in + j1 * h1, rs, M, astar, *omega, mm, E[il], R1n + 4 * j1, R1loc);
          CKerr_Radial2nd(rr2,M,astar,*omega,mm,E[il],R1loc,R1loc+4);
          long j3 = (long)floor((rs - r_out) / h3); if (j3 < 0) j3 = 0; if (j3 > NSW - 1) j3 = NSW - 1;
          step_to(r_out + j3 * h3, rs, M, astar, *omega, mm, E[il], R3n + 4 * j3, R3loc);
          CKerr_Radial2nd(rr2,M,astar,*omega,mm,E[il],R3loc,R3loc+4);
        }

        /* Accumulate the complex products source x R: Coef[0,1] (Z^H) with R3, Coef[2,3] (Z^inf) with R1. */
        Coef[4*il  ] += A0[0]*R3loc[0]-A0[1]*R3loc[1]-A1[0]*R3loc[2]+A1[1]*R3loc[3]+A2[0]*R3loc[4]-A2[1]*R3loc[5];
        Coef[4*il+1] += A0[1]*R3loc[0]+A0[0]*R3loc[1]-A1[1]*R3loc[2]-A1[0]*R3loc[3]+A2[1]*R3loc[4]+A2[0]*R3loc[5];
        Coef[4*il+2] += A0[0]*R1loc[0]-A0[1]*R1loc[1]-A1[0]*R1loc[2]+A1[1]*R1loc[3]+A2[0]*R1loc[4]-A2[1]*R1loc[5];
        Coef[4*il+3] += A0[1]*R1loc[0]+A0[0]*R1loc[1]-A1[1]*R1loc[2]-A1[0]*R1loc[3]+A2[1]*R1loc[4]+A2[0]*R1loc[5];
      }
    }

    /* 5. Normalize: divide the complex sums by the Wronskian aleph (written as xi_norm * beta_norm to avoid
     *    overflow, as in CKerr_RadialFunc) and by Nr*Nt, which turns the sum over samples into the torus
     *    average that the Fourier integral over (psi_r, psi_theta) requires. */
    temp[0] = Coef[4*il  ]; temp[1] = Coef[4*il+1];
    xi_norm =  fmax(aleph[0], aleph[1]);
    beta_norm[0] = aleph[0]/xi_norm;
    beta_norm[1] = aleph[1]/xi_norm;
    Coef[4*il  ] = (temp[0]*beta_norm[0]+temp[1]*beta_norm[1])/((beta_norm[0]*beta_norm[0]+beta_norm[1]*beta_norm[1]) * xi_norm)/(double)(Nr*Nt);
    Coef[4*il+1] = (temp[1]*beta_norm[0]-temp[0]*beta_norm[1])/((beta_norm[0]*beta_norm[0]+beta_norm[1]*beta_norm[1]) * xi_norm)/(double)(Nr*Nt);
    temp[0] = Coef[4*il+2]; temp[1] = Coef[4*il+3];
    Coef[4*il+2] = (temp[0]*beta_norm[0]+temp[1]*beta_norm[1])/((beta_norm[0]*beta_norm[0]+beta_norm[1]*beta_norm[1]) * xi_norm)/(double)(Nr*Nt);
    Coef[4*il+3] = (temp[1]*beta_norm[0]-temp[0]*beta_norm[1])/((beta_norm[0]*beta_norm[0]+beta_norm[1]*beta_norm[1]) * xi_norm)/(double)(Nr*Nt);
  }

  free(SS); free(dSS); free(R1n); free(R3n);
  free((char*)b);                                              /* allocated by CKerr_SbCoefN */
  return(1);
}

