#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define CKERR_ERRTOL_ANG 1.0e-8
#define CKERR_NSTEP_ANGLE 4096
#define CKERR_NSTEP_ANGLE2 256
#define CKERR_NSTEP_RINTEG 1024

/* Computes the contributions to the stress-energy tensor, int dOmega dr T_{??},
 * for a particle at specified xu[]={r,theta,phi,u_r,u_theta,u_phi}.  Returns
 * the vector of 6 components:
 * Re T_{nn}, Im T_{nn}, Re T_{nmbar}, Im T_{nmbar}, Re T_{mbarmbar}, Im T_{mbarmbar}.
 * The Im T_{nn} component is always zero, but we keep it for consistency.
 */
void CKerr_xu2intT(double *xu, double *T, double M, double astar) {

  double Delta, Sigma, a, costheta, sintheta;
  double gc_tt, gc_tphi, gc_phiphi, u_t, uc_t;
  double udotn, udotmbar[2]; /* udotmbar is complex, need Re & Im components */

  /* Parameters */
  a = M*astar;
  costheta = cos(xu[1]);
  sintheta = sin(xu[1]);
  Delta = xu[0]*(xu[0]-2*M) + a*a;
  Sigma = xu[0]*xu[0]+a*a*costheta*costheta;
  gc_tt = -( (xu[0]*xu[0]+a*a)*(xu[0]*xu[0]+a*a) - Delta*a*a*sintheta*sintheta )/(Delta*Sigma);
  gc_tphi = -2*a*M*xu[0]/Delta/Sigma;
  gc_phiphi = (Delta - a*a*sintheta*sintheta) / (Delta*Sigma*sintheta*sintheta);

  /* Energy */
  u_t = gc_tphi*xu[5] - sqrt(gc_tphi*xu[5]*gc_tphi*xu[5] - gc_tt*(1. + gc_phiphi*xu[5]*xu[5] + Delta/Sigma*xu[3]*xu[3] + xu[4]*xu[4]/Sigma));
  u_t /= -gc_tt;

  /* Contravariant 4-velocity t-component */
  uc_t = gc_tt*u_t + gc_tphi*xu[5];

  /* Dot products with n and mbar */
  udotn = Delta/2./Sigma * ( (xu[0]*xu[0]+a*a)/Delta*u_t + a/Delta*xu[5] - xu[3] );
  udotmbar[0] = ( a*a*sintheta*costheta*u_t + xu[0]*xu[4] + a*costheta/sintheta*xu[5] ) / sqrt(2.) / Sigma;
  udotmbar[1] = ( -a*xu[0]*sintheta*u_t + a*costheta*xu[4] -xu[0]/sintheta*xu[5] ) / sqrt(2.) / Sigma;

  /* Stress-energy components are (u.?)(u.?)/(Sigma uc_t) */
  T[0] = udotn*udotn/(Sigma*uc_t);
  T[1] = 0.;
  T[2] = udotn*udotmbar[0]/(Sigma*uc_t);
  T[3] = udotn*udotmbar[1]/(Sigma*uc_t);
  T[4] = (udotmbar[0]*udotmbar[0]-udotmbar[1]*udotmbar[1])/(Sigma*uc_t);
  T[5] = 2.*udotmbar[0]*udotmbar[1]/(Sigma*uc_t);

  return;
}

/* Take a torus specified by the M-inverse matrix and origin, and a Fourier mode of the torus
 * (qr,qt,mm).  Takes Nr samples in the r-direction and Nt in the
 * theta-direction.  Uses nl theta-modes, i.e. 'l' runs from lmin=max(|mm|,|s|)
 * to lmax=lmin+nl-1.
 *
 * Returns the coefficients c1 and c3 of R1 and R3 on the ingoing and outgoing waves, in
 * the format:  [total length = 4*nl]
 *
 * Coef[] = {Re c1, Im c1, Re c3, Im c3 (for lmin) ...
 *           Re c1, Im c1, Re c3, Im c3 (for lmax)}.
 *
 * and the angular frequency of the mode (in *omega)
 *
 * Also returns the eigenvalues E[nl].
 */
int CKerr_RadialFunc(double *Minv, double *xuorig, double M, double astar, long qr, long qt, long mm,
  long Nr, long Nt, long nl, double *Coef, double *omega, double *E) {

  /* Functions used from other files */
  void CKerr_Minv2Omega(double*,double*);
  int CKerr_J2EQL(double*,double*,double,double);
  int CKerr_EQL2J(double*,double*,double,double,double*);
  long CKerr_SbCoefN(long,long,long,double*,double**,double,double);
  void CKerr_AngleStep(double*,double*,double*,double,double,double,double,long);
  void CKerr_SpheroidalYSLM(double,long,long,long,long,double*,double*,double*);
  void CKerr_IngoingR1(double,double,double,double,long,double,double*);
  void CKerr_OutgoingR3(double,double,double,double,long,double,double*);
  void CKerr_RadialStepMany(double*,double*,double,double,double,long,double,double*,double*,long);
  void CKerr_Radial2nd(double*,double,double,double,long,double,double*,double*);

  /* Variable declarations */
  double a, chi;
  double TOmega[3], xu[6], intT[6], xutemp[6];
  double *b, *Sclm, *dSclm, Spheroidal[3];
  long nb, lmin, il, ir, it, j;
  double phase, cosphase, sinphase;
  double costheta, sintheta;
  double absrho, phaserho, K, Delta, Sigma, rho_iKD[2];
  double temp[2], A0[2], A1[2], A2[2], tempc[2];
  double deriv_KD;
  double r_in, r_out, rhp, rtemp[4];
  double R1_r_in[4], R3_r_out[4], R1_r_out[4], aleph[2];
  double R1loc[6], R3loc[6];
  double xi_norm, beta_norm[2];

  /* Basic setup */
  a = astar*M;
  lmin = abs(mm)<2? 2: abs(mm);

  /* First obtain the angular frequency of the emitted wave */
  CKerr_Minv2Omega(Minv,TOmega);
  *omega = qr*TOmega[0] + qt*TOmega[1] + mm*TOmega[2];
  chi = a*(*omega);

#if 0
fprintf(stderr, "nkm nl indices: %ld %ld %ld %ld\n", qr, qt, mm, nl);
#endif

#if 0
  /* Get periapsis & apoapsis --> ancillary[1], ancillary[2].
   * Then get inner & outer radii to which we integrate the radial Teukolsky equation:
   * this guarantees that we don't unnecessarily repeat the time-consuming integrations from
   * the horizon and from infinity for each point on the orbit.
   */
  CKerr_J2EQL(J,EQL,M,astar);
  CKerr_EQL2J(EQL,J1,M,astar,ancillary);
  printf("incl = %9.5lf (geom) %9.5lf (Carter), rminmax = %13.9lf %13.9lf\n",
    180./M_PI*ancillary[0], 180./M_PI*atan2(sqrt(EQL[1]),EQL[2]), ancillary[1], ancillary[2]);
#endif

  /* Get inner & outer radii to which we integrate the Teukolsky equation (this way we don't have to
   * integrate from horizon and infinity at each timestep)
   */
  CKerr_AngleStep(xuorig,xu,Minv,M_PI,0,M,astar,CKERR_NSTEP_ANGLE);
  rhp = M+M*sqrt(1-astar*astar);
  r_in = 0.99*xuorig[0]+0.01*rhp;
  r_out = 1.01*xu[0];

  /* Now obtain the spheroidal harmonics */
  Sclm = (double*)malloc((size_t)(nl*sizeof(double)));
  dSclm = (double*)malloc((size_t)(nl*sizeof(double)));
  nb = CKerr_SbCoefN(-2,mm,nl,E,&b,chi,CKERR_ERRTOL_ANG);

#if 0
fprintf(stderr, "%19.12le %19.12le %19.12le %19.12le %19.12le (%ld)\n", TOmega[0], TOmega[1], TOmega[2], r_in, r_out, nb);
#endif
  /* Loop over l-modes */
  for(il=0; il<nl; il++) {

    /* Get ingoing solution at r_in and the outgoing solution at r_out */
    CKerr_IngoingR1(r_in,M,astar,*omega,mm,E[il],R1_r_in);
    CKerr_OutgoingR3(r_out,M,astar,*omega,mm,E[il],R3_r_out);

    /* Get the Wronskian coefficient.  Requires getting ingoing solution at r_out. */
    rtemp[0]=r_in; rtemp[1]=0; rtemp[2]=r_out; rtemp[3]=0;
    CKerr_RadialStepMany(rtemp,rtemp+2,M,astar,*omega,mm,E[il],R1_r_in,R1_r_out,CKERR_NSTEP_RINTEG*4);
    aleph[0] = R3_r_out[0]*R1_r_out[2] - R3_r_out[1]*R1_r_out[3] - R1_r_out[0]*R3_r_out[2] + R1_r_out[1]*R3_r_out[3];
    aleph[1] = R3_r_out[0]*R1_r_out[3] + R3_r_out[1]*R1_r_out[2] - R1_r_out[0]*R3_r_out[3] - R1_r_out[1]*R3_r_out[2];
    aleph[0] /= r_out*(r_out-2*M)+a*a;
    aleph[1] /= r_out*(r_out-2*M)+a*a;

#if 0
  fprintf(stderr, "il=%ld aleph = %12.5le,%12.5le\n", il,aleph[0],aleph[1]);
#endif

    /* Initialize modes */
    Coef[4*il  ] = Coef[4*il+1] = Coef[4*il+2] = Coef[4*il+3] = 0.;

    /* Now build sources using a grid on the torus. */
    for(ir=0;ir<Nr;ir++) {
      for(it=0;it<Nt;it++) {
        /* Get current position */
        if (it==0) {
          CKerr_AngleStep(xuorig,xu,Minv,2.*M_PI*ir/(double)Nr, 2.*M_PI*it/(double)Nt, M, astar, CKERR_NSTEP_ANGLE);
        } else {
          for(j=0;j<6;j++) xutemp[j] = xu[j];
          CKerr_AngleStep(xutemp,xu,Minv,0., 2.*M_PI/(double)Nt, M, astar, CKERR_NSTEP_ANGLE2);
        }

        /* Must multiply all sources by the complex exponential, which includes the phase in the
         * Fourier transform as well as in the longitude of the particle at psi_phi=0.
         */
        phase = 2.*M_PI*(qr*ir/(double)Nr+qt*it/(double)Nt) - mm*xu[2];
        cosphase = cos(phase); sinphase = sin(phase);

        /* Get integral of stress-energy tensor, shift phase */
        CKerr_xu2intT(xu,intT,M,astar);
        for(j=0;j<6;j+=2) {
          temp[0] = cosphase*intT[j] - sinphase*intT[j+1];
          temp[1] = sinphase*intT[j] + cosphase*intT[j+1];
          intT[j] = temp[0];
          intT[j+1] = temp[1];
        }

        /* Get spheroidal harmonics.  Spheroidal = {S, L2dS, L1dL2dS} */
        costheta = cos(xu[1]); sintheta = sin(xu[1]);
        CKerr_SpheroidalYSLM(xu[1],-2,mm,nl,nb,b,Sclm,dSclm);
        Spheroidal[0] = Sclm[il];
        Spheroidal[1] = dSclm[il] + (-mm/sintheta + chi*sintheta + 2*costheta/sintheta)*Sclm[il];
        Spheroidal[2] = 2.*(-mm/sintheta + chi*sintheta + costheta/sintheta)*dSclm[il]
                        + (-chi*chi*(costheta*costheta-sintheta*sintheta) - 2*mm*chi + 2*mm*mm/sintheta/sintheta
                           - 6*mm*costheta/sintheta/sintheta - 2. + 4/sintheta/sintheta - E[il])*Sclm[il];

        /* Obtain the A2-coefficient */
        absrho = 1./sqrt(xu[0]*xu[0]+a*a*costheta*costheta);
        phaserho = M_PI - atan2(-a*costheta,xu[0]);
        temp[0] = -Spheroidal[0]/absrho/absrho*cos(4*phaserho);
        temp[1] =  Spheroidal[0]/absrho/absrho*sin(4*phaserho);
        A2[0] = temp[0]*intT[4] - temp[1]*intT[5];
        A2[1] = temp[0]*intT[5] + temp[1]*intT[4];

        /* The A1-coefficient */
        K = *omega*(xu[0]*xu[0]+a*a) - a*mm;
        Delta = xu[0]*(xu[0]-2.*M)+a*a;
        Sigma = xu[0]*xu[0]+a*a*costheta*costheta;
        rho_iKD[0] = absrho*cos(phaserho);
        rho_iKD[1] = absrho*sin(phaserho) - K/Delta;
        A1[0] = -2*( A2[0]*rho_iKD[0] - A2[1]*rho_iKD[1] ); /* mbar mbar part */
        A1[1] = -2*( A2[1]*rho_iKD[0] + A2[0]*rho_iKD[1] );

        temp[0] =  Spheroidal[1]/absrho/absrho/absrho*cos(3*phaserho) + 2.*a*a*sintheta*costheta*Spheroidal[0]/Sigma/absrho/absrho/absrho*cos(3*phaserho);
        temp[1] = -Spheroidal[1]/absrho/absrho/absrho*sin(3*phaserho) - 2.*a*a*sintheta*costheta*Spheroidal[0]/Sigma/absrho/absrho/absrho*sin(3*phaserho);
        temp[0] *= -sqrt(8.)/Delta;
        temp[1] *= -sqrt(8.)/Delta;
        A1[0] += temp[0]*intT[2] - temp[1]*intT[3]; /* n mbar part */
        A1[1] += temp[0]*intT[3] + temp[1]*intT[2];

        /* The A0-coefficient */
        deriv_KD = (Delta*2.*(*omega)*xu[0] - K*2.*(xu[0]-M))/Delta/Delta;
        tempc[0] = K*K/Delta/Delta - 2.*absrho*sin(phaserho)*K/Delta;
        tempc[1] = 2*absrho*cos(phaserho)*K/Delta + deriv_KD;
        temp[0] = ( tempc[0]*cos(4*phaserho) + tempc[1]*sin(4*phaserho))/absrho/absrho * Spheroidal[0];
        temp[1] = (-tempc[0]*sin(4*phaserho) + tempc[1]*cos(4*phaserho))/absrho/absrho * Spheroidal[0];
        A0[0] = temp[0]*intT[4] - temp[1]*intT[5]; /* mbar mbar part */
        A0[1] = temp[0]*intT[5] + temp[1]*intT[4];

        tempc[0] = 2*xu[0]/Sigma * ( Spheroidal[1] - 2*a*a*sintheta*costheta/Sigma*Spheroidal[0] ) * sqrt(8.)/Delta/absrho/absrho/absrho;
        tempc[1] = K/Delta   * ( Spheroidal[1] + 2*a*a*sintheta*costheta/Sigma*Spheroidal[0] ) * sqrt(8.)/Delta/absrho/absrho/absrho;
        temp[0] = -(tempc[0]*cos(3*phaserho) + tempc[1]*sin(3*phaserho));
        temp[1] = -(tempc[1]*cos(3*phaserho) - tempc[0]*sin(3*phaserho));
        A0[0] += temp[0]*intT[2] - temp[1]*intT[3]; /* n mbar part */
        A0[1] += temp[0]*intT[3] + temp[1]*intT[2];

        tempc[0] = -2*a*absrho*sintheta*sin(phaserho)*Spheroidal[1] + Spheroidal[2];
        tempc[1] =  2*a*absrho*sintheta*cos(phaserho)*Spheroidal[1];
        temp[0] = -2*Sigma/absrho/absrho/Delta/Delta*( tempc[0]*cos(2*phaserho) + tempc[1]*sin(2*phaserho) );
        temp[1] = -2*Sigma/absrho/absrho/Delta/Delta*( tempc[1]*cos(2*phaserho) - tempc[0]*sin(2*phaserho) );
        A0[0] += temp[0]*intT[0] - temp[1]*intT[1]; /* n n part */
        A0[1] += temp[0]*intT[1] + temp[1]*intT[0];

        /* The radial functions and their first and second derivatives */
        rtemp[0]=r_in; rtemp[1]=0; rtemp[2]=xu[0]; rtemp[3]=0;
        CKerr_RadialStepMany(rtemp,rtemp+2,M,astar,*omega,mm,E[il],R1_r_in,R1loc,CKERR_NSTEP_RINTEG);
        CKerr_Radial2nd(rtemp+2,M,astar,*omega,mm,E[il],R1loc,R1loc+4);
        rtemp[0]=r_out; rtemp[1]=0; rtemp[2]=xu[0]; rtemp[3]=0;
        CKerr_RadialStepMany(rtemp,rtemp+2,M,astar,*omega,mm,E[il],R3_r_out,R3loc,CKERR_NSTEP_RINTEG);
        CKerr_Radial2nd(rtemp+2,M,astar,*omega,mm,E[il],R3loc,R3loc+4);

        /* Increment the ingoing and outgoing coefficients */
        Coef[4*il  ] += A0[0]*R3loc[0]-A0[1]*R3loc[1]-A1[0]*R3loc[2]+A1[1]*R3loc[3]+A2[0]*R3loc[4]-A2[1]*R3loc[5];
        Coef[4*il+1] += A0[1]*R3loc[0]+A0[0]*R3loc[1]-A1[1]*R3loc[2]-A1[0]*R3loc[3]+A2[1]*R3loc[4]+A2[0]*R3loc[5];
        Coef[4*il+2] += A0[0]*R1loc[0]-A0[1]*R1loc[1]-A1[0]*R1loc[2]+A1[1]*R1loc[3]+A2[0]*R1loc[4]-A2[1]*R1loc[5];
        Coef[4*il+3] += A0[1]*R1loc[0]+A0[0]*R1loc[1]-A1[1]*R1loc[2]-A1[0]*R1loc[3]+A2[1]*R1loc[4]+A2[0]*R1loc[5];

#if 0
        fprintf(stderr, " ***  %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le\n", xu[0], xu[1], xu[2], xu[3], xu[4], xu[5]);
        fprintf(stderr, "      %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le\n", intT[0], intT[1], intT[2], intT[3], intT[4], intT[5]);
        fprintf(stderr, "      %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le\n", A0[0], A0[1], A1[0], A1[1], A2[0], A2[1]);
        fprintf(stderr, "      %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le\n", R1loc[0], R1loc[1], R1loc[2], R1loc[3], R1loc[4], R1loc[5]);
        fprintf(stderr, "      %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le\n", R3loc[0], R3loc[1], R3loc[2], R3loc[3], R3loc[4], R3loc[5]);
#endif
      } /* end for(it) */
    } /* end for(ir) */

    /* Phase-rotate coefficients by aleph */ //TODO: Rewrite to 
    // temp[0] = Coef[4*il  ]; temp[1] = Coef[4*il+1];
    // Coef[4*il  ] = (temp[0]*aleph[0]+temp[1]*aleph[1])/(aleph[0]*aleph[0]+aleph[1]*aleph[1])/(double)(Nr*Nt);
    // Coef[4*il+1] = (temp[1]*aleph[0]-temp[0]*aleph[1])/(aleph[0]*aleph[0]+aleph[1]*aleph[1])/(double)(Nr*Nt);
    // temp[0] = Coef[4*il+2]; temp[1] = Coef[4*il+3];
    // Coef[4*il+2] = (temp[0]*aleph[0]+temp[1]*aleph[1])/(aleph[0]*aleph[0]+aleph[1]*aleph[1])/(double)(Nr*Nt);
    // Coef[4*il+3] = (temp[1]*aleph[0]-temp[0]*aleph[1])/(aleph[0]*aleph[0]+aleph[1]*aleph[1])/(double)(Nr*Nt);

    /* Phase-rotate coefficients by aleph */ // Use this for higher resonance modes
    temp[0] = Coef[4*il  ]; temp[1] = Coef[4*il+1];
    xi_norm =  fmax(aleph[0], aleph[1]);
    beta_norm[0] = aleph[0]/xi_norm;
    beta_norm[1] = aleph[1]/xi_norm;
    Coef[4*il  ] = (temp[0]*beta_norm[0]+temp[1]*beta_norm[1])/((beta_norm[0]*beta_norm[0]+beta_norm[1]*beta_norm[1]) * xi_norm)/(double)(Nr*Nt);
    Coef[4*il+1] = (temp[1]*beta_norm[0]-temp[0]*beta_norm[1])/((beta_norm[0]*beta_norm[0]+beta_norm[1]*beta_norm[1]) * xi_norm)/(double)(Nr*Nt);
    temp[0] = Coef[4*il+2]; temp[1] = Coef[4*il+3];
    Coef[4*il+2] = (temp[0]*beta_norm[0]+temp[1]*beta_norm[1])/((beta_norm[0]*beta_norm[0]+beta_norm[1]*beta_norm[1]) * xi_norm)/(double)(Nr*Nt);
    Coef[4*il+3] = (temp[1]*beta_norm[0]-temp[0]*beta_norm[1])/((beta_norm[0]*beta_norm[0]+beta_norm[1]*beta_norm[1]) * xi_norm)/(double)(Nr*Nt);
#if 0
    fprintf(stderr, "xi_norm=%12.5le beta_norm = %12.5le,%12.5le\n", xi_norm,beta_norm[0],beta_norm[1]);
#endif
  } /* end for(il) */

  free((char*)b);
  free((char*)Sclm);
  free((char*)dSclm);
  return(1);
}

/* Returns the alpha-coefficient, Eq.4.17 of Hughes 2000, for radiation flux into the hole.
 */
double CKerr_HorizonFluxNorm(double M, double astar, double omega, long mm, double eigenlm) {
  double a, rhp, varpi, OmegaH, epsilon, numer, C2, lambda;

  a = M*astar;
  rhp = M+sqrt(M*M-a*a);
  OmegaH = astar/2./rhp;
  varpi = omega-mm*OmegaH;
  lambda = eigenlm - 2*a*mm*omega + a*a*omega*omega - 2.;
  epsilon = sqrt(M*M-a*a)/4./M/rhp;
  numer = 256.*pow(2.*M*rhp,5.)*varpi*(varpi*varpi+4*epsilon*epsilon)*(varpi*varpi+16*epsilon*epsilon)*omega*omega*omega;
  C2 = ((lambda+2)*(lambda+2)+4*a*omega*(mm-a*omega)) * (lambda*lambda+36*mm*a*omega-36*a*a*omega*omega)
       + (2*lambda+3)*(96*a*a*omega*omega-48*mm*a*omega)
       + 144*omega*omega*(M*M-a*a);
  return(numer/C2);
}

/* Returns rate of quanta emission, dE/dt/omega, in the qr,qt,mm mode give the
 * vector of coefficients Coef[4*nl], summed over the nl modes.  The -->infinity
 * and -->horizon parts are returned separately.
 */
void CKerr_QuantaEmitted(double M, double astar, long qr, long qt, long mm, long nl,
  double *E, double *Coef, double omega, double *QuantaInf, double *QuantaHor) {

  long il;
  double alpha;

  *QuantaInf = *QuantaHor = 0.;

  for(il=0;il<nl;il++) {
    alpha = CKerr_HorizonFluxNorm(M,astar,omega,mm,E[il]);
    *QuantaInf += (Coef[4*il+2]*Coef[4*il+2]+Coef[4*il+3]*Coef[4*il+3])/(2.*omega*omega*omega);
    *QuantaHor += (Coef[4*il  ]*Coef[4*il  ]+Coef[4*il+1]*Coef[4*il+1])/(2.*omega*omega*omega)*alpha;
  }

}

/* Computes the resonance amplitude |S(m)|^2 for a specified Lindblad resonance (black hole
 * M,astar; radius of perturber orbit r0, and resonance type mm (+ for ILR, - for OLR).
 * Also requires nl, the number of latitude modes to use; and jr1, the radial action of the
 * test particle (which should go -->0, but we will use this to try to estimate a derivative
 * with respect to the eccentricity, so this should be small but nonzero).
 * Also returns r1 if desired.
 *
 * Returns -1 if jr1 is too large (orbit with that radial excitation plunges into the hole).
 */
double CKerr_LindbladResonanceStrength(double M, double astar, double r0, long mm, long nl,
  double jr1, double *r1) {

  double CKerr_FindLindbladResonance(double,double,double,double);
  int CKerr_getData_CircEq(double,double,double,double*);
  void CKerr_GWScatMatrix(double,double,double,long,double,double*,double*);
  void CKerr_AngleStep(double*,double*,double*,double,double,double,double,long);
  int CKerr_Minverse(double*,double*,double,double);
  int CKerr_TorusOrigin(double*,double*,double,double);
  void CKerr_Minv2Omega(double*,double*);

  long il;
  int msign;
  double myr1, info0[7], info1[7], actionJ0[3], actionJ1[3];
  double omega0, omega1, *E0, *E1, *C0, *C1, *C0p, *C1p;
  double Minv[9], xuorig[6], xu_apocenter[6], tpOmega[3];
  double c[20], P[2], temp[2], prefactor, epsilon, dP[2];

  /* Get resonance location and orbital data */
  msign = mm>0? 1: -1; /* ILR vs OLR */
  myr1 = CKerr_FindLindbladResonance(M,astar,r0,(double)mm);
  CKerr_getData_CircEq(M,astar,r0,info0);
  CKerr_getData_CircEq(M,astar,myr1,info1);
  actionJ0[0] = 0.;
  actionJ0[1] = 0.;
  actionJ0[2] = info0[0];
  actionJ1[0] = jr1;
  actionJ1[1] = 0.;
  actionJ1[2] = info1[0];

  /* Allocate separation eigenvalues and radiated wave amplitudes */
  E0 = (double*)malloc((size_t)(nl*10*sizeof(double)));
  E1 = E0 + nl;
  C0 = E1 + nl;
  C1 = C0 + 4*nl;

  /* Secondary waves.  Since the secondary is in uniform motion, computing the sources
   * at only one position is sufficient.  The origin is already known and we don't need
   * to take time re-computing it.
   */
  CKerr_Minverse(actionJ0,Minv,M,astar);
  xuorig[0] = r0; xuorig[1] = M_PI/2.; xuorig[2] = 0.;
  xuorig[3] = 0.; xuorig[4] = 0.;      xuorig[5] = info0[0];
  CKerr_RadialFunc(Minv,xuorig,M,astar,0,0,abs(mm),1,1,nl,C0,&omega0,E0);

  /* Test particle waves.  Now we need only 1 position in the theta-direction but
   * several in the r-direction; here we take 4.  Also get epsilon, the amplitude
   * of the radial epicycle.
   */
  if (!CKerr_Minverse(actionJ1,Minv,M,astar)) {
    free((char*)E0);
    return(-1);
  }
  CKerr_TorusOrigin(actionJ1,xuorig,M,astar);
  CKerr_RadialFunc(Minv,xuorig,M,astar,-msign,0,abs(mm),4,1,nl,C1,&omega1,E1);
  CKerr_AngleStep(xuorig,xu_apocenter,Minv,M_PI,0,M,astar,CKERR_NSTEP_ANGLE);
  epsilon = (xu_apocenter[0] - xuorig[0])/2.;
  CKerr_Minv2Omega(Minv,tpOmega);
#if 0
  /* check */
  fprintf(stderr, "omega = %14.12lf, %14.12lf; test particle = %14.12lf,%14.12lf jr1=%12.5le epsilon=%12.5le\n", omega0, omega1, tpOmega[0], tpOmega[2], jr1, epsilon);
#endif

  /* Add up contribution from each l-mode */
  P[0] = P[1] = 0.;
  for(il=0;il<nl;il++) {
    /* Scattering parameters; note: c_{13} = c[0]+ic[1], alpha = c[16]. */
    CKerr_GWScatMatrix(M,astar,omega0,abs(mm),E0[il],c,c+16);
    /* Get pointers to this l-mode in the emitted waveforms */
    C0p = C0 + 4*il;
    C1p = C1 + 4*il;

    if (mm>0) {
      /* ILR Case */
      temp[0] =  c[0]*C1p[2] + c[1]*C1p[3] + c[16]*C1p[0];
      temp[1] =  c[1]*C1p[2] - c[0]*C1p[3] - c[16]*C1p[1];
      dP[0] = C0p[0]*temp[0] - C0p[1]*temp[1];
      dP[1] = C0p[0]*temp[1] + C0p[1]*temp[0];
    } else {
      /* OLR Case */
      temp[0] = -c[0]*C1p[0] + c[1]*C1p[1] + C1p[2];
      temp[1] =  c[1]*C1p[0] + c[0]*C1p[1] - C1p[3];
      dP[0] = C0p[2]*temp[0] - C0p[3]*temp[1];
      dP[1] = C0p[2]*temp[1] + C0p[3]*temp[0];
    }
    dP[0] /= -omega0*omega0;
    dP[1] /= -omega0*omega0;

#if 0
    fprintf(stderr, "Bdown(0) = %12.5le,%12.5le; Bout(1) = %12.5le,%12.5le sep = %9.5lf ", C0p[0], C0p[1], C1p[2], C1p[3], E0[il]);
    fprintf(stderr, "c13 = %12.5le,%12.5le alpha = %12.5le; DP = %12.5le,%12.5le eps=%12.5le\n", c[0], c[1], c[16], dP[0], dP[1], epsilon);
#endif
    P[0] += dP[0]; P[1] += dP[1];
  }

  /* Conversion from P --> S.  info1[5] is the epicyclic impedance. */
  prefactor = 2./omega0/info1[5]/epsilon;

  if (r1!=NULL) *r1=myr1;
  free((char*)E0);
  return(prefactor*prefactor*(P[0]*P[0]+P[1]*P[1]));
}

#undef CKERR_ERRTOL_ANG
#undef CKERR_NSTEP_ANGLE
#undef CKERR_NSTEP_ANGLE2
#undef CKERR_NSTEP_RINTEG
