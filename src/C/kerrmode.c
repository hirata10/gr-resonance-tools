#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define CKERR_NSTEP_YSLM 28
#define CKERR_MODE_OVERFLOW 1e49
#define CKERR_NSTEP_TORTOISE 1024
#define CKERR_INNER_BDY (-30.0)
#define CKERR_RMODE_STEP 0.003
#define CKERR_RMODE_STEP_EXTR 0.03
#define CKERR_OUTER_BDY (1200.0)

/* Returns vector of spin-weighted spherical harmonics: Y^s_{lm}(theta),
 * for m = -l ... +l.  If requested, also gives dY/dtheta.
 */
void CKerr_Yslm(double theta, long l, long s, double *Y, double *dYdtheta) {

  long i, N, c, m1, j1, j2, j3;
  double *R1, *R1old;
  double theta1;

  N = 2*l+1; /* size of matrices */
  c = 2*l*(l+1); /* center entry of matrix */

  /* These matrices are 1D N*N arrays with row-major ordering */
  R1 = (double*)malloc((size_t)(N*N*sizeof(double))); /* Allocation */
  R1old = (double*)malloc((size_t)(N*N*sizeof(double)));

  /* Put latitutde in -pi .. +pi */
  theta -= 2.*M_PI*floor(theta/2./M_PI+0.5);

  /* Build small rotation matrix first, through 2nd order. */
  theta1 = theta*pow(0.5,CKERR_NSTEP_YSLM);
  for(i=0;i<N*N;i++) R1[i] = 0.;
  /* Diagonal entries */
  for(m1=-l;m1<=l;m1++) R1[(N+1)*m1+c] = -( l*(l+1) - m1*m1 ) *theta1*theta1/4.;
  /* Single off-diagonal entries */
  for(m1=-l;m1<l;m1++) {
    R1[N*(m1)+(m1+1)+c] = sqrt((l+m1+1)*(l-m1)) * theta1/2.;
    R1[N*(m1+1)+(m1)+c] = -sqrt((l+m1+1)*(l-m1)) * theta1/2.;
  }
  /* Double off-diagonal entries -- via multiplication */
  for(i=0;i<N-2;i++) {
    R1[N*i+i+2] = R1[N*(i+2)+i] = R1[N*i+i+1]*R1[N*(i+1)+i+2]/2.;
  }

  /* Sequentially double the rotation angle */
  for(i=0;i<CKERR_NSTEP_YSLM;i++) {
    for(j3=0;j3<N*N;j3++) R1old[j3] = R1[j3];
    /* Matrix multiplication: R1 = 2*R1old + R1old*R1old */
    for(j1=0;j1<N;j1++) for(j2=0;j2<=j1 && j1+j2<N;j2++) {
      R1[N*j1+j2] = 2. * R1old[N*j1+j2];
      for(j3=0;j3<N;j3++) R1[N*j1+j2] += R1old[N*j1+j3] * R1old[N*j3+j2];
      /* Use symmetry relations to get the other triangle */
      R1[N*j2+j1] = R1[N*j1+j2] * ((2*N+j1-j2)%2? -1: 1);
    }
    /* The other quadrants also follow by symmetry: R_{m1,m2} = R_{-m2,-m1} */
    for(j1=1;j1<N;j1++) for(j2=N-j1;j2<N;j2++) R1[N*j1+j2] = R1[N*(N-1-j2)+(N-1-j1)];
  }

  /* Add the identity back */
  for(j1=0;j1<N;j1++) R1[N*j1+j1] += 1;

#if 0
  /* For debugging purposes: write R1 */
  for(j1=0;j1<N;j1++) {
    for(j2=0;j2<N;j2++) printf(" %15.12lf", R1[N*j1+j2]);
    printf("\n");
  }
#endif

  for(j1=0;j1<N;j1++) {
    Y[j1] = R1[N*(l-s)+j1] * sqrt(N/4./M_PI); /* Output the desired row */
    if ((j1+l)%2==1) Y[j1]=-Y[j1];
  }

  /* Colatitude derivatives */
  if (dYdtheta!=NULL) for(j1=0;j1<N;j1++) {
    dYdtheta[j1] = ( sqrt((l-s+1)*(l+s))*R1[N*(l-s+1)+j1] - sqrt((l-s)*(l+s+1))*R1[N*(l-s-1)+j1] ) / 2. * sqrt(N/4./M_PI);
    if ((j1+l)%2==1) dYdtheta[j1]=-dYdtheta[j1];
  }

  free((char*)R1); /* Cleanup */
  free((char*)R1old);
  return;
}

/* Returns the bjchi coefficients for the spheroidal harmonics in terms of the spin-weighted
 * spherical harmonics.  The number of vertical basis modes is nb (lmax = |m,s|+nb-1), and the 
 * number of usable modes is returned (i.e. leakage <errtol).  The eigenvalues are returned to
 * E[nb], and the basis components are returned to b[nb*nb] (format: b[nb*i+j] is the component
 * of the jth vertical mode starting at j=0 for l=|m,s| in the ith eigenvector).
 */
long CKerr_SbCoef(long s, long m, long nb, double *E, double *b, double chi, double errtol) {

  int flag=1;
  int ji=0;
  double *band0, *band1, *band2, *A, err=1.;
  long i, i2, J, lmin, p1, p2, n_use;
  double phi, cosphi, sinphi, temp, norm;
  double *bp1old, *bp2old;

  lmin = abs(m);
  if (lmin<abs(s)) lmin=abs(s);

  if (nb<2) {fprintf(stderr, "Error: CKerr_SbCoef: nb must be >=2.\n"); return(0);}

  /* The components of the eigenvalue system, negative of Hughes 2000 Eq A4 so that we will have
   * a positive definite system.
   */
  band0 = (double*)malloc((size_t)(nb*sizeof(double)));
  band1 = (double*)malloc((size_t)(nb*sizeof(double)));
  band2 = (double*)malloc((size_t)(nb*sizeof(double)));
  bp1old = (double*)malloc((size_t)(nb*sizeof(double)));
  bp2old = (double*)malloc((size_t)(nb*sizeof(double)));
  A = (double*)malloc((size_t)(nb*nb*sizeof(double)));

  for(i=0;i<nb;i++) band0[i] = band1[i] = band2[i] = 0.;

  /* Build the nonrotating component */
  for(i=0;i<nb;i++) {
    J = lmin+i;
    band0[i] = J*(J+1.);
  }

  /* Build the cos theta component */
  for(i=0;i<nb;i++) {
    J = lmin+i;
    /* 2s * chi * sqrt{(2J+1)/(2J+3)} * <J,1;m,0|J+1,m> * <J,1;-s,0|J+1,-s> */
    band1[i] += 2. * s * chi / sqrt(1.+1./(J+0.5)) * sqrt((J+1.-m)*(J+1.+m)*(J+1.-s)*(J+1.+s)) / (2*J+1.) / (J+1.);
    /* 2s * chi * <J,1;m,0|J,m> * <J,1;-s,0|J,-s> */
    if (J>=1)
      band0[i] += 2. * s * chi * (-m*s/(J*(J+1.)));
  }

  /* Build the cos^2theta component */
  for(i=0;i<nb;i++) {
    J = lmin+i;
    /* -chi^2 * [1/3*delta_{JJ'} + 2/3*sqrt{(2J+1)/(2J'+1)}*<J,2;m,0|J',m>*<J,2;-s,0|J',-s> */
    if (J>=1)
      band0[i] += -chi*chi/3. - 2*chi*chi/3. * (3*m*m-J*(J+1.)) * (3*s*s-J*(J+1.)) / (J*(J+1.)*(2*J-1.)*(2*J+3.));
    if (J>=1)
      band1[i] -= 2*chi*chi * (-m*s) * sqrt((J+1.-m)*(J+1.+m)*(J+1.-s)*(J+1.+s)) / (J*(2*J+1.)*(J+1.)*(J+2.)) / sqrt(1.+1./(J+0.5));
    band2[i] -= 2.*chi*chi * sqrt((J+1.-m)*(J+1.+m)*(J+1.-s)*(J+1.+s)) * sqrt((J+2.-m)*(J+2.+m)*(J+2.-s)*(J+2.+s))
                / ((J+2.)*(2*J+1.)*(2*J+2.)*(2*J+3.)) / sqrt(1.+2./(J+0.5));
  }

  /* From here on, we build the A-matrix, which is in the basis determined by b.
   * Successive Jacobi rotations will be performed on both the basis and the A-matrix.
   * The original matrix is A_orig = b^T A b, where b_{i,i2} = b[i*nb+i2].
   */
  for(i=0;i<nb;i++) {
    for(i2=0;i2<nb;i2++) {
      A[i*nb+i2] = i2-i==2? band2[i]: i2-i==1? band1[i]: i2-i==0? band0[i]: i2-i==-1? band1[i2]: i2-i==-2? band2[i2]: 0;
      b[i*nb+i2] = i==i2? 1: 0;
    }
  }

  while(err>errtol) {

    /* Cycle through pivot elements: rotate the p1-p2 plane. */
    for(p1=1;p1<nb;p1++) for(p2=0;p2<p1;p2++) {
      if (A[p2*nb+p2]>A[p1*nb+p1]) {
        phi = 0.5*atan2(2*A[p1*nb+p2], A[p2*nb+p2]-A[p1*nb+p1]);
      } else {
        phi = 0.5*atan2(-2*A[p1*nb+p2], -A[p2*nb+p2]+A[p1*nb+p1]);
      }
      cosphi = cos(phi); sinphi = sin(phi);

      /* Now rotate the basis vectors */
      for(i2=0;i2<nb;i2++) {
        bp1old[i2] = b[p1*nb+i2];
        bp2old[i2] = b[p2*nb+i2];
      }
      for(i2=0;i2<nb;i2++) {
        b[p1*nb+i2] =  bp1old[i2] * cosphi - bp2old[i2] * sinphi;
        b[p2*nb+i2] =  bp1old[i2] * sinphi + bp2old[i2] * cosphi;
      }
      /* ... and the A-matrix */
      for(i2=0;i2<nb;i2++) {
        bp1old[i2] = A[p1*nb+i2]; 
        bp2old[i2] = A[p2*nb+i2];  
      }
      for(i2=0;i2<nb;i2++) {
        A[p1*nb+i2] =  bp1old[i2] * cosphi - bp2old[i2] * sinphi;
        A[p2*nb+i2] =  bp1old[i2] * sinphi + bp2old[i2] * cosphi;
      }
      for(i2=0;i2<nb;i2++) {
        bp1old[i2] = A[i2*nb+p1]; 
        bp2old[i2] = A[i2*nb+p2];  
      }
      for(i2=0;i2<nb;i2++) {
        A[i2*nb+p1] =  bp1old[i2] * cosphi - bp2old[i2] * sinphi;
        A[i2*nb+p2] =  bp1old[i2] * sinphi + bp2old[i2] * cosphi;
      }
    }

    err = 0.;
    for(i=1;i<nb;i++) for(i2=0;i2<i;i2++) err = err>fabs(A[i*nb+i2])? err: fabs(A[i*nb+i2]);

    ji++; if (ji>20) {fprintf(stderr, "Too many Jacobi iterations.\n"); break;}
  }

  /* If necessary, sort the eigenvalues using a bubble sort */
  while(flag) {
    flag=0;
    for(i=0;i<nb-1;i++) if (A[i*nb+i]>A[(i+1)*nb+i+1]+errtol) {
      flag=1;
      p1=i; p2=i+1;
      /* Now swap the basis vectors */
      for(i2=0;i2<nb;i2++) {
        temp = b[p1*nb+i2];
        b[p1*nb+i2] = b[p2*nb+i2];
        b[p2*nb+i2] = temp;
      }
      /* ... and the A-matrix */
      for(i2=0;i2<nb;i2++) {
        temp = A[p1*nb+i2]; 
        A[p1*nb+i2] = A[p2*nb+i2];  
        A[p2*nb+i2] = temp;
      }
      for(i2=0;i2<nb;i2++) {
        temp = A[i2*nb+p1]; 
        A[i2*nb+p1] = A[i2*nb+p2];  
        A[i2*nb+p2] = temp;
      }
    }
  }

  /* Determine how many of the eigenmodes are accurate by investigating the 'missing' part of
   * the eigenvalue operator: we have truncated it so measure norm of (E-E_trunc)|v>.
   */
  for(i=0;i<nb;i++) {
    n_use = i+1;
    norm = (band1[nb-1]*b[i*nb+nb-1]+band2[nb-2]*b[i*nb+nb-2])*(band1[nb-1]*b[i*nb+nb-1]+band2[nb-2]*b[i*nb+nb-2])
           + band2[nb-1]*band2[nb-1]*b[i*nb+nb-1]*b[i*nb+nb-1];
    norm = sqrt(norm);

    E[i] = A[i*nb+i];
    if (norm>errtol) break;
  }  

  free((char*)band0);
  free((char*)band1);
  free((char*)band2);
  free((char*)bp1old);
  free((char*)bp2old);
  free((char*)A);

  /* Normalize the eigenvectors by sqrt(2pi) to get normalization consistent with the literature */
  for(i=0;i<nb;i++) {
    for(i2=0;i2<nb;i2++) {
      b[i*nb+i2] *= sqrt(2.*M_PI);
    }
  }

  return(n_use);
}

/* Sets up a table of spheroidal harmonics for a given s, m, and number of usable harmonics.
 * Returns the number of basis modes actually used (nb), as well as a vector of energies, and a table
 * of b-coefficients in the format: b[nb*i+j] is the component of the jth vertical mode starting
 * at j=0 for l=|m,s| in the ith eigenvector.  j = 0..nb-1.
 *
 * b is allocated by this function since the desired size is not known beforehand.  The proper way
 * to call this is to declare:
 * double *my_b;
 * nb = CKerr_SbCoefN( ..., &b, ... );  <-- note you must save nb since there is no other way
 *                                          to know the number of basis modes!
 * ... {lines that use b} ...
 * free((char*)b);
 */
long CKerr_SbCoefN(long s, long m, long nu, double *E, double **b, double chi, double errtol) {
  long i, nb, nu_ret=0;
  double *E_temp;

  *b = (double*)malloc(sizeof(double));
  for(nb=nu+1; nu_ret<=nu; nb++) {
    E_temp = (double*)malloc((size_t)(nb*sizeof(double)));
    free((char*)(*b));
    *b = (double*)malloc((size_t)(nb*nb*sizeof(double)));
    nu_ret = CKerr_SbCoef(s,m,nb,E_temp,*b,chi,errtol);
    for(i=0;i<nu;i++) E[i] = E_temp[i];
    free((char*)E_temp);
  }
  nb--;

  return(nb);
}

/* Returns the spheroidal harmonics Slm and their theta-derivatives, taking in the number of modes (nu)
 * and the b-coefficients, and returning a vector of S[nu], where S[j]: l=min(|m|,|s|)+j.  The
 * theta-derivative is also returned.
 */
void CKerr_SpheroidalYSLM(double theta, long s, long m, long nu, long nb, double *b, double *S, double *dSdtheta) {

  long lmin, lmax, l, Nmax, i, j;
  double *Ytemp, *dYtemp, *Ym, *dYm;

  /* Set multipoles, array sizes & allocation */
  lmin = abs(m); if (lmin<abs(s)) lmin=abs(s);
  lmax = lmin + nb - 1;
  Nmax = 2*lmax + 1;
  Ytemp = (double*)malloc((size_t)(Nmax*sizeof(double)));
  dYtemp = (double*)malloc((size_t)(Nmax*sizeof(double)));
  Ym = (double*)malloc((size_t)(nb*sizeof(double)));
  dYm = (double*)malloc((size_t)(nb*sizeof(double)));

  /* Build vector of spherical harmonics */
  for(l=lmin; l<=lmax; l++) {
    CKerr_Yslm(theta, l, s, Ytemp, dYtemp);
    Ym[l-lmin] = Ytemp[l+m];
    dYm[l-lmin] = dYtemp[l+m];
  }

  /* Combine with the b-coefficients to get the spheroidal harmonics */
  for(i=0;i<nu;i++) {
    S[i] = dSdtheta[i] = 0.;
    for(j=0;j<nb;j++) {
      S[i] += b[i*nb+j] * Ym[j];
      dSdtheta[i] += b[i*nb+j] * dYm[j];
    }
  }

  free((char*)Ytemp); free((char*)dYtemp);
  free((char*)Ym); free((char*)dYm);
  return;
}

/* Tortoise coordinate transformations: r --> rstar */
double CKerr_r2rstar(double r, double M, double astar) {
  double a, rm, rp, rstar, ratio;

  a = M*astar;
  rm = M*(1.-sqrt(1.-astar*astar));
  rp = M*(1.+sqrt(1.-astar*astar));
  ratio = 1./sqrt(1.-astar*astar);

  if (r<=rm) return(-CKERR_MODE_OVERFLOW);

  rstar = r + ratio*rp*log((r-rp)/(2.*M)) - ratio*rm*log((r-rm)/(2.*M));
  return(rstar);
}

/* Tortoise coordinate transformation: rstar --> r.  If requested, returns Delta,
 * since at rstar << -M we may have r~rp, Delta~0 and numerical error would result
 * from computing Delta(r).
 */
double CKerr_rstar2r(double rstar, double M, double astar, double *Delta) {

  int i;
  double a, rm, rp, r, u, ratio, rstar_trial, du, D;

  a = M*astar;
  rm = M*(1.-sqrt(1.-astar*astar));
  rp = M*(1.+sqrt(1.-astar*astar));
  ratio = 1./sqrt(1.-astar*astar);

  /* Newton-Raphsom search as a function of u = r-rp.
   * Since the test function is concave down, an iteration must take us to below
   * the true u unless it hits the bound.
   */
  u = rstar<M? M: rstar;
  for(i=0;i<CKERR_NSTEP_TORTOISE;i++) {
    r = rp+u;
    rstar_trial = r + ratio*rp*log(u/(2.*M)) - ratio*rm*log((r-rm)/(2.*M)) - rstar;
    D = u*(r-rm);
    du = -rstar_trial * D / (r*r+a*a);
    if (du<-0.999*u) {
      du=-0.999*u; /* Bound to prevent run to u<0, which is unphysical. */
      if (i==CKERR_NSTEP_TORTOISE-1) fprintf(stderr, "Warning: CKerr_rstar2r exits without converging.\n");
    }
    u += du;
    if (u<=0) {
      if (Delta!=NULL) *Delta = 2.*M*exp(rstar - rp + ratio*rm*log((rp-rm)/(2.*M)))/ratio/rp;
      return(rp);
    }
  }

  /* Outputs: Delta, if requested, and r.  Note Delta is computed directly using u. */
  if (Delta!=NULL) *Delta = u*(rp-rm+u);
  return(rp+u);
}

/* Potential for the radial Teukolsky equation in terms of rstar (NOT r!)
 * Requires as input the frequency omega and eigenvalue eigenlm, as well as
 * azimuthal quantum number mm.
 * Returns both the real and imaginary parts.
 */
void CKerr_Vr(double rstar, double M, double astar, double omega, long mm, double eigenlm,
  double *ReV, double *ImV) {

  double r, Delta, a, lambda, K;

  a = astar*M;
  r = CKerr_rstar2r(rstar, M, astar, &Delta);

  lambda = eigenlm - 2*a*mm*omega + a*a*omega*omega - 2.;
  K = (r*r+a*a)*omega - mm*a;

  *ReV = -K*K/Delta + lambda;
  *ImV = -4*(r-M)*K/Delta + 8*omega*r;
  return;
}

/* Returns the ingoing radial solution R1 for a wave interior to any matter sources.
 * Integrates to r=rf.
 * Normalized to Delta^2 e^{-i varpi rstar} as rstar --> -infty.
 * Returns length-4 vector: Re(R1), Im(R1), Re(dR1/dr), Im(dR1/dr).
 * [Note: The integrator works in rstar space, but converts to r at the end.]
 */
void CKerr_IngoingR1(double rf, double M, double astar, double omega, long mm, double eigenlm,
  double *R1) {

  int j,k;
  long i, N;
  double rstar_final, rstar_initial, drstar, rstar;
  double varpi, Gamma, a, OmegaH, rhp, Delta, r=0;
  double rk4frac[] = {0, 0.5, 0.5, 1};
  double rk4wt[] = {1./6., 1./3., 1./3., 1./6.};
  double dR1[4], R1test[4], R1prime[4];
  double V[2], coef0, coef1;

  /* Parameters */
  a = astar*M;
  rhp = M*(1.+sqrt(1.-astar*astar));
  OmegaH = astar/2./(1.+sqrt(1.-astar*astar));
  varpi = omega - mm*OmegaH;
  Gamma = sqrt(1.-astar*astar)/rhp;
  rstar_initial = CKERR_INNER_BDY / Gamma;
  rstar_final = CKerr_r2rstar(rf,M,astar);
  N = 2 + (long)floor((rstar_final-rstar_initial)/CKERR_RMODE_STEP*Gamma);
  if (N<=0) fprintf(stderr, "Warning: CKerr_IngoingR1: N=%ld\n", N);
  drstar = (rstar_final-rstar_initial)/N;

  /* Initial conditions */
  CKerr_rstar2r(rstar_initial,M,astar,&Delta);
  R1[0] = Delta*Delta*cos(varpi*rstar_initial);
  R1[1] = -Delta*Delta*sin(varpi*rstar_initial);
  R1[2] = 2.*Gamma*R1[0] + varpi*R1[1];
  R1[3] = 2.*Gamma*R1[1] - varpi*R1[0];
  rstar = rstar_initial;

  /* Integrate by RK4 algorithm */
  for(i=0;i<N;i++) {

    for(k=0;k<4;k++) dR1[k] = R1prime[k] = 0.;
    for(j=0;j<4;j++) {
      rstar = rstar_initial + drstar*(i + rk4frac[j]);
      for(k=0;k<4;k++) R1test[k] = R1[k] + drstar * R1prime[k] * rk4frac[j];

      /* Now get potential and dif eq coefs */
      r = CKerr_rstar2r(rstar,M,astar,&Delta);
      CKerr_Vr(rstar,M,astar,omega,mm,eigenlm,V,V+1);
      coef0 = (r*r+a*a)*(r*r+a*a)/Delta;
      coef1 = 2.*r - 4./Delta*(r-M)*(r*r+a*a);

      /* Derivatives */
      R1prime[0] = R1test[2];
      R1prime[1] = R1test[3];
      R1prime[2] = (V[0]*R1test[0]-V[1]*R1test[1]-coef1*R1test[2])/coef0;
      R1prime[3] = (V[1]*R1test[0]+V[0]*R1test[1]-coef1*R1test[3])/coef0;

      for(k=0;k<4;k++) dR1[k] += rk4wt[j] * drstar * R1prime[k];
    }

    for(k=0;k<4;k++) R1[k] += dR1[k];
  }

  /* Convert derivatives from d/dr* --> d/dr by multiplying by dr* / dr=(r^2+a^2)/Delta.
   * Note that the RK4 integrator has already saved the "correct" final r and Delta.
   */
  R1[2] *= (r*r+a*a)/Delta;
  R1[3] *= (r*r+a*a)/Delta;
}

/* Generates the asymptotic form of the vacuum R3 solution (purely outgoing at r=infty).
 * Takes in complex r[2] (Re, Im) and outputs R3[4] (ReR3, ImR3, ReR3', ImR3').
 */
void CKerr_AsympR3(double *r, double M, double astar, double omega, long mm,
  double eigenlm, double *R3) {

  double rstar[2], rsqa[2], Delta[2], drstardr[2];
  double a, rm, rp, ratio;
  double rinv[2], rinv2[2], rinv3[2];
  double lambda, C1[2], C2[2], lnS1[4];

  a = M*astar;
  rm = M*(1.-sqrt(1.-astar*astar));
  rp = M*(1.+sqrt(1.-astar*astar));
  ratio = 1./sqrt(1.-astar*astar);

  /* Real and imag parts of rstar; the branch cut is taken on the real axis at r<rp,
   * i.e. far away from the regime where the asymptotic solution is valid.
   */
  rstar[0] = r[0] + ratio*rp*log(sqrt((r[0]-rp)*(r[0]-rp)+r[1]*r[1])/(2.*M))
                  - ratio*rm*log(sqrt((r[0]-rm)*(r[0]-rm)+r[1]*r[1])/(2.*M));
  rstar[1] = r[1] + ratio*rp*atan2(r[1],r[0]-rp) - ratio*rm*atan2(r[1],r[0]-rm);

  /* ... and the derivative of r*, (r^2+a^2)/Delta */
  rsqa[0] = r[0]*r[0]-r[1]*r[1]+a*a;
  rsqa[1] = 2*r[0]*r[1];
  Delta[0] = rsqa[0] - 2*M*r[0];
  Delta[1] = rsqa[1] - 2*M*r[1];
  drstardr[0] = (rsqa[0]*Delta[0] + rsqa[1]*Delta[1])/(Delta[0]*Delta[0]+Delta[1]*Delta[1]);
  drstardr[1] = (rsqa[1]*Delta[0] - rsqa[0]*Delta[1])/(Delta[0]*Delta[0]+Delta[1]*Delta[1]);

  /* The inverse powers of 1/r */
  rinv[0] = r[0]/(r[0]*r[0]+r[1]*r[1]);
  rinv[1] = -r[1]/(r[0]*r[0]+r[1]*r[1]);
  rinv2[0] = rinv[0]*rinv[0] - rinv[1]*rinv[1];
  rinv2[1] = 2*rinv[0]*rinv[1];
  rinv3[0] = rinv[0]*rinv2[0] - rinv[1]*rinv2[1];
  rinv3[1] = rinv[0]*rinv2[1] + rinv[1]*rinv2[0];

  /* The coefficients of the asymptotic series -- Press & Teukolsky 1973 ApJ 185,649, Eq D15 */
  lambda = eigenlm - 2*a*mm*omega + a*a*omega*omega - 2.;
  C1[0] = 0.;
  C1[1] = lambda/(2.*omega) + a*mm;
  C2[0] = 1.5*a*a - 1.5*a*mm/omega - lambda/(4.*omega*omega);
  C2[1] = a*mm*M + 1.5*M/omega;

  /* The argument of the exponential in P&T Eq D15 and its r-derivative */
  lnS1[0] = -omega*rstar[1] + 3*log(sqrt(r[0]*r[0]+r[1]*r[1])/M)
              + (C1[0]*rinv[0]-C1[1]*rinv[1]) + (C2[0]*rinv2[0]-C2[1]*rinv2[1]);
  lnS1[1] = omega*rstar[0] + 3*atan2(r[1],r[0])
              + (C1[0]*rinv[1]+C1[1]*rinv[0]) + (C2[0]*rinv2[1]+C2[1]*rinv2[0]);
  lnS1[2] = -omega*drstardr[1] + 3.*rinv[0]
              - (C1[0]*rinv2[0]-C1[1]*rinv2[1]) - 2.*(C2[0]*rinv3[0]-C2[1]*rinv3[1]);
  lnS1[3] = omega*drstardr[0] + 3.*rinv[1]
              - (C1[0]*rinv2[1]+C1[1]*rinv2[0]) - 2.*(C2[0]*rinv3[1]+C2[1]*rinv3[0]);

  /* Return S1: The complex exponential */
  R3[0] = exp(lnS1[0]) * cos(lnS1[1]);
  R3[1] = exp(lnS1[0]) * sin(lnS1[1]);
  R3[2] = R3[0]*lnS1[2] - R3[1]*lnS1[3];
  R3[3] = R3[0]*lnS1[3] + R3[1]*lnS1[2];
  return;
}

/* Generates the second derivative d^2RR/dr^2 given RR[4] (ReR, ImR, ReR', ImR')
 * for the radial Teukolsky equation.  This function allows r to be complex.
 */
void CKerr_Radial2nd(double *r, double M, double astar, double omega, long mm,
  double eigenlm, double *RR, double *d2Rdr2) {

  double a, lambda, rsqa[2], Delta[2], K[2], K__D[2], Kp[2], V[2], DRpp[2];
  double absDelta2;

  /* Parameters */
  a = M*astar;
  lambda = eigenlm - 2*a*mm*omega + a*a*omega*omega - 2.;

  /* Useful functions */
  rsqa[0] = r[0]*r[0]-r[1]*r[1]+a*a;
  rsqa[1] = 2*r[0]*r[1];
  Delta[0] = rsqa[0] - 2*M*r[0];
  Delta[1] = rsqa[1] - 2*M*r[1];
  K[0] = omega*rsqa[0] - a*mm;
  K[1] = omega*rsqa[1];

  /* K / Delta */
  absDelta2 = Delta[0]*Delta[0]+Delta[1]*Delta[1];
  K__D[0] = (K[0]*Delta[0]+K[1]*Delta[1])/absDelta2;
  K__D[1] = (K[1]*Delta[0]-K[0]*Delta[1])/absDelta2;
  /* K + 4i(r-M) */
  Kp[0] = K[0] - 4*r[1];
  Kp[1] = K[1] + 4*(r[0]-M);
  /* Potential */
  V[0] = -K__D[0]*Kp[0] +K__D[1]*Kp[1] - 8*omega*r[1] + lambda;
  V[1] = -K__D[0]*Kp[1] -K__D[1]*Kp[0] + 8*omega*r[0];
  /* Delta * d^2R/dr^2 */
  DRpp[0] = V[0]*RR[0]-V[1]*RR[1] + 2.*( (r[0]-M)*RR[2] - r[1]*RR[3] );
  DRpp[1] = V[0]*RR[1]+V[1]*RR[0] + 2.*( (r[0]-M)*RR[3] + r[1]*RR[2] );

  /* ... and finally, the second derivative! */
  d2Rdr2[0] = (DRpp[0]*Delta[0]+DRpp[1]*Delta[1])/absDelta2;
  d2Rdr2[1] = (DRpp[1]*Delta[0]-DRpp[0]*Delta[1])/absDelta2;
  return;
}

/* Takes an initial condition for the radial Teukolsky equation and uses the RK4 method
 * to integrate it from r --> r+h (1 step).
 */
void CKerr_RadialRKStep(double *r, double *h, double M, double astar, double omega, long mm,
  double eigenlm, double *RR_old, double *RR_new) {

  int j,k;
  double rk4frac[] = {0, 0.5, 0.5, 1};
  double rk4wt[] = {1./6., 1./3., 1./3., 1./6.};
  double RRtest[4], hRRprime[4];
  double rtest[2], d2Rdr2[2];

  for(k=0;k<4;k++) {
    hRRprime[k] = 0.;
    RR_new[k] = RR_old[k];
  }
  for(j=0;j<4;j++) {
    /* RK4 is based on derivatives computed at 4 trial points -- loop over these */
    rtest[0] = r[0] + h[0]*rk4frac[j];
    rtest[1] = r[1] + h[1]*rk4frac[j];
    for(k=0;k<4;k++) RRtest[k] = RR_old[k] + hRRprime[k]*rk4frac[j];
    CKerr_Radial2nd(rtest, M, astar, omega, mm, eigenlm, RRtest, d2Rdr2);

    /* Get h times derivative of RR, then increment output */
    hRRprime[0] = h[0]*RRtest[2] - h[1]*RRtest[3];
    hRRprime[1] = h[0]*RRtest[3] + h[1]*RRtest[2];
    hRRprime[2] = h[0]*d2Rdr2[0] - h[1]*d2Rdr2[1];
    hRRprime[3] = h[0]*d2Rdr2[1] + h[1]*d2Rdr2[0];

    for(k=0;k<4;k++) RR_new[k] += rk4wt[j] * hRRprime[k];
  }
  return;
}

/* Modified midpoint step for the radial Teukolsky equation, using N sub-steps.
 * Requires N>=2.
 */
void CKerr_RadialMidptStep(double *r, double *h, double M, double astar, double omega, long mm,
  double eigenlm, double *RR_old, double *RR_new, long N) {

  int j;
  long i;
  double subh[2], *RRtemp, RRderiv[4], rtemp[2], dtemp[2];

  if (N<2) {fprintf(stderr, "Error: CKerr_RadialMidptStep: N=%ld<2\n", N); exit(1);}
  subh[0] = h[0]/N; subh[1] = h[1]/N;
  RRtemp = (double*)malloc((size_t)(4*(N+1)*sizeof(double)));
  for(j=0;j<4;j++) RRtemp[j] = RR_old[j];

  /* Cycle through midpoint steps, NR C 2nd ed 16.3.2 */
  for(i=0; i<=N; i++) {
    rtemp[0] = r[0]+i*subh[0]; rtemp[1] = r[1]+i*subh[1];

    /* Obtain the derivative function, and then multiply by subh */
    CKerr_Radial2nd(rtemp, M, astar, omega, mm, eigenlm, RRtemp+4*i, dtemp);
    RRderiv[0] = RRtemp[4*i+2] * subh[0] - RRtemp[4*i+3] * subh[1];
    RRderiv[1] = RRtemp[4*i+3] * subh[0] + RRtemp[4*i+2] * subh[1];
    RRderiv[2] = dtemp[0] * subh[0] - dtemp[1] * subh[1];
    RRderiv[3] = dtemp[1] * subh[0] + dtemp[0] * subh[1];
    
    if (i==0) { /* initial step */
      for(j=0;j<4;j++) RRtemp[4+j] = RRtemp[j] + RRderiv[j];
    } else if (i==N) { /* final step: output */
      for(j=0;j<4;j++) RR_new[j] = 0.5*(RRtemp[4*N-4+j]+RRtemp[4*N+j]+RRderiv[j]);
    } else { /* "interior" step case */
      for(j=0;j<4;j++) RRtemp[4*i+4+j] = RRtemp[4*i-4+j] + 2*RRderiv[j];
    }
  }

  free((char*)RRtemp);
  return;
}

/* Extrapolated modified-midpoint step */
void CKerr_RadialMidptExtrStep(double *r, double *h, double M, double astar, double omega, long mm,
  double eigenlm, double *RR_old, double *RR_new) {

  int j;
  double temp_RR_new[16];

  /* Possible steps and weights:
   *  (16,8,4) --> (64,-20,1)/45
   *  (16,8,6) --> (1792,-880,243)/1155
   *  (16,8,6,4) --> (28672,-17600,6561,-308)/17325
   */

  /* Current extrapolation based on N = 4, 6, 8, and 16 steps */
  CKerr_RadialMidptStep(r,h,M,astar,omega,mm,eigenlm,RR_old,temp_RR_new   , 4);
  CKerr_RadialMidptStep(r,h,M,astar,omega,mm,eigenlm,RR_old,temp_RR_new+ 4, 6);
  CKerr_RadialMidptStep(r,h,M,astar,omega,mm,eigenlm,RR_old,temp_RR_new+ 8, 8);
  CKerr_RadialMidptStep(r,h,M,astar,omega,mm,eigenlm,RR_old,temp_RR_new+12,16);

  for(j=0;j<4;j++) {
    RR_new[j] =  (28672.*temp_RR_new[12+j] -17600.*temp_RR_new[8+j] +6561.*temp_RR_new[4+j] -308.*temp_RR_new[j])/17325.;
  }
  return;
}

/* Returns the outgoing radial solution R3 for a wave exterior to any matter sources.
 * Integrates to r=rf.
 * Normalized to r^3 e^{i omega rstar} as r --> +infty.
 * Returns length-4 vector: Re(R1), Im(R1), Re(dR1/dr), Im(dR1/dr).
 * Integrates along an arced trajectory in the complex r-plane (see notes for motivation).
 */
void CKerr_OutgoingR3(double rf, double M, double astar, double omega, long mm, double eigenlm,
  double *R3) {

  long i, k, N;
  double Re_init, Re_f, dRe;
  double r[2], rprev[2], h[2];
  double RR[4], RRprev[4];
  double rp, rsep;

  rp = M*(1. + sqrt(1.-astar*astar));
  rsep = 3.*M;
  if (rsep<1./fabs(omega)) rsep = 1./fabs(omega);

  /* Set up grid in real part of r.  Integrate to Re_f first, i.e. if rf<rsep then we
   * arrive at the real axis at rsep.
   */
  Re_f = rf>rsep? rf: rsep;
  Re_init = fabs(omega)<M? CKERR_OUTER_BDY/fabs(omega): CKERR_OUTER_BDY*M;

  /* Initial point from asymptotic solution */
  r[0] = Re_init;
  r[1] = 3*log(r[0]/Re_f)/omega; if (r[1]>r[0]) r[1]=r[0]; if (r[1]<-r[0]) r[1]=-r[0];
  CKerr_AsympR3(r,M,astar,omega,mm,eigenlm,RR);

  /* Step inward */
  while (r[0]>Re_f) {
    rprev[0]=r[0]; rprev[1]=r[1];
    r[0] -= CKERR_RMODE_STEP_EXTR * (r[0]>1./fabs(omega)? 1./fabs(omega): r[0]);
    r[1] = 3*log(r[0]/Re_f)/omega; if (r[1]>r[0]) r[1]=r[0]; if (r[1]<-r[0]) r[1]=-r[0];
    h[0] = r[0]-rprev[0]; h[1] = r[1]-rprev[1];
    for(k=0;k<4;k++) RRprev[k] = RR[k];
#if 0
    CKerr_RadialRKStep(rprev,h,M,astar,omega,mm,eigenlm,RRprev,RR);
#endif
    CKerr_RadialMidptExtrStep(rprev,h,M,astar,omega,mm,eigenlm,RRprev,RR);
  }
  rprev[0]=r[0]; rprev[1]=r[1];
  r[0] = Re_f; r[1] = 0.;
  h[0] = r[0]-rprev[0]; h[1] = r[1]-rprev[1];
  for(k=0;k<4;k++) RRprev[k] = RR[k];
  CKerr_RadialMidptExtrStep(rprev,h,M,astar,omega,mm,eigenlm,RRprev,RR);

  /* If we are going inside rsep, integrate the rest of the way. */
  if (rf<rsep) {
    N = 2+(long)floor(log((rsep-rp)/(rf-rp))/CKERR_RMODE_STEP);
    if (N<=0) fprintf(stderr, "Warning: CKerr_OutgoingR3: [#2] N=%ld\n", N);
    for(i=0;i<N;i++) {
      rprev[0]=r[0]; rprev[1]=r[1];
      r[0] = rp + (rsep-rp)*pow((rf-rp)/(rsep-rp),(i+1)/(double)N); r[1] = 0.;
      h[0] = r[0]-rprev[0]; h[1] = r[1]-rprev[1];
      for(k=0;k<4;k++) RRprev[k] = RR[k];
      CKerr_RadialRKStep(rprev,h,M,astar,omega,mm,eigenlm,RRprev,RR);
    }
  }

  for(k=0;k<4;k++) R3[k]=RR[k];
  return;
}

/* Takes an initial condition for the radial Teukolsky equation and uses the extrapolation method
 * to integrate it from r_old --> r_new (Nstep steps).
 */
void CKerr_RadialStepMany(double *r_old, double *r_new, double M, double astar, double omega, long mm,
  double eigenlm, double *RR_old, double *RR_new, long Nstep) {

  long i;
  double h[2], r[2], RR[4], RRX[4];

  h[0] = (r_new[0]-r_old[0])/(double)Nstep;
  h[1] = (r_new[1]-r_old[1])/(double)Nstep;

  RR[0] = RR_old[0];
  RR[1] = RR_old[1];
  RR[2] = RR_old[2];
  RR[3] = RR_old[3];
  for(i=0;i<Nstep;i++) {
    r[0] = r_old[0]+i*h[0];
    r[1] = r_old[1]+i*h[1];
    CKerr_RadialMidptExtrStep(r,h,M,astar,omega,mm,eigenlm,RR,RRX);
    RR[0] = RRX[0];
    RR[1] = RRX[1];
    RR[2] = RRX[2];
    RR[3] = RRX[3];
  }
  RR_new[0] = RR[0];
  RR_new[1] = RR[1];
  RR_new[2] = RR[2];
  RR_new[3] = RR[3];
}

/* Obtains the scattering matrix, cscat[16]: (Re,Im) components of c13, c14, c23, c24
 * and then the inverse: (Re, Im) components of c31, c32, c41, c42.
 *
 * Also returns auxiliary data, if requested:
 * aux[0] = alpha (downward wave flux normalizaton)
 * aux[1] = alpha' (upward wave flux normalization)
 * aux[2] = power reflection coefficient (past inf --> future inf)
 * aux[3] = power transmission coefficient (past inf --> future hor)
 */
void CKerr_GWScatMatrix(double M, double astar, double omega, long mm, double eigenlm,
  double *cscat, double *aux) {

  int j;
  double rnorm;
  double a, rhp, varpi, OmegaH, epsilon, numer, C2, lambda;
  double alpha, alphaprime;
  double rmatch, R1m[4], R3m[4], aleph[2];
  double R1old[4], R1new[4], R3asym[4], r_old[2], r_new[2], beta[2], detc[2], idetc[2];

  /* First we can compute the Teukolsky-Starobinsky constant */
  a = M*astar;
  rhp = M+sqrt(M*M-a*a);
  OmegaH = astar/2./rhp;
  varpi = omega-mm*OmegaH;
  lambda = eigenlm - 2*a*mm*omega + a*a*omega*omega - 2.;
  epsilon = sqrt(M*M-a*a)/4./M/rhp;
  C2 = ((lambda+2)*(lambda+2)+4*a*omega*(mm-a*omega)) * (lambda*lambda+36*mm*a*omega-36*a*a*omega*omega)
       + (2*lambda+3)*(96*a*a*omega*omega-48*mm*a*omega)
       + 144*omega*omega*(M*M-a*a);

  /* Auxiliary constants */
  numer = 256.*pow(2.*M*rhp,5.)*varpi*(varpi*varpi+4*epsilon*epsilon)*(varpi*varpi+16*epsilon*epsilon)*omega*omega*omega;
  alpha = numer/C2;
  if (aux!=NULL) aux[0] = alpha;
  alphaprime = (varpi*varpi+16*epsilon*epsilon)*omega*omega*omega
               /pow(2.*M*rhp,3.)/varpi/(varpi*varpi+epsilon*epsilon)/(varpi*varpi+epsilon*epsilon);
  if (aux!=NULL) aux[1] = alphaprime;

  /* Now start integrating solutions: first get ingoing and outgoing modes to meet at some radius */
  rmatch = fabs(omega)>0.4/M? 2.5*M: 1./fabs(omega);
  CKerr_IngoingR1(2.5*M,M,astar,omega,mm,eigenlm,R1m);
  if (fabs(omega)<0.4/M) {
    for(j=0;j<4;j++) R1old[j] = R1m[j];
    r_old[0] = 2.5*M; r_old[1] = 0.;
    while (r_old[0] < rmatch) {
      r_new[0] = r_old[0] * (1. + CKERR_RMODE_STEP_EXTR); r_new[1] = 0.;
      CKerr_RadialStepMany(r_old,r_new,M,astar,omega,mm,eigenlm,R1old,R1new,2);
      for(j=0;j<2;j++) r_old[j]=r_new[j];
      for(j=0;j<4;j++) R1old[j]=R1new[j];
    }
    r_new[0] = rmatch; r_new[1] = 0.;
    CKerr_RadialStepMany(r_old,r_new,M,astar,omega,mm,eigenlm,R1old,R1m,2);
  }

  /* Outgoing mode */
  CKerr_OutgoingR3(rmatch,M,astar,omega,mm,eigenlm,R3m);
  aleph[0] = R3m[0]*R1m[2] - R3m[1]*R1m[3] - R1m[0]*R3m[2] + R1m[1]*R3m[3];
  aleph[1] = R3m[0]*R1m[3] + R3m[1]*R1m[2] - R1m[0]*R3m[3] - R1m[1]*R3m[2];
  aleph[0] /= rmatch*(rmatch-2*M)+a*a;
  aleph[1] /= rmatch*(rmatch-2*M)+a*a;

  /* Save c_{14} */
  cscat[2] =  aleph[1]/(2.*omega);
  cscat[3] = -aleph[0]/(2.*omega);

  /* Next we will get the W_{41}/Delta Wronskian.  Match at some very large radius. */
  rnorm = fabs(R1m[0]) + fabs(R1m[1]) + fabs(R1m[2]*rmatch) + fabs(R1m[3]*rmatch);
  for(j=0;j<4;j++) R1old[j] = R1m[j] / rnorm;
  r_old[0] = rmatch; r_old[1] = 0.;
  while (r_old[0] < CKERR_OUTER_BDY/fabs(omega) ) {
    r_new[0] = r_old[0] + CKERR_RMODE_STEP_EXTR*(r_old[0]>1./fabs(omega)? 1./fabs(omega): r_old[0]); r_new[1] = 0.;
    CKerr_RadialStepMany(r_old,r_new,M,astar,omega,mm,eigenlm,R1old,R1new,2);
    for(j=0;j<2;j++) r_old[j]=r_new[j];
    for(j=0;j<4;j++) R1old[j]=R1new[j];
  }
  CKerr_AsympR3(r_new,M,astar,omega,mm,eigenlm,R3asym);

  /* Obtain c_{13} by division */
  cscat[0] = (R1new[0]*R3asym[0]+R1new[1]*R3asym[1])/(R3asym[0]*R3asym[0]+R3asym[1]*R3asym[1]);
  cscat[1] = (R1new[1]*R3asym[0]-R1new[0]*R3asym[1])/(R3asym[0]*R3asym[0]+R3asym[1]*R3asym[1]);
  cscat[0] *= rnorm;
  cscat[1] *= rnorm;

  /* The determinant */
  beta[0] = -2*M*sqrt(1.-astar*astar);
  beta[1] = -a*mm + 2.*M*omega*rhp;
  detc[0] = -beta[1]/omega;
  detc[1] = beta[0]/omega;

  /* The lower elements: c24, then c23 */
  cscat[6] = -(detc[0]*cscat[0]+detc[1]*cscat[1])/alpha;
  cscat[7] = -(-detc[0]*cscat[1]+detc[1]*cscat[0])/alpha;
  cscat[4] = -pow(2.*omega,8)/C2/alpha * (detc[0]*cscat[2]+detc[1]*cscat[3]);
  cscat[5] = -pow(2.*omega,8)/C2/alpha * (-detc[0]*cscat[3]+detc[1]*cscat[2]);

  /* Now the inverse elements */
  idetc[0] = detc[0]/(detc[0]*detc[0]+detc[1]*detc[1]);
  idetc[1] = -detc[1]/(detc[0]*detc[0]+detc[1]*detc[1]);
  cscat[ 8] =   cscat[6]*idetc[0] - cscat[7]*idetc[1];
  cscat[ 9] =   cscat[6]*idetc[1] + cscat[7]*idetc[0];
  cscat[10] = -(cscat[2]*idetc[0] - cscat[3]*idetc[1]);
  cscat[11] = -(cscat[2]*idetc[1] + cscat[3]*idetc[0]);
  cscat[12] = -(cscat[4]*idetc[0] - cscat[5]*idetc[1]);
  cscat[13] = -(cscat[4]*idetc[1] + cscat[5]*idetc[0]);
  cscat[14] =   cscat[0]*idetc[0] - cscat[1]*idetc[1];
  cscat[15] =   cscat[0]*idetc[1] + cscat[1]*idetc[0];

  /* Scattering probabilities, if requested */
  if (aux!=NULL) {
    /* Wave coming in from past null infinity (R1 solution) */
    aux[2] = C2/pow(2.*omega,8) * (cscat[0]*cscat[0]+cscat[1]*cscat[1]) / (cscat[2]*cscat[2]+cscat[3]*cscat[3]);
    aux[3] = C2/pow(2.*omega,8) * alpha / (cscat[2]*cscat[2]+cscat[3]*cscat[3]);
  }

}

#undef CKERR_NSTEP_YSLM
#undef CKERR_MODE_OVERFLOW
#undef CKERR_NSTEP_TORTOISE
#undef CKERR_INNER_BDY
#undef CKERR_RMODE_STEP
#undef CKERR_RMODE_STEP_EXTR
#undef CKERR_OUTER_BDY
