#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define CKERR_NBISECT_ITER 56
#define CKERR_NBISECT_ITER2 32
#define CKERR_RESOLUTION_J_INTEG 64
#define CKERR_J_TOL 1e-9
#define CKERR_ACTION_MAX 1e49
#define CKERR_NPOINT_DERIV 3
#define CKERR_DACTION_DERIV 5e-4

/* Routine to convert (E,Q,L) --> (Jr,Jtheta,Jphi).
 * Returns 1 if successful, 0 if no orbit with those parameters (energy too low or Q<0), and
 * 2 if no orbit with these parameters (energy too high).
 * If the ancillary data is not null and there is a valid orbit, report it:
 *
 * ancillary[0] = inclination
 * ancillary[1] = pericenter
 * ancillary[2] = apocenter
 */
int CKerr_EQL2J(double *EQL, double *J, double M, double astar, double *ancillary) {
  long i;
  double zminus2, zplusinv2, dz, value, a;
  double zminus, zplusinv, zmzp;
  double Jt_integ, y, cosy, siny, arg, dy, secharg;
  double K;
  double c[5], disc, rinfl, rinfl2, rpeak, rdip, r, rh;
  double rturn[4];
  double rmid, drmax, Jr, ur;
  double denom;
  int inflsign = 0;
  int returnflag = 1;

  a = M*astar;

  /* Reject unbound orbits */
  if (EQL[0]<=-1) return(0);
  if (EQL[0]>=1) return(2);

  /* Azimuthal action is easy */
  J[2] = EQL[2];

  /* Vertical action: first get z-limits */
  if (EQL[1]<0) return(0);
  zminus = dz = 0.5;
  for(i=1;i<=CKERR_NBISECT_ITER;i++) {
    dz /= 2;
    zminus2 = zminus*zminus;
    value = EQL[1]*(1-zminus2) - a*a*(1.-EQL[0]*EQL[0])*zminus2*(1-zminus2) - EQL[2]*EQL[2]*zminus2;
    zminus += value>0? dz: -dz;
  }
  zplusinv = dz = 0.5;  
  for(i=1;i<=CKERR_NBISECT_ITER;i++) {
    dz /= 2;
    zplusinv2 = zplusinv*zplusinv;
    value = EQL[1]*(zplusinv2-1)*zplusinv2 - a*a*(1.-EQL[0]*EQL[0])*(zplusinv2-1) - EQL[2]*EQL[2]*zplusinv2;
    zplusinv += value>0? dz: -dz;
  }
  zmzp = zminus*zplusinv;
  /* Now do integral */
  Jt_integ = 0.;
  for(i=0;i<CKERR_RESOLUTION_J_INTEG*CKERR_NBISECT_ITER;i++) {
    arg = (i+0.5)/(double)CKERR_RESOLUTION_J_INTEG*log(2.);
    y = M_PI/2.*tanh(arg);
    secharg = 2.*exp(-arg)/(1.+exp(-2*arg));
    dy = M_PI/2.*secharg*secharg/(double)CKERR_RESOLUTION_J_INTEG*log(2.);
    cosy = cos(y); siny = sin(y);
    denom = 1.-zminus*zminus*siny*siny;
    if (denom>0) Jt_integ += dy*sqrt(1.-zmzp*siny*siny)*(1. - (1.-zminus*zminus)*siny*siny/denom);
#if 0
    Jt_integ += cosy*cosy*sqrt(1.-zmzp*siny*siny)/(1.-zminus*zminus*siny*siny)*dy;
#endif
  }
  Jt_integ /= M_PI/2.;
#ifdef CKERR_DIAG
  fprintf(stderr, "zminus=%12.5le zplusinv=%12.5le Jt_integ=%19.12le\n", zminus, zplusinv, Jt_integ);
#endif
  J[1] = Jt_integ * zminus*sqrt(EQL[1]);

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

#ifdef CKERR_DIAG
  fprintf(stderr, "r_infl=%12.5le,%12.5le r_dip,peak=%12.5le,%12.5le roots = %12.5le %12.5le %12.5le %12.5le\n",
    rinfl2, rinfl, rdip, rpeak, rturn[0], rturn[1], rturn[2], rturn[3]);
#endif

  /* Now we may finally do the action integral!  Substitution:
   * r = rmid + drmax*tanh(arg)
   * arg = -infty .. +infty, then *=2, and /=2pi to get Jr.
   */
  rmid = 0.5*(rturn[3]+rturn[2]);
  drmax = 0.5*(rturn[3]-rturn[2]);
  Jr = 0.;
  if (returnflag>=1 && EQL[0]<1)
    for(i=-CKERR_RESOLUTION_J_INTEG*CKERR_NBISECT_ITER; i<CKERR_RESOLUTION_J_INTEG*CKERR_NBISECT_ITER; i++) {
      arg = (i+0.5)/(double)CKERR_RESOLUTION_J_INTEG*log(2.);
      secharg = 2.*exp(-arg)/(1.+exp(-2*arg));
      r = rmid + drmax*tanh(arg);
      ur = (((c[4]*r+c[3])*r+c[2])*r+c[1])*r+c[0];
      ur = ur>0? sqrt(ur)/(r*(r-2.*M)+a*a): 0;
      Jr += secharg*secharg/(double)CKERR_RESOLUTION_J_INTEG*log(2.) * drmax * ur;
    }
  Jr /= M_PI;
  J[0] = Jr;

  if (ancillary!=NULL) {
    ancillary[0] = zminus<1? asin(zminus): M_PI/2.; if (EQL[2]<0) ancillary[0]=M_PI-ancillary[0]; /* inclination */
    ancillary[1] = rturn[2]; /* pericenter */
    ancillary[2] = rturn[3]; /* apocenter */
  }

  return(returnflag);
#ifdef CKERR_DIAG
#undef CKERR_DIAG
#endif
}

/* Determines the minimum energy for a given (Q,L) via a bisection search.
 */
double CKerr_Emin(double Q, double L, double M, double astar) {

  long i; double delta;  
  double EQL[3], J[3];

  EQL[1] = Q; EQL[2] = L;
  EQL[0] = 0.775;

  delta = 0.225;
  for(i=0; i<CKERR_NBISECT_ITER; i++) {
    delta /= 2.;
    EQL[0] += CKerr_EQL2J(EQL,J,M,astar,NULL)>=1? -delta: delta;
  }
  return(EQL[0]);
}

/* Determines the maximum energy for a given (Q,L) via a bisection search.
 * Returns the maximum radial action to Jrmax if desired (not NULL).
 */
double CKerr_Emax(double Q, double L, double M, double astar, double *Jrmax) {

  long i; double delta;  
  double EQL[3], J[3];

  EQL[1] = Q; EQL[2] = L;

  EQL[0] = CKerr_Emin(Q,L,M,astar);
  if (EQL[0]<0.95) {
    delta = (1.-EQL[0]) * 0.1;
    while(CKerr_EQL2J(EQL,J,M,astar,NULL)<2) {
      EQL[0] += delta;
    }
    EQL[0] -= 0.5*delta;
  }
  else {
    EQL[0] = 0.975; delta = 0.025;
  }

  /* EQL[0] = 0.775; delta = 0.225; */

  for(i=0; i<CKERR_NBISECT_ITER; i++) {
    delta /= 2.;
    EQL[0] += CKerr_EQL2J(EQL,J,M,astar,NULL)>=2? -delta: delta;
  }

  if (Jrmax!=NULL) {
    *Jrmax = EQL[0]<1? J[0]: CKERR_ACTION_MAX;
    EQL[0] -= CKERR_J_TOL;
    if (CKerr_EQL2J(EQL,J,M,astar,NULL)==0) *Jrmax=0;
    EQL[0] += CKERR_J_TOL;
  }
  return(EQL[0]);
}

/* Returns the energy for given (Q,L,Jr) via a bisection search.  Returns a number outside -1..+1
 * if the desired action does not exist.  If successful and requested (!=NULL), returns Jtheta.
 */
double CKerr_QLJr2E(double Q, double L, double Jr, double M, double astar, double *Jtheta) {
  
  long i,flag; double delta;
  double EQL[3], J[3];
  double Emin, Emax, Jrmax;

  Emin = CKerr_Emin(Q,L,M,astar);
  Emax = CKerr_Emax(Q,L,M,astar,&Jrmax);

  if (Jr<0) return(-2);
  if (Jr>Jrmax) return(2);

  EQL[1] = Q; EQL[2] = L;

  /* Bisection to get the energy */
  EQL[0] = (Emin+Emax)/2.;
  delta = (Emax-Emin)/2.;
  for(i=0; i<CKERR_NBISECT_ITER; i++) {
    delta /= 2.;
    flag=CKerr_EQL2J(EQL,J,M,astar,NULL);
    if (flag==0) J[0]=-1;
    if (flag==2) J[0]=2*Jr+M;
    EQL[0] += J[0]>Jr? -delta: delta;
  }
  flag=CKerr_EQL2J(EQL,J,M,astar,NULL);
  if (flag==2) return(2);

  /* Vertical action if requested */
  if (Jtheta!=NULL) *Jtheta=J[1];

  /* Return a failure if we are below zero action. */
  EQL[0] += CKERR_J_TOL;
  flag=CKerr_EQL2J(EQL,J,M,astar,NULL);
  if (flag==0) return(-2);
  EQL[0] -= CKERR_J_TOL;
  return(EQL[0]);
}

/* Takes in an action J and sets the energy, Carter constant, and angular momentum.
 * Returns 1 (successful) or 0 (failed: Jr too large).
 */
int CKerr_J2EQL(double *J, double *EQL, double M, double astar) {

  int i;
  int id=0,iu=0;
  double current_E, current_Jtheta, delta;

  /* Angular momentum is azimuthal action */
  EQL[2] = J[2];

  /* The rest of this is a 2D nonlinear equation solver.  We first guess the Carter constant,
   * and then iterate until the desired J_theta is obtained.
   */
  EQL[1] = J[1]*J[1]; /* Guess */

  /* Final bisection in negative powers of 2 */
  delta = sqrt(2.);
  for(i=0;i<CKERR_NBISECT_ITER;i++) {
    current_E = CKerr_QLJr2E(EQL[1],J[2],J[0],M,astar,&current_Jtheta);
    if (current_E<-1) current_Jtheta = -CKERR_ACTION_MAX;
    if (current_E>1) current_Jtheta = CKERR_ACTION_MAX;
    EQL[1] *= current_Jtheta<J[1]? delta: 1./delta;
    if (current_Jtheta<J[1]) {iu++;} else {id++;}
    if (iu*id>0) break;
  }
  for(i=0;i<CKERR_NBISECT_ITER;i++) {
    current_E = CKerr_QLJr2E(EQL[1],J[2],J[0],M,astar,&current_Jtheta);
    if (current_E<-1) current_Jtheta = -CKERR_ACTION_MAX;
    if (current_E>1) current_Jtheta = CKERR_ACTION_MAX;
#if 0
    if (fabs(current_E)>1) current_Jtheta = -CKERR_ACTION_MAX;
#endif
    EQL[1] *= current_Jtheta<J[1]? delta: 1./delta;
    delta=sqrt(delta);
  }

  EQL[0] = current_E;
  if (fabs(EQL[0])>1) return(0);
  return(1);
}

/* Computes the derivatives of the energy and Carter constant with respect to each
 * of the actions.  Returns 1 if successful, 0 if no such orbit.
 */
int CKerr_dEQdJ(double *J, double *EQL, double *dEdJ, double *dQdJ, double M, double astar) {

  int dirflag;
  long i,j;
  double dj;
  double Jcurrent[3], EQLcurrent[3];
  double newE[CKERR_NPOINT_DERIV], newQ[CKERR_NPOINT_DERIV], weights[CKERR_NPOINT_DERIV];

  /* Tell us if there's a failure. */
  if (CKerr_J2EQL(J,EQL,M,astar)==0) return(0);

  dj = CKERR_DACTION_DERIV * (J[0]+J[1]+fabs(J[2]));

  /* Derivatives with respect to phi-action */
  for(j=0;j<3;j++) Jcurrent[j] = J[j];
  newE[0] = EQL[0];
  newQ[0] = EQL[1];
  for(i=1;i<CKERR_NPOINT_DERIV;i++) {
    Jcurrent[2] = J[2]+i*dj;
    CKerr_J2EQL(Jcurrent,EQLcurrent,M,astar);
    newE[i] = EQLcurrent[0];
    newQ[i] = EQLcurrent[1];
  }

  /* Get numerical derivative of E & Q.
   * Uses the rule that the derivative of an N-1th order polynomial with roots at 0...N-1 (except i
   * where it evaluates to 1) at 0 is [(-1)^{N-2} (N-1)!/i] / [prod_{j=0..N-1,j!=i} (i-j)]
   */
  for(i=1;i<CKERR_NPOINT_DERIV;i++) {
    weights[i] = 1./i;
    for(j=0;j<CKERR_NPOINT_DERIV;j++) {
      weights[i] *= -1.;
      if (j>0) weights[i] *= j;
      if (j!=i) weights[i] /= i-j;
    }
  }
  dEdJ[2] = dQdJ[2] = 0.;
  for(i=1;i<CKERR_NPOINT_DERIV;i++) {
    dEdJ[2] += (newE[i] - newE[0])*weights[i]/dj;
    dQdJ[2] += (newQ[i] - newQ[0])*weights[i]/dj;
  }

  /* Derivatives with respect to theta-action */
  for(j=0;j<3;j++) Jcurrent[j] = J[j];
  newE[0] = EQL[0];
  newQ[0] = EQL[1];
  for(i=1;i<CKERR_NPOINT_DERIV;i++) {
    Jcurrent[1] = J[1]+i*dj;
    CKerr_J2EQL(Jcurrent,EQLcurrent,M,astar);
    newE[i] = EQLcurrent[0];
    newQ[i] = EQLcurrent[1];
  }
  dEdJ[1] = dQdJ[1] = 0.;
  for(i=1;i<CKERR_NPOINT_DERIV;i++) {
    dEdJ[1] += (newE[i] - newE[0])*weights[i]/dj;
    dQdJ[1] += (newQ[i] - newQ[0])*weights[i]/dj;
  }

  /* For the radial action, must set dj appropriately. */
  for(j=0;j<3;j++) Jcurrent[j] = J[j];
  Jcurrent[0] = 2*J[0];
  dirflag = CKerr_J2EQL(J,EQL,M,astar);
  if (dirflag==0) {
    if (dj*CKERR_NPOINT_DERIV>J[0]) dj = -J[0]/CKERR_NPOINT_DERIV;
  }
  dj /= 3.;
  /* Derivatives with respect to r-action */
  for(j=0;j<3;j++) Jcurrent[j] = J[j];
  newE[0] = EQL[0];
  newQ[0] = EQL[1];
  for(i=1;i<CKERR_NPOINT_DERIV;i++) {
    Jcurrent[0] = J[0]+i*dj;
    CKerr_J2EQL(Jcurrent,EQLcurrent,M,astar);
    newE[i] = EQLcurrent[0];
    newQ[i] = EQLcurrent[1];
  }
  dEdJ[0] = dQdJ[0] = 0.;
  for(i=1;i<CKERR_NPOINT_DERIV;i++) {
    dEdJ[0] += (newE[i] - newE[0])*weights[i]/dj;
    dQdJ[0] += (newQ[i] - newQ[0])*weights[i]/dj;
  }

  return(1);
}

/* Computes the M-inverse matrix and returns it in row-order (0,0) (0,1)
 * (0,2) ... (2,2).  Uses 9 slots.  This routine uses the numerical
 * differentiation method.
 * Returns 1 if successful, 0 if no such orbit.
 */
int CKerr_MinverseND(double *J, double *Minv, double M, double astar) {

  int i;
  double EQL[3], dEdJ[3], dQdJ[3], Mm[9], detM;

  if (CKerr_dEQdJ(J,EQL,dEdJ,dQdJ,M,astar)==0) return(0);

  /* Build the M-matrix */
  Mm[0] = dEdJ[0]; Mm[1] = dQdJ[0]; Mm[2] = 0;
  Mm[3] = dEdJ[1]; Mm[4] = dQdJ[1]; Mm[5] = 0;
  Mm[6] = dEdJ[2]; Mm[7] = dQdJ[2]; Mm[8] = 1;

  /* The determinant -- use last column */
  detM = Mm[0]*Mm[4]-Mm[1]*Mm[3];

  /* Inverses */
  Minv[0] = Mm[4]*Mm[8]-Mm[5]*Mm[7];
  Minv[3] = Mm[5]*Mm[6]-Mm[3]*Mm[8];
  Minv[6] = Mm[3]*Mm[7]-Mm[4]*Mm[6];
  Minv[1] = Mm[7]*Mm[2]-Mm[8]*Mm[1];
  Minv[4] = Mm[8]*Mm[0]-Mm[6]*Mm[2];
  Minv[7] = Mm[6]*Mm[1]-Mm[7]*Mm[0];
  Minv[2] = Mm[1]*Mm[5]-Mm[2]*Mm[4];
  Minv[5] = Mm[2]*Mm[3]-Mm[0]*Mm[5];
  Minv[8] = Mm[0]*Mm[4]-Mm[1]*Mm[3];
  for(i=0;i<9;i++) Minv[i]/=detM;

  return(1);
}

/* Computes the M-inverse matrix and returns it in row-order (0,0) (0,1)
 * (0,2) ... (2,2).  Uses 9 slots.  This routine uses analytic derivatives
 * inside the action integral and is recommended for precision work.
 * Returns 1 if successful, 0 if no such orbit.
 */
int CKerr_Minverse(double *J, double *Minv, double M, double astar) {

  int flag;
  long i;
  double EQL[3], newJ[3], ancillary[3];
  double zminus, a;
  double zminus__sqrtQ;
  double zplusinv, zplusinv2, dz, value, dz__u_z__1zz, arg, z, secharg;
  double c[5], sumroot, prodroot;
  double r, dr__u_r, Delta, factor;

  /* Get ~E,Q,~L and ancillary data */
  flag = CKerr_J2EQL(J,EQL,M,astar);
  if (flag==0) return(0);
  CKerr_EQL2J(EQL,newJ,M,astar,ancillary);
  a = M*astar;

#if 0
fprintf(stderr, "Minverse: EQL: %15.13lf %15.13lf %15.13lf J: %15.13lf %15.13lf %15.13lf\n", EQL[0], EQL[1], EQL[2], newJ[0], newJ[1], newJ[2]);
#endif

  /* Derivatives of Jphi */
  Minv[2] = 0.; Minv[5] = 0.; Minv[8] = 1.;

  /* Derivatives of Jtheta */
  zminus = sin(ancillary[0]);
  zplusinv = dz = 0.5;
  for(i=1;i<=CKERR_NBISECT_ITER;i++) {
    dz /= 2;
    zplusinv2 = zplusinv*zplusinv;
    value = EQL[1]*(zplusinv2-1)*zplusinv2 - a*a*(1.-EQL[0]*EQL[0])*(zplusinv2-1) - EQL[2]*EQL[2]*zplusinv2;
    zplusinv += value>0? dz: -dz;
  }
  zminus__sqrtQ = zminus>1e-8? zminus/sqrt(EQL[1]): 1./sqrt(EQL[2]*EQL[2]+a*a*(1.-EQL[0]*EQL[0])); /* fix Taylor expansion for low inclinations */

  /* Now do integral */
  Minv[1] = Minv[4] = Minv[7] = 0.;
  for(i=0;i<CKERR_RESOLUTION_J_INTEG*CKERR_NBISECT_ITER;i++) {
    arg = (i+0.5)/(double)CKERR_RESOLUTION_J_INTEG*log(2.);
    z = zminus * tanh(arg);
    secharg = 2.*exp(-arg)/(1.+exp(-2*arg));
    dz__u_z__1zz = zminus__sqrtQ/(double)CKERR_RESOLUTION_J_INTEG*log(2.) / sqrt(1.-z*z*zplusinv*zplusinv) * secharg;
    Minv[1] += dz__u_z__1zz * 2.*a*a*EQL[0]*z*z;
    Minv[4] += dz__u_z__1zz;
    Minv[7] += dz__u_z__1zz * (-2.)*EQL[2]*z*z/(1.-z*z);
  }
  Minv[1] /= M_PI; Minv[4] /= M_PI; Minv[7] /= M_PI;  

  /* Derivatives of Jr */
  /* First get additional roots of P(r) */
  c[0] = -a*a*EQL[1];                                             /* negative or 0 */
  c[1] = 2*M*(EQL[1] + (EQL[2]-a*EQL[0])*(EQL[2]-a*EQL[0]));      /* positive or 0 */
  c[2] = -a*a*(1.-EQL[0]*EQL[0]) - EQL[2]*EQL[2] - EQL[1];        /* negative      */
  c[3] = 2.*M;                                                    /* positive      */
  c[4] = EQL[0]*EQL[0]-1.;                                        /* negative      */
  sumroot = -c[3]/c[4] - ancillary[1] - ancillary[2];
  prodroot = c[0]/c[4] / ancillary[1] / ancillary[2];

  /* Now do integral */
  Minv[0] = Minv[3] = Minv[6] = 0.;
  factor = 1./sqrt(1.-EQL[0]*EQL[0]);
  for(i=-CKERR_RESOLUTION_J_INTEG*CKERR_NBISECT_ITER;i<CKERR_RESOLUTION_J_INTEG*CKERR_NBISECT_ITER;i++) {
    arg = (i+0.5)/(double)CKERR_RESOLUTION_J_INTEG*log(2.);
    r = (ancillary[2]+ancillary[1])/2. + (ancillary[2]-ancillary[1])/2.*tanh(arg);
    secharg = 2.*exp(-arg)/(1.+exp(-2*arg));
    Delta = r*(r-2.*M)+a*a;
    dr__u_r = Delta*secharg/(double)CKERR_RESOLUTION_J_INTEG*log(2.) / sqrt(r*(r-sumroot)+prodroot) * factor;
    Minv[0] += dr__u_r * ( 2.*((r*r+a*a)*EQL[0]-a*EQL[2])*(r*r+a*a)/Delta/Delta + 2*a*(EQL[2]-a*EQL[0])/Delta );
    Minv[3] += dr__u_r * ( -1./Delta );
    Minv[6] += dr__u_r * ( -2.*a*((r*r+a*a)*EQL[0]-a*EQL[2])/Delta/Delta - 2.*(EQL[2]-a*EQL[0])/Delta );
  }
  Minv[0] /= 2.*M_PI; Minv[3] /= 2.*M_PI; Minv[6] /= 2.*M_PI;

  return(flag);
}

/* Returns the phase-space coordinates (r,theta,phi,u_r,u_theta,u_phi) of the point
 * on the torus with actions J and at zero angles psi=(0,0,0).
 * Returns 1 if successful and 0 if the torus does not exist.
 */
int CKerr_TorusOrigin(double *J, double *xu, double M, double astar) {

  double EQL[3], newJ[3];
  double ancillary[3];

  if(CKerr_J2EQL(J,EQL,M,astar)==0) return(0);
  CKerr_EQL2J(EQL,newJ,M,astar,ancillary);

  /* Set the parameters */
  xu[0] = ancillary[1];      /* r */
  xu[1] = M_PI/2.;           /* theta */
  xu[2] = 0.;                /* phi */
  xu[3] = 0.;                /* u_r */
  xu[4] = sqrt(EQL[1]);      /* u_theta */
  xu[5] = J[2];              /* u_phi */

  return(1);
}

/* Takes in the M^{-1} matrix and returns the fundamental frequencies.
 * Sizes: Minv[9], Omega[3].
 */
void CKerr_Minv2Omega(double *Minv, double *Omega) {

  double detMinv;

  detMinv = Minv[0]*Minv[4]-Minv[1]*Minv[3];
  Omega[0] = Minv[4]/detMinv;
  Omega[1] = -Minv[3]/detMinv;
  Omega[2] = (Minv[3]*Minv[7]-Minv[4]*Minv[6])/detMinv;
}

/* Returns the flow vector for the Hamiltonian and Carter constant, I dH and I dQ,
 * in 6D phase space.
 */
void CKerr_FlowVectors(double *xu, double *xudot_H, double *xudot_Q, double M, double astar) {

  int i;
  double Sigma, Delta, a, gc_tt, gc_tphi, gc_phiphi;
  double u_t, uc_t, uc_r, uc_theta, uc_phi;
  double costheta, sintheta, Q, K;
  double QEcoef;

  /* Parameters */
  a = astar*M;
  costheta = cos(xu[1]);
  sintheta = sin(xu[1]);
  Delta = xu[0]*(xu[0]-2.*M)+a*a;
  Sigma = xu[0]*xu[0]+a*a*costheta*costheta;
  gc_tt = -( (xu[0]*xu[0]+a*a)*(xu[0]*xu[0]+a*a) - Delta*a*a*sintheta*sintheta )/(Delta*Sigma);
  gc_tphi = -2*a*M*xu[0]/Delta/Sigma;
  gc_phiphi = (Delta - a*a*sintheta*sintheta) / (Delta*Sigma*sintheta*sintheta);

  /* Hamiltonian trajectories.
   * Covariant velocities
   */
  u_t = gc_tphi*xu[5] - sqrt(gc_tphi*xu[5]*gc_tphi*xu[5] - gc_tt*(1. + gc_phiphi*xu[5]*xu[5] + Delta/Sigma*xu[3]*xu[3] + xu[4]*xu[4]/Sigma));
  u_t /= -gc_tt;

  /* Contravariant velocities */
  uc_t = gc_tt*u_t + gc_tphi*xu[5];
  uc_phi = gc_tphi*u_t + gc_phiphi*xu[5];
  uc_r = Delta*xu[3]/Sigma;
  uc_theta = xu[4]/Sigma;

  /* Carter constants */
  Q = xu[4]*xu[4] + a*a*(1.-u_t*u_t)*costheta*costheta + xu[5]*xu[5]*costheta*costheta/sintheta/sintheta;
  K = Q + (xu[5]+a*u_t)*(xu[5]+a*u_t);

  /* Position derivatives */
  xudot_H[0] = uc_r / uc_t;
  xudot_H[1] = uc_theta / uc_t;
  xudot_H[2] = uc_phi / uc_t;

  /* Momentum derivatives -- use conservation laws */
  xudot_H[3] = ( 2.*((xu[0]*xu[0]+a*a)*(-u_t) - a*xu[5])*xu[0]*(-u_t)/(Delta*Sigma)
               - xu[0]/Sigma
               - (xu[0]-M)*(xu[0]*xu[0]+K)/Delta/Sigma
               - 2.*(xu[0]-M)/Sigma*xu[3]*xu[3] ) / uc_t;
  xudot_H[4] = ( a*a*(1.-u_t*u_t)*costheta*sintheta + xu[5]*xu[5]*costheta/sintheta/sintheta/sintheta )/(Sigma*uc_t);
  xudot_H[5] = 0.;

  /* Derivatives of the Carter constant Q */
  QEcoef = -2.*a*a*(-u_t)*costheta*costheta;
  for(i=0;i<6;i++) xudot_Q[i] = QEcoef * xudot_H[i];
  xudot_Q[2] += 2.*xu[5]*costheta*costheta/sintheta/sintheta;
  xudot_Q[1] += 2.*xu[4];
  xudot_Q[4] += 2.*a*a*(1.-u_t*u_t)*costheta*sintheta + 2.*xu[5]*xu[5]*costheta/sintheta/sintheta/sintheta;
}

/* Returns the flow vector in the psi_r and psi_theta directions (I dJr and I dJtheta)
 * in 6D phase space.  Requires the M-inverse matrix Minv[9], which can be obtained using
 * CKerr_Minverse, in order to convert H and Q flow-vectors into angle flows.
 */
void CKerr_AngleFlowVectors(double *xu, double *xudot_r, double *xudot_theta, double *Minv, double M, double astar) {

  int i;
  double xudot_H[6], xudot_Q[6], xudot_L[6];

  /* Generate (H,Q,L) flow vectors */
  CKerr_FlowVectors(xu,xudot_H,xudot_Q,M,astar);
  for(i=0;i<6;i++) xudot_L[i] = 0.; xudot_L[2] = 1.;

  /* Angle flow vectors */
  for(i=0;i<6;i++) xudot_r[i] = Minv[0]*xudot_H[i] + Minv[3]*xudot_Q[i] + Minv[6]*xudot_L[i];
  for(i=0;i<6;i++) xudot_theta[i] = Minv[1]*xudot_H[i] + Minv[4]*xudot_Q[i] + Minv[7]*xudot_L[i];
}

/* Steps a point forward by an amount dpsi_r,dpsi_t using an N-step RK4 algorithm.
 */
void CKerr_AngleStep(double *xu_in, double *xu_out, double *Minv, double dpsi_r, double dpsi_t, double M, double astar, long N) {

  long i,ns;
  double xudr[6], xudt[6], h[6], xtemp[6];

  for(i=0;i<6;i++) xu_out[i] = xu_in[i];  /* Copy old position to new */
  dpsi_r/=N; dpsi_t/=N;                   /* Divide step size by N    */

  for(ns=0;ns<N;ns++) {
    /* 1st RK4 deriv */
    CKerr_AngleFlowVectors(xu_out,xudr,xudt,Minv,M,astar);
    for(i=0;i<6;i++) h[i] = (dpsi_r*xudr[i]+dpsi_t*xudt[i])/6.;

    /* 2nd RK4 deriv */
    for(i=0;i<6;i++) xtemp[i] = xu_out[i]+(dpsi_r*xudr[i]+dpsi_t*xudt[i])/2.;
    CKerr_AngleFlowVectors(xtemp,xudr,xudt,Minv,M,astar);
    for(i=0;i<6;i++) h[i] += (dpsi_r*xudr[i]+dpsi_t*xudt[i])/3.;

    /* 3rd RK4 deriv */
    for(i=0;i<6;i++) xtemp[i] = xu_out[i]+(dpsi_r*xudr[i]+dpsi_t*xudt[i])/2.;
    CKerr_AngleFlowVectors(xtemp,xudr,xudt,Minv,M,astar);
    for(i=0;i<6;i++) h[i] += (dpsi_r*xudr[i]+dpsi_t*xudt[i])/3.;

    /* 4th RK4 deriv */
    for(i=0;i<6;i++) xtemp[i] = xu_out[i]+(dpsi_r*xudr[i]+dpsi_t*xudt[i]);
    CKerr_AngleFlowVectors(xtemp,xudr,xudt,Minv,M,astar);
    for(i=0;i<6;i++) h[i] += (dpsi_r*xudr[i]+dpsi_t*xudt[i])/6.;

    /* Step */
    for(i=0;i<6;i++) xu_out[i] += h[i];
  }
}

/* Returns the angular momentum of a circular equatorial orbit with the specified
 * radius.  Formulas from Chandrasekhar ch.7.
 */
double CKerr_getL_CircEq(double M, double astar, double r) {
  double u, x, Q, a, E;

  a = M*astar;
  u = 1./r;
  Q = 1.-3*M*u-2.*a*u*sqrt(M*u);
  x = (-a*sqrt(u)+sqrt(M)) / sqrt(u*Q);
  E = (1.-2*M*u-a*u*sqrt(M*u))/sqrt(Q);

  return(a*E+x);
}

/* Returns dynamical data on circular equatorial orbit with the specified
 * radius.  Formulas from Chandrasekhar ch.7 and Paper I.
 * Returned data:
 *
 ***** Constants of the motion ****
 * info[0] = L = angular momentum
 * info[1] = E = energy
 ***** Advance rates ****
 * info[2] = w^t = conversion from proper to coordinate time
 * info[3] = Omega = orbital velocity seen at infinity
 ***** Epicyclic data ****
 * info[4] = kappa = epicyclic frequency
 * info[5] = Z = epicyclic impedance
 ***** Other differential data ***
 * info[6] = Omega' = d(Omega)/dr = angular velocity gradient
 *
 * The function returns 1 if successful and 0 otherwise (inside isco)
 */
int CKerr_getData_CircEq(double M, double astar, double r, double *info) {
  double u, x, Q, a, E, L, Delta;
  double gctt, gctphi, gcphiphi, wt, wphi, Omega;
  double C002, kappa, Z, deriv2;

  if (r<M+M*sqrt(1-astar*astar)) return(0);

  /* Basic data */
  a = M*astar;
  u = 1./r;
  Q = 1.-3*M*u+2.*a*u*sqrt(M*u); if (Q<=0) return(0);
  x = (-a*sqrt(u)+sqrt(M)) / sqrt(u*Q);
  E = info[1] = (1.-2*M*u+a*u*sqrt(M*u))/sqrt(Q);
  L = info[0] = a*E+x;

  /* Advance rates */
  Delta = r*(r-2*M)+a*a;
  gctt = -(r*r+a*a)*(r*r+a*a)/r/r/Delta + a*a/r/r;
  gctphi = -2*a*M/Delta/r;
  gcphiphi = (1.-2*M/r)/Delta;
  wt = info[2] = gctt*(-E) + gctphi*L;
  wphi = gctphi*(-E) + gcphiphi*L;
  Omega = info[3] = wphi/wt;
  if (gctt*E*E -2.*gctphi*E*L +gcphiphi*L*L>0) return(0);

  /* Epicycles */  
  /* kappa from radial equation */
  deriv2 = 4*E*E*(3*r*r+a*a) - 4*a*E*L - 12*r*(r-M) -2*(a*a+x*x);
  if (deriv2>=0) return(0);
  kappa = info[4] = sqrt(-deriv2/2.) / (r*r*wt);
  C002 = Delta/r/r/wt;
  Z = info[5] = kappa/C002;

  /* Angular velocity gradient -- derivative of Chandrasekhar 7.121 */
  info[6] = -1.5*Omega/r/(1.+a*Omega)/(1.+a*Omega);
  return(1);
}

/* Returns the location of the order-mm Lindblad resonance around the
 * specified black hole.  positive mm finds ILRs; negative mm finds the
 * OLR with order |mm|.
 *
 * Only integer mm gives Lindblad resonances, but we preserve the option to
 * insert fractional mm to search for higher-order eccentricity resonances --
 * e.g. the 3:1 resonance involving two factors of the test particle eccentricity
 * can be found by setting mm = 1.5.
 *
 * However may be ill-behaved if |mm|<1.
 */
double CKerr_FindLindbladResonance(double M, double astar, double r0, double mm) {

  int i, Dsign, msign, count1, count2;
  double Omega_s, r1, info[7], jumpfactor;

  /* Basic data: secondary orbital rate, sign of resonance */
  CKerr_getData_CircEq(M,astar,r0,info);
  Omega_s = info[3];
  msign = mm>0? 1: -1;

  /* First stage of bisection search -- find right octave.
   */
  r1 = r0;
  jumpfactor = 2.;
  count1 = count2 = 0;
  while (count1*count2==0) {
    if (CKerr_getData_CircEq(M,astar,r1,info)==0) {
      Dsign = -1;
    } else {
      Dsign = fabs(mm)*(info[3]-Omega_s)-msign*info[4]>0? -1: 1;
    }
    if (Dsign==1) {
      r1 /= jumpfactor; count1++;
    } else {
      r1 *= jumpfactor; count2++;
    }
  }

  /* Now zero in on the solution */
  for(i=0; i<CKERR_NBISECT_ITER; i++) {
    if (CKerr_getData_CircEq(M,astar,r1,info)==0) {
      Dsign = -1;
    } else {
      Dsign = fabs(mm)*(info[3]-Omega_s)-msign*info[4]>0? -1: 1;
    }
    r1 *= Dsign<0? jumpfactor: 1./jumpfactor;
    jumpfactor = sqrt(jumpfactor);
  }
  return(r1);
}

/* Obtains constants (E,Q,L) for a trajectory of specified Carter inclination,
 * pericenter, and apocenter.  Returns 1 if successful, 0 if no solution.
 */
int CKerr_FindEQL_IRR(double Ic, double rp, double ra, double *EQL, double M, double astar) {

  double L__F, Q__F2;
  double cp[3], ca[3], Delta, r2a2;
  double cd[3], fe, disc;
  double a, E, F;
  int i;

  a = M*astar;

  /* The basic program is to find energy E and F=sqrt(L^2+Q).  The turning
   * points give conic section restrictions in the (E,F) plane, and we will find
   * their intersection.  First we need the mapping of F-->L and F^2-->Q.
   */
  L__F = cos(Ic);
  Q__F2 = sin(Ic)*sin(Ic);

  /* Pericenter hyperbola.  Format: cp[0]*E^2 + 2*cp[1]*E*F + cp[2]*F^2 = 1. */
  Delta = rp*rp+a*a-2*rp*M;
  r2a2 = rp*rp+a*a;
  cp[0] = r2a2*r2a2 - Delta*a*a;
  cp[1] = -a*r2a2*L__F + a*Delta*L__F;
  cp[2] = a*a*L__F*L__F - Delta;
  for(i=0;i<=2;i++) cp[i] /= Delta*rp*rp;

  /* Apocenter hyperbola.  Format: ca[0]*E^2 + 2*ca[1]*E*F + ca[2]*F^2 = 1. */
  Delta = ra*ra+a*a-2*ra*M;
  r2a2 = ra*ra+a*a;
  ca[0] = r2a2*r2a2 - Delta*a*a;
  ca[1] = -a*r2a2*L__F + a*Delta*L__F;
  ca[2] = a*a*L__F*L__F - Delta;
  for(i=0;i<=2;i++) ca[i] /= Delta*ra*ra;

  /* Difference equations to get ratio of F/E */
  for(i=0;i<=2;i++) cd[i] = ca[i]-cp[i];
  disc = cd[1]*cd[1]-cd[0]*cd[2];
  if (disc<0) return(0);
#if 0
  fe = -cd[1]/cd[2] + sqrt(disc)/fabs(cd[2]);
  fprintf(stderr, "cd[] = %12.5le,%12.5le,%12.5le\n", cd[0],cd[1],cd[2]);
#endif
  fe = -cd[0]/(cd[1]+sqrt(disc));

  /* Now the pericenter hyperbola equation is (cp[0]+2*cp[1]*fe+cp[2]*fe*fe)*E^2 = 1 */
  if (cp[0]+2*cp[1]*fe+cp[2]*fe*fe<0) return(0);
  E = 1./sqrt(cp[0]+2*cp[1]*fe+cp[2]*fe*fe);
  F = fe*E;
  EQL[0] = E;
  EQL[1] = Q__F2*F*F;
  EQL[2] = L__F*F;
  return(1);
}

/* Obtains constants (E,Q,L) for a trajectory of specified Carter inclination,
 * and a circular orbit.  Returns 1 if successful, 0 if no solution.
 */
int CKerr_FindEQL_IRCirc(double Ic, double r, double *EQL, double M, double astar) {

  double L__F, Q__F2;
  double cp[3], ca[3], Delta, r2a2;
  double cd[3], fe, disc;
  double a, E, F;
  int i;

  a = M*astar;

  /* The basic program is to find energy E and F=sqrt(L^2+Q).  The turning
   * points give conic section restrictions in the (E,F) plane, and we will find
   * their intersection.  First we need the mapping of F-->L and F^2-->Q.
   */
  L__F = cos(Ic);
  Q__F2 = sin(Ic)*sin(Ic);

  /* Pericenter hyperbola.  Format: cp[0]*E^2 + 2*cp[1]*E*F + cp[2]*F^2 = 1. */
  Delta = r*r+a*a-2*r*M;
  r2a2 = r*r+a*a;
  cp[0] = r2a2*r2a2 - Delta*a*a;
  cp[1] = -a*r2a2*L__F + a*Delta*L__F;
  cp[2] = a*a*L__F*L__F - Delta;
  for(i=0;i<=2;i++) cp[i] /= Delta*r*r;

  /* Need to get derivative of c[0]*E^2 + 2*c[1]*E*F + c[2]*F^2 = c[3] with respect to radius. */
  ca[0] = 2*r2a2*(2*r) - 2*(r-M)*a*a;
  ca[1] = -2*a*r*L__F + a*2*(r-M)*L__F;
  ca[2] = -2*(r-M);
  for(i=0;i<=2;i++) ca[i] /= 2*(r-M)*r*r + 2*r*Delta;

  /* Difference equations to get ratio of F/E */
  for(i=0;i<=2;i++) cd[i] = ca[i]-cp[i];
  disc = cd[1]*cd[1]-cd[0]*cd[2];
  if (disc<0) return(0);
  fe = -cd[0]/(cd[1]+sqrt(disc));

  /* Now the pericenter hyperbola equation is (cp[0]+2*cp[1]*fe+cp[2]*fe*fe)*E^2 = 1 */
  if (cp[0]+2*cp[1]*fe+cp[2]*fe*fe<0) return(0);
  E = 1./sqrt(cp[0]+2*cp[1]*fe+cp[2]*fe*fe);
  F = fe*E;
  EQL[0] = E;
  EQL[1] = Q__F2*F*F;
  EQL[2] = L__F*F;
  return(1);
}

/* Finds the orbit at which the ratio of azimuthal to vertical frequencies is fratio,
 * for (Carter) inclination Ic.  Also returns the actions J[3].
 * Currently requires fratio >= 2, astar > 0.
 * Returns 1 if the orbit exists, 0 otherwise.
 */
int CKerr_FindResCirc(int fratio, double Ic, double *J, double M, double astar) {

  double r, dlnr, Minv[9], Omega[3], EQL[3], ancillary[3], EQL2[3], c[5];
  int i, f, flag, flag2, fl;
  double a = M*astar;

  /* Bisection search for the correct radius */
  r = 2*M;
  dlnr = 0.5*log(2.);
  for(i=0; i<CKERR_NBISECT_ITER2; i++) {

    /* flag == 0 if either (1) inside ISCO or (2) frequency ratio above fratio:1 */
    flag = flag2 = CKerr_FindEQL_IRCirc(Ic,r,EQL,M,astar);

    if (EQL[0]>=1 || EQL[0]<0.5) flag=0;
    if (flag==1) {
      
      c[0] = -a*a*EQL[1];                                             /* negative or 0 */
      c[1] = 2*M*(EQL[1] + (EQL[2]-a*EQL[0])*(EQL[2]-a*EQL[0]));      /* positive or 0 */
      c[2] = -a*a*(1.-EQL[0]*EQL[0]) - EQL[2]*EQL[2] - EQL[1];        /* negative      */
      c[3] = 2.*M;                                                    /* positive      */
      c[4] = EQL[0]*EQL[0]-1.;                                        /* negative      */

      if (2*c[2] + 6*c[3]*r + 12*c[4]*r*r >= 0) flag=2;
    }

    if (flag==1) {
      EQL[0] += CKERR_J_TOL;
      f = CKerr_EQL2J(EQL,J,M,astar,ancillary);
      if (f==1) {
        fl = CKerr_Minverse(J,Minv,M,astar);
        CKerr_Minv2Omega(Minv,Omega);
        if (Omega[2]/Omega[1]>fratio) flag=0;
      } else {
        flag = 0;
      }
    }

#if 0
fprintf(stderr, "r=%15.13lf (%15.13lf,%15.13lf) flag=%1d,%1d,%1d ratio=%15.13lf\n", r, ancillary[1], ancillary[2], flag2, f, fl, Omega[2]/Omega[1]);
fprintf(stderr, "     J[] = %15.12lf,%15.12lf,%15.12lf Omega[] = %15.12lf,%15.12lf,%15.12lf\n", J[0], J[1], J[2], Omega[0], Omega[1], Omega[2]);
fprintf(stderr, "     E = %15.13lf range = %15.13lf,%15.13lf\n", EQL[0], CKerr_Emin(EQL[1],EQL[2],1,astar), CKerr_Emax(EQL[1],EQL[2],1,astar,NULL));
    fprintf(stderr, "    %15.12lf %15.12lf %15.12lf\n", EQL[0], EQL[1], EQL[2]);
    fprintf(stderr, "<%1d>", CKerr_J2EQL(J,EQL2,M,astar));
    fprintf(stderr, " %15.12lf %15.12lf %15.12lf\n", EQL2[0], EQL2[1], EQL2[2]);
#endif
    r *= exp(flag==1? -dlnr: dlnr);
    dlnr/=2.;
  }

  flag = CKerr_FindEQL_IRCirc(Ic,r,EQL,M,astar);

fprintf(stderr, "%10.8lf %10.8lf %1d %10.8lf\n", Ic, r, flag, Omega[2]/Omega[1]);

  return(flag);
}

/* Finds the ISCO radius for a general inclination orbit, as a function of (Carter)
 * inclination Ic.
 */
double CKerr_FindISCO(double Ic, double M, double astar) {

  double a, r, EQL[3], dlnr, c[5];
  int i, flag;

  a = M*astar;

  /* Bisection search for the correct radius */
  r = 3*M;
  dlnr = 0.5*log(3.);
  for(i=0; i<CKERR_NBISECT_ITER; i++) {
    flag = CKerr_FindEQL_IRCirc(Ic,r,EQL,M,astar);
    if (EQL[0]>=1 || EQL[0]<0.5) flag=0;
    if (flag==1) {

      c[0] = -a*a*EQL[1];                                             /* negative or 0 */
      c[1] = 2*M*(EQL[1] + (EQL[2]-a*EQL[0])*(EQL[2]-a*EQL[0]));      /* positive or 0 */
      c[2] = -a*a*(1.-EQL[0]*EQL[0]) - EQL[2]*EQL[2] - EQL[1];        /* negative      */
      c[3] = 2.*M;                                                    /* positive      */
      c[4] = EQL[0]*EQL[0]-1.;                                        /* negative      */

      if (2*c[2] + 6*c[3]*r + 12*c[4]*r*r >= 0) flag=2;
    }

    r *= exp(flag==1? -dlnr: dlnr);
    dlnr/=2.;
  }
  return(r);
}

#undef CKERR_NBISECT_ITER
#undef CKERR_NBISECT_ITER2
#undef CKERR_RESOLUTION_J_INTEG
#undef CKERR_J_TOL
#undef CKERR_ACTION_MAX
#undef CKERR_NPOINT_DERIV
#undef CKERR_DACTION_DERIV
