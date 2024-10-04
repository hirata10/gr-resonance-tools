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


  fprintf(stderr, "r_infl=%12.5le,%12.5le r_dip,peak=%12.5le,%12.5le roots = %12.5le %12.5le %12.5le %12.5le\n",
    rinfl2, rinfl, rdip, rpeak, rturn[0], rturn[1], rturn[2], rturn[3]);


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

int radial_roots(double *EQL, double *J, double M, double astar, double *ancillary) {
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

  fprintf(stderr, "r_infl = %12.5le,%12.5le r_dip,peak = %12.5le,%12.5le inner and outer horizon = %12.5le %12.5le roots = %12.5le %12.5le %12.5le %12.5le\n",
    rinfl2, rinfl, rdip, rpeak, r_h[0], r_h[1], rturn[0], rturn[1], rturn[2], rturn[3]);
  


}

int main(int argc, char **argv) {
  double EQL[3]; //Command line input for EQL
  double J[3];
  double ancillary[3];
  double M = 1.0;
  double a = 0.5;

  //Set values to 0.9892 12.309 4.105 on command line
  //Set values to 0.9891 24.3007 2.4172 on command line

  sscanf(argv[1], "%lg", &EQL[0]);
  sscanf(argv[1], "%lg", &EQL[1]);
  sscanf(argv[1], "%lg", &EQL[2]);

  radial_roots(EQL, J, M, a, ancillary);

  CKerr_EQL2J(EQL, J, M, a, ancillary);

  printf("ancillary: %12.5le %12.5le %12.5le \n", ancillary[0], ancillary[1], ancillary[2]);
  printf("J's: %12.5le %12.5le %12.5le \n", J[0], J[1], J[2]);

  return(0);
  
}