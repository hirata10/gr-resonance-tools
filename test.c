#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "CKerr.h"

void action_test(void) {
  int i, j, flag;
  double EQL[3], J[]={0,0,0}, ancillary[3], M=1, astar=0.998;99046064;
  double newEQL[3], dEdJ[3], dQdJ[3];
  double djr=1e-4;
  double Jrmax;
  double Minv[9], xu[6], Omega[3], dxu_r[6], dxu_t[6];

  scanf("%lg %lg %lg", EQL, EQL+1, EQL+2);
  flag = CKerr_EQL2J(EQL,J,M,astar,ancillary);
  printf("[%1d] %12.5le %12.5le %12.5le | %12.5le %12.5le %12.5le [%12.10lf %12.10lf] ",
    flag, J[0], J[1], J[2], ancillary[0], ancillary[1], ancillary[2],
    CKerr_Emin(EQL[1], EQL[2], M, astar),
    CKerr_Emax(EQL[1], EQL[2], M, astar, &Jrmax));

  printf("%12.5le\n> ", Jrmax);

  for(djr=1e-6; djr<101; djr*=10) printf("%12.5le ", CKerr_QLJr2E(EQL[1], EQL[2], djr, M, astar, NULL));
  printf("\n");

  flag = CKerr_J2EQL(J, newEQL, M, astar);
  printf("inverse: %18.11le %18.11le %18.11le : %18.11le %18.11le %18.11le (%1d)\n",
    J[0], J[1], J[2], newEQL[0], newEQL[1], newEQL[2], flag);

#if 0
  CKerr_dEQdJ(J, newEQL, dEdJ, dQdJ, M, astar);
  printf("dEdJ = [%18.11le %18.11le %18.11le] dQdJ = [%18.11le %18.11le %18.11le]\n",
    dEdJ[0], dEdJ[1], dEdJ[2], dQdJ[0], dQdJ[1], dQdJ[2]);
  CKerr_MinverseND(J,Minv,M,astar);
  for(i=0;i<9;i++) {
    printf(" %18.11le", Minv[i]);
    if (i%3==2) printf("\n");
  }
#endif
  CKerr_Minverse(J,Minv,M,astar);
  for(i=0;i<9;i++) {
    printf(" %18.11le", Minv[i]);
    if (i%3==2) printf("\n");
  }
  CKerr_Minv2Omega(Minv,Omega);
  printf("Omega = %18.11le, %18.11le, %18.11le\n", Omega[0], Omega[1], Omega[2]);

  CKerr_TorusOrigin(J,xu,M,astar);

#if 0
  for(j=0;j<=32;j++) {
    for(i=0;i<6;i++) printf(" %17.13lf", xu[i]); printf("\n");
    if (j<32) CKerr_AngleStep(xu, xu, Minv, M_PI/16., 0., M, astar, 50);
//    CKerr_AngleFlowVectors(xu, dxu_r, dxu_t, Minv, M, astar);
  }
  for(j=1;j<=32;j++) {
    CKerr_AngleStep(xu, xu, Minv, 0., M_PI/16., M, astar, 50);
    for(i=0;i<6;i++) printf(" %17.13lf", xu[i]); printf("\n");
  }
#endif
  return;
}

void mode_test(void) {
  long l,mm,i,j;
  double chi;
  double *Y, *dYdt, theta;
  long nb, nu=10;
  double E[25], *b, S[25], dS[25];

#if 0
  scanf("%ld %lg", &l, &theta);
  Y = (double*)malloc((size_t)((2*l+1)*sizeof(double)));
  dYdt = (double*)malloc((size_t)((2*l+1)*sizeof(double)));
  CKerr_Yslm(theta, l, 0, Y, dYdt);
  for(i=0;i<2*l+1;i++) printf(" %13.9lf", Y[i]);
  printf("\n");
  for(i=0;i<2*l+1;i++) printf(" %13.9lf", dYdt[i]);
  printf("\n");
  printf("\n");
  CKerr_Yslm(theta, l, 2, Y, dYdt);
  for(i=0;i<2*l+1;i++) printf(" %13.9lf", Y[i]);
  printf("\n");
  for(i=0;i<2*l+1;i++) printf(" %13.9lf", dYdt[i]);
  printf("\n");
  free((char*)Y);
  free((char*)dYdt);
#endif

  scanf("%ld %lg",&mm,&chi);
  nb=CKerr_SbCoefN(-2,mm,nu,E,&b,chi,1e-11);
#if 1
  printf("(%2ld,%2ld)\n", nu, nb);
  for(i=0;i<nu;i++) {
    printf("%16.11lf :", E[i]);
    for(j=0;j<nb;j++) printf(" %9.6lf", b[i*nb+j]);
    printf("\n");
  }
  printf("\n");
#endif

#if 1
  for(theta=0; theta<M_PI+1e-3; theta+=M_PI/6.) {
    CKerr_SpheroidalYSLM(theta, -2, mm, nu, nb, b, S, dS);
    printf("%8.6lf", theta);
    for(i=0;i<nu;i++) {
      printf(" %9.6lf", S[i]);
    }
    printf("\n");
  }
#endif
  free((char*)b);
}

void coord_test(void) {
  double M = 1.;
  double r, rstar, astar, Delta, ReV, ImV;
  double R1[4], R3[4], rr[2], R1new[4], R1comp[4], h[2];

  scanf("%lg %lg", &r, &astar);
  rstar = CKerr_r2rstar(r,M,astar);
  printf("%19.12le, %19.12le, %19.12le\n", r, rstar, M+M*sqrt(1.-astar*astar));
  r = CKerr_rstar2r(rstar,M,astar,&Delta);
  printf("Delta = %19.12le\n", Delta);

  for(rstar = -100.; rstar<160.1; rstar+=10.0) {
    r = CKerr_rstar2r(rstar,M,astar,&Delta);
    rr[0] = r; rr[1] = 0;r>2? 4.*log(0.5*r): 0;
//    CKerr_IngoingR1(r,M,astar,0.5,2,6,R1);
    CKerr_OutgoingR3(r,M,astar,0.5,2,6,R3);

    h[0] = 0.01; h[1] = 0.02;
    CKerr_AsympR3(rr,M,astar,0.5,2,6,R1);
    rr[0] += h[0]; rr[1] += h[1];
    CKerr_AsympR3(rr,M,astar,0.5,2,6,R1comp);
    rr[0] -= h[0]; rr[1] -= h[1];
    CKerr_RadialRKStep(rr,h,M,astar,0.5,2,6,R1,R1new);

//    printf("%19.12le %19.12le %19.12le %19.12le %19.12le %19.12le\n", rr[0], rr[1], R1[0], R1[1], R1[2], R1[3]);
//    printf("                                        %19.12le %19.12le %19.12le %19.12le\n", R3[0], R3[1], R3[2], R3[3]);
//    printf("                                        %19.12le %19.12le %19.12le %19.12le\n", R1comp[0], R1comp[1], R1comp[2], R1comp[3]);
//    printf("                                        %19.12le %19.12le %19.12le %19.12le\n", R1new[0], R1new[1], R1new[2], R1new[3]);

    printf("%19.12le %19.12le %19.12le %12.9lf %12.9lf %12.9lf\n", rr[0], rr[0]*rr[0]-2*M*rr[0]+M*M*astar*astar,
      sqrt(R3[0]*R3[0]+R3[1]*R3[1]), atan2(R3[1],R3[0]),
      (R3[2]*R3[0]+R3[3]*R3[1])/(R3[0]*R3[0]+R3[1]*R3[1])*(1.-2*M*rr[0]/(rr[0]*rr[0]+M*M*astar*astar)),
      (R3[3]*R3[0]-R3[2]*R3[1])/(R3[0]*R3[0]+R3[1]*R3[1])*(1.-2*M*rr[0]/(rr[0]*rr[0]+M*M*astar*astar)));

#if 0
 sqrt(R3[2]*R3[2]+R3[3]*R3[3])*(1.-2*M*rr[0]/(rr[0]*rr[0]+M*M*astar*astar)), atan2(R3[3],R3[2]));
#endif

//    printf("%19.12le %19.12le %19.12le %12.9lf %19.12le %12.9lf\n", rr[0], rr[1], sqrt(R1[0]*R1[0]+R1[1]*R1[1]), atan2(R1[1],R1[0]), sqrt(R1[2]*R1[2]+R1[3]*R1[3]), atan2(R1[3],R1[2]));
//    printf("%12.6lf %16.13lf %19.12le %12.9lf %19.12le %12.9lf\n", rstar, r, sqrt(R1[0]*R1[0]+R1[1]*R1[1]), atan2(R1[1],R1[0]), sqrt(R1[2]*R1[2]+R1[3]*R1[3]), atan2(R1[3],R1[2]));
#if 0
    CKerr_Vr(rstar,M,astar,0.5,2,6,&ReV,&ImV);
    printf("%12.6lf %16.13lf %16.8lf %16.8lf\n", rstar, r, ReV, ImV);
#endif
  }

}

void gwscat_test(void) {

  double M=1, astar=0.9999;
  long mm, ll, il, nb;
  double omega, E[100], *b;
  double c[22];

  scanf("%ld %ld", &ll, &mm);

  il = ll-(abs(mm)>2? abs(mm): 2);
  for(omega=0.8; omega<1.005; omega+=0.01) {
    nb = CKerr_SbCoefN(-2,mm,il+1,E,&b,astar*M*omega,1e-9);
    free((char*)b);
    CKerr_GWScatMatrix(M,astar,omega,mm,E[il],c,c+16);
    printf("%17.10lE %17.10lE %17.10lE %17.10lE\n", omega, E[il], c[18], c[19]);
  }
}

void gwem_test(void) {

  long il,nl;
  double M=1, astar=0.9;
  double J[3], omega, Coef[625], E[25];
  double Minv[9], xuorig[6];
  long qr=0, qt=0, mm, qrmin, qrmax;
  long lmax=8;
  double QI,QH,EI,EH,LI,LH;
  double c[22];

  scanf("%lg %lg %lg", J, J+1, J+2);

  CKerr_Minverse(J,Minv,M,astar);
  CKerr_TorusOrigin(J,xuorig,M,astar);

EI=EH=LI=LH=0;

#if 1
for(mm=0; mm<=lmax; mm++) for(qt=-lmax; qt<=lmax; qt++) if (mm!=0 || qt>=0) {
qrmax = lmax;
qrmin = mm==0 && qt==0? 1: -qrmax;
for(qr=qrmin; qr<=qrmax; qr++) {
#else
for(mm=0; mm<=lmax; mm++) for(qt=-lmax+abs(mm); qt<=lmax-abs(mm); qt++) if (mm!=0 || qt>=0) {
qrmax = lmax-abs(mm)-abs(qt);
qrmin = mm==0 && qt==0? 1: -qrmax;
for(qr=qrmin; qr<=qrmax; qr++) {
#endif

  nl=lmax-mm+1; if (mm<=2) nl=lmax-1;

  CKerr_RadialFunc(Minv, xuorig, M, astar, qr, qt, mm, 32, 32, nl, Coef, &omega, E);
// double *J, double M, double astar, long qr, long qt, long mm, long Nr, long Nt, long nl, double *Coef, double *omega

  for(il=0;il<nl;il++) {
#if 0

    CKerr_GWScatMatrix(M,astar,omega,mm,E[il],c,c+16);

    printf("%9.5lf %12.5le (%12.5le %12.5le): ", E[il], CKerr_HorizonFluxNorm(M,astar,omega,mm,E[il]), c[16], c[17]);
    printf("%12.5le %12.5le %12.5le %12.5le out = %17.10le\n", Coef[4*il], Coef[4*il+1], Coef[4*il+2], Coef[4*il+3],
      (Coef[4*il+2]*Coef[4*il+2]+Coef[4*il+3]*Coef[4*il+3])/(omega*omega));
    printf("  (%11.4le %11.4le, %11.4le %11.4le)   (%11.4le %11.4le, %11.4le %11.4le)\n",
      c[0], c[1], c[2], c[3], c[8], c[9], c[10], c[11]);
    printf("  (%11.4le %11.4le, %11.4le %11.4le)   (%11.4le %11.4le, %11.4le %11.4le)\n",
      c[4], c[5], c[6], c[7], c[12], c[13], c[14], c[15]);
    printf(" scat coefs: from past infinity: %16.9lE %16.9lE\n", c[18], c[19]);
#endif
  }

  CKerr_QuantaEmitted(M,astar,qr,qt,mm,nl,E,Coef,omega,&QI,&QH);
#if 0
  printf("omega = %12.5le QI = %12.5le QH = %12.5le\n", omega, QI, QH);
#endif
  printf("%3ld %3ld %3ld %12.9lf %16.9le %16.9le %16.9le %16.9le\n", mm, qt, qr, omega, 2*QI*omega, 2*QH*omega, 2*QI*mm, 2*QH*mm);

  EI += 2*QI*omega;
  EH += 2*QH*omega;
  LI += 2*QI*mm;
  LH += 2*QH*mm;
}}

  printf("Edot %16.9le %16.9le %16.9le Ldot %16.9le %16.9le %16.9le\n", EI, EH, EI+EH, LI, LH, LI+LH);

  return;
}

void resonance_test(void) {
  double M=1, astar;
  double r, rILR, rOLR; long mm, lmin, lmax, ii;
  double Dprime, D1, D2;
  double info[7], infoI[7], infoO[7], S_i, S_o, r1_i, r1_o;
  double rrtable[] = {50, 30, 20, 12, 8, 6, 4};
  double r_init;

  scanf("%ld %lg %lg", &mm, &astar, &r_init);
  lmin = 2;
  if (abs(mm)>2) lmin=abs(mm);

  /* information: L,E,w^t,Omega,kappa,Z,Omega' */
    printf("         ");
    printf("--- SECONDARY --- ");
    printf         ("         ");
    printf("                  ");
    printf("                  ");
    printf("                  ");
    printf("--- RESONANCE --- \n");
    printf("        r         ");
    printf("      Omega       ");
    printf("        r         ");
    printf("      Omega       ");
    printf("      kappa       ");
    printf("        Z         ");
    printf("     dt/dtau      ");
    printf("        D'        ");
    printf("        S         ");
    printf("      astar       \n");
  r = 20;
#if 0
  for(r=r_init; r<250.01; r+=(r<4.9? 0.5: r<9.9? 1: r<19.1? 2: r<49.9? 10: r<99.9? 25: 50)) {
#else
  for(astar=-0.95; astar<=0.96; astar+=0.05) {
#endif
    if (CKerr_getData_CircEq(M,astar,r,info)) {
      rILR = CKerr_FindLindbladResonance(M,astar,r,mm);
      CKerr_getData_CircEq(M,astar,rILR,infoI);

      lmax = r>40.01? 20: r>20.01? 30: r>7.01? 40: 50;
      /* old jr = 3e-7*sqrt(r) */
      S_i = CKerr_LindbladResonanceStrength(M,astar,r,mm,lmax-lmin+1,4e-6*(r-rILR)*(r-rILR)*infoI[5],&r1_i);

      CKerr_getData_CircEq(M,astar,rILR+0.001*M,infoO);
      D1 = mm*infoO[3] - infoO[4];
      CKerr_getData_CircEq(M,astar,rILR-0.001*M,infoO);
      D2 = mm*infoO[3] - infoO[4];
      Dprime = (D1-D2)/(0.002*M) * (mm>0? 1: -1);

      printf("%17.10lE %17.10lE %17.10lE %17.10lE %17.10lE %17.10lE %17.10lE %17.10lE %17.10lE %17.10lE\n",
        r, info[3], rILR, infoI[3], infoI[4], infoI[5], infoI[2], Dprime, S_i, astar);
    }
  }

#if 0
  /* information: L,E,w^t,Omega,kappa,Z,Omega' */
    printf("        r         ");
    printf("        L         ");
    printf("        E         ");
    printf("       w^t        ");
    printf("      Omega       ");
    printf("      kappa       ");
    printf("        Z         ");
    printf("    dOmega/dr    \n");
  for(r=6.5; r<20.01; r+=0.5) {
    if (CKerr_getData_CircEq(M,astar,r,info))
    printf("%17.10lE %17.10lE %17.10lE %17.10lE %17.10lE %17.10lE %17.10lE %17.10lE\n", r, info[0], info[1], info[2], info[3], info[4], info[5], info[6]);
  }
#endif
}

void params_test(void) {
  double EQL[3], J[3], ancillary[3];
  double i,rp,ra;
  double astar = 0.9;

  scanf("%lg %lg %lg", &i, &rp, &ra);

  printf("[%d] ", CKerr_FindEQL_IRR(i,rp,ra,EQL,1,astar));
  printf(" EQL = %14.11lf %14.11lf %14.11lf\n", EQL[0], EQL[1], EQL[2]);
  CKerr_EQL2J(EQL,J,1,astar,ancillary);
  printf("     J   = %14.11lf %14.11lf %14.11lf\n", J[0], J[1], J[2]);
  printf("    anc  = %14.11lf %14.11lf %14.11lf\n", ancillary[0], ancillary[1], ancillary[2]);
  printf("I_Carter = %14.11lf\n", atan2(sqrt(EQL[1]), EQL[2]));
}

void gwkick(void) {
  double J[3], Minv[9], xuorig[6];
  double Ic; int fr;
  double M = 1, astar = 0.998;
  long lmax = 20;
  long mmax = 20;
  long kmax = 10;
  long Ntheta = 64;
  long qt,mm,nl,il,itheta,nb;
  long Delk, Delm;
  double *psiem, *psiptr, *psiptr2, *Coef, omega, *b, *S;
  double E[100], Om[3], QI, QH, Jj[3], ANCu[3], ANCd[3];
  double sintheta;
  double EQLu[3], EQLd[3], Ju[3], Jd[3], Omu[3], Omd[3], delta=1e-4;
  double r, rdot, alphaddot, jtdot=0;
  double Etot=0, Ltot=0, E_inf=0, E_em=0, px_em=0, py_em=0;

  scanf("%lg %d", &Ic, &fr);

  printf("[%d]\n", CKerr_FindResCirc(fr,Ic,J,M,astar));
  printf("     J   = %14.11lf %14.11lf %14.11lf\n", J[0], J[1], J[2]);

  CKerr_Minverse(J,Minv,M,astar);
  CKerr_TorusOrigin(J,xuorig,M,astar);
  CKerr_Minv2Omega(Minv,Om);
  r = xuorig[0];

  /* Emitted waveform is psi[Delm*(m+mmax) + Delk*(k+kmax) + 2*itheta] (re), +1 (im) */
  Delk = 2*Ntheta;
  Delm = Delk * (2*kmax+1);
  psiem = (double*)malloc((size_t)(2*(2*mmax+1)*(2*kmax+1)*Ntheta*sizeof(double)));
  Coef = (double*)malloc((size_t)(4*(2*lmax+1)*sizeof(double)));
  S = (double*)malloc((size_t)(2*(lmax+1)*sizeof(double)));

  /* Obtain emitted waveform at each angle: theta = (itheta+0.5)/Ntheta*M_PI */
  for(mm=-mmax; mm<=mmax; mm++) for(qt=-kmax;qt<=kmax;qt++) if (mm*fr+qt!=0) {

    nl=lmax-abs(mm)+1; if (abs(mm)<=2) nl=lmax-1;
    CKerr_RadialFunc(Minv, xuorig, M, astar, 0, qt, mm, 2, 3*kmax, nl, Coef, &omega, E);

    nb = CKerr_SbCoefN(-2,mm,nl,E,&b,M*astar*mm,1e-9);
    CKerr_QuantaEmitted(M,astar,0,qt,mm,nl,E,Coef,omega,&QI,&QH);
    Ltot += (QI+QH)*mm;
    Etot += (QI+QH)*omega;
    E_inf += QI*omega;
    jtdot += (QI+QH)*qt;

printf("# %3ld %3ld (%3ld %3ld)\n", mm, qt, nl, nb);

    for(itheta=0; itheta<Ntheta; itheta++) {
      CKerr_SpheroidalYSLM((itheta+0.5)/Ntheta*M_PI, -2, mm, nl, nb, b, S, S+nl);
      psiptr = psiem + Delm*(mm+mmax) + Delk*(qt+kmax) + 2*itheta;
      psiptr[0] = psiptr[1] = 0.;
      for(il=0;il<nl;il++) {
        psiptr[0] += S[il] * Coef[4*il+2];
        psiptr[1] += S[il] * Coef[4*il+3];
      }
    }

    free((char*)b);
  }
  else {
    /* Fill in zeroes for mm=qt=0 mode */
    for(itheta=0; itheta<Ntheta; itheta++) {
      psiptr = psiem + Delm*(mm+mmax) + Delk*(qt+kmax) + 2*itheta;
      psiptr[0] = psiptr[1] = 0.;
    }
  }

  printf("Mode sums: Ltot = %15.8le Etot = %15.8le Einf = %15.8le\n", Ltot, Etot, E_inf);

  /* Add up parameters for emitted wave forms */
  for(mm=-mmax; mm<=mmax; mm++) for(qt=-kmax;qt<=kmax;qt++) if (mm*fr+qt!=0) for(itheta=0; itheta<Ntheta;itheta++) {
    omega = mm*Om[2]+qt*Om[1];
    psiptr = psiem + Delm*(mm+mmax) + Delk*(qt+kmax) + 2*itheta;
    sintheta = sin((itheta+0.5)*M_PI/Ntheta);
    E_em += (psiptr[0]*psiptr[0]+psiptr[1]*psiptr[1])/(2.*omega*omega) * M_PI/(double)Ntheta * sintheta;

    if (mm>-mmax && qt<=kmax-fr) {
      psiptr2 = psiem + Delm*(mm-1+mmax) + Delk*(qt+fr+kmax) + 2*itheta;
      px_em += (psiptr[0]*psiptr2[0]+psiptr[1]*psiptr2[1])/(2.*omega*omega) * M_PI/(double)Ntheta/2. * sintheta * sintheta;
      py_em += (psiptr[0]*psiptr2[1]-psiptr[1]*psiptr2[0])/(2.*omega*omega) * M_PI/(double)Ntheta/2. * sintheta * sintheta;
    }
#if 0
    printf("%3ld %3ld %15.11lf %10.7lf %15.8le %15.8le\n", mm, qt, mm*Om[2]+qt*Om[1], (itheta+0.5)/Ntheta*M_PI, psiptr[0], psiptr[1]);
#endif
  }

  /* Report emitted power & force */
  printf("E_em = %15.8le px_em = %15.8le py_em = %15.8le\n", E_em, px_em, py_em);

  /* Report sink rate */
  CKerr_FindEQL_IRCirc(Ic,r+1e-4,EQLu,M,astar);
  CKerr_FindEQL_IRCirc(Ic,r-1e-4,EQLd,M,astar);
  rdot = -Ltot / (EQLu[2]-EQLd[2]) * 2e-4;
  printf("r = %15.8le rdot = %15.8le\n", r, rdot);

  /* Angular acceleration */
  CKerr_EQL2J(EQLu,J,M,astar,NULL);
  CKerr_Minverse(J,Minv,M,astar);
  CKerr_Minv2Omega(Minv,Omu);
  CKerr_EQL2J(EQLd,J,M,astar,NULL);
  CKerr_Minverse(J,Minv,M,astar);
  CKerr_Minv2Omega(Minv,Omd);
  alphaddot = ( (Omu[2]-Omd[2]) - 2*(Omu[1]-Omd[1]) ) / 2e-4 * rdot;
  printf("alphaddot = %15.8le\n", alphaddot);

  /* Net kick; should be multiplied by mu^1.5 */
  printf("pkick = %15.8le\n", sqrt(px_em*px_em+py_em*py_em) * sqrt(2.*M_PI/fabs(alphaddot)) );

  /* Adiabatic invariant method for angular acceleration */
  Ju[0] = J[0];
  Ju[1] = J[1] + delta*jtdot;
  Ju[2] = J[2] + delta*Ltot;
  Jd[0] = J[0];
  Jd[1] = J[1] - delta*jtdot;
  Jd[2] = J[2] - delta*Ltot;
  CKerr_Minverse(Ju,Minv,M,astar);
  CKerr_Minv2Omega(Minv,Omu);
  CKerr_J2EQL(Ju,EQLu,M,astar);
  CKerr_EQL2J(EQLu,Jj,M,astar,ANCu);
  CKerr_Minverse(Jd,Minv,M,astar);
  CKerr_Minv2Omega(Minv,Omd);
  CKerr_J2EQL(Jd,EQLd,M,astar);
  CKerr_EQL2J(EQLd,Jj,M,astar,ANCd);
  alphaddot = -( (Omu[2]-Omd[2]) - 2*(Omu[1]-Omd[1]) ) / (2.*delta);
  printf("rdot = %15.8le idot=%15.8le\n",
    -(ANCu[1]+ANCu[2]-ANCd[1]-ANCd[2]) /(4.*delta),
    -(atan2(sqrt(EQLu[1]),EQLu[2]) - atan2(sqrt(EQLd[1]),EQLd[2]))/(2.*delta)
  );
  printf("jdot = %15.8le,%15.8le,%15.8le alphaddot = %15.8le\n", 0., -jtdot, -Ltot, alphaddot);

  /* Net kick; should be multiplied by mu^1.5 */
  printf("pkick = %15.8le\n", sqrt(px_em*px_em+py_em*py_em) * sqrt(2.*M_PI/fabs(alphaddot)) );

  free((char*)psiem);
  free((char*)Coef);
  free((char*)S);
}

int main(void) {

#if 0
  double Ic, J[3];
  for(Ic=25.0/180.*M_PI; Ic<6./18.*M_PI+1e-5; Ic+=M_PI/180.*0.5) {
    CKerr_FindResCirc(2,Ic,J,1,.98);
  }
  return(0);
  for(Ic=0; Ic<M_PI*1.00001; Ic+= M_PI/180.) printf("%10.8lf %10.8lf\n", Ic, CKerr_FindISCO(Ic, 1, 0.98));
  return(0);
#endif

#if 0
  double R1[4],R3[4];
  CKerr_IngoingR1(500,1,0,0.0002,2,6,R1);
  CKerr_OutgoingR3(500,1,0,0.0002,2,6,R3);
  printf("R1(500) = %12.5le,%12.5le,%12.5le,%12.5le R3(500) = %12.5le,%12.5le,%12.5le,%12.5le\n", R1[0], R1[1], R1[2], R1[3], R3[0], R3[1], R3[2], R3[3]);
  return(0);
#endif
//  gwkick();
//return(0);
//  params_test();
  resonance_test();
//  gwscat_test();
//  gwem_test();
//  coord_test();
//  mode_test();
//  action_test();
  return(0);
}
