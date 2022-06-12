from libc.stdlib cimport malloc, free
cimport gr_wrapper
import numpy as np



# The documentation here covers primarily things concerning data management.
# If you were looking for math documentation, see "CKerr.h"

# Depends on libc.stdlib for everything and numpy for just CKerr_SbCoefN




# All the function wrappers from kerrtraj.c



#Inputs: an array of 3 doubles (arrJ), a double (M), a double from 0 to 1 (astar)
#Outputs: an array of 3 doubles (EQL), checker int (result)
#int CKerr_J2EQL(double *J, double *EQL, double M, double astar);

def ckerr_j2eql(arrJ, M, astar):
    assert len(arrJ) == 3

    cdef double *arrJ_cpy = <double *> malloc(len(arrJ)*sizeof(double))
    for i in range(len(arrJ)):
        arrJ_cpy[i] = <double> arrJ[i]

    cdef int result = -1 
    cdef double *EQL = <double *> malloc (3 * sizeof(double))
    try:
        result = gr_wrapper.CKerr_J2EQL(arrJ_cpy, EQL, <double> M, <double> astar)
    finally:
        output = [EQL[i] for i in range(3)]
        free(EQL)
        free(arrJ_cpy)

    return result, output





#Inputs: an array of 3 doubles (arrEQL), a double (M), a double from 0 to 1 (astar)
#Outputs: an array of 3 doubles (J), an array of 3 doubles (ancilliary), checker int (result)
#int CKerr_EQL2J(double *EQL, double *J, double M, double astar, double *ancillary);
def ckerr_eql2j(arrEQL, M, astar):
    assert len(arrEQL) == 3

    cdef double *arrEQL_cpy = <double *> malloc(len(arrEQL)*sizeof(double))
    for i in range(len(arrEQL)):
        arrEQL_cpy[i] = <double> arrEQL[i]

    cdef int result = -1 
    cdef double *J = <double *> malloc (3 * sizeof(double))
    cdef double *ancilliary = <double *> malloc (3 * sizeof(double))
    try:
        result = gr_wrapper.CKerr_EQL2J(arrEQL_cpy, J, <double> M, <double> astar, ancilliary)
    finally:
        output1 = [J[i] for i in range(3)]
        output2 = [ancilliary[i] for i in range(3)]
        free(J)
        free(ancilliary)
        free(arrEQL_cpy)

    return result, output1, output2






#Inputs: an array of 3 doubles (arrJ), a double (M), a double from 0 to 1 (astar)
#Outputs: an array of 9 doubles (Minv), checker int (result)
#int CKerr_Minverse(double *J, double *Minv, double M, double astar);
def ckerr_minverse(arrJ, M, astar):
    assert len(arrJ) == 3

    cdef double *arrJ_cpy = <double *> malloc(len(arrJ)*sizeof(double))
    for i in range(len(arrJ)):
        arrJ_cpy[i] = <double> arrJ[i]

    cdef int result = -1 
    cdef double *Minv = <double *> malloc (9 * sizeof(double))
    try:
        result = gr_wrapper.CKerr_Minverse(arrJ_cpy, Minv, <double> M, <double> astar)
    finally:
        output = [Minv[i] for i in range(9)]
        free(Minv)
        free(arrJ_cpy)
    
    return result, output



#Inputs: an array of 9 doubles (arrMinv)
#Outputs: an array of 3 doubles (Omega)
# void CKerr_Minv2Omega(double *Minv, double *Omega);
def ckerr_minv2omega(arrMinv):
    assert len(arrMinv) == 9

    cdef double *arrMinv_cpy = <double *> malloc(len(arrMinv)*sizeof(double))
    for i in range(len(arrMinv)):
        arrMinv_cpy[i] = <double> arrMinv[i]

    cdef double *Omega = <double *> malloc (3 * sizeof(double))

    try:
        gr_wrapper.CKerr_Minv2Omega(arrMinv_cpy, Omega)
    finally:
        output = [Omega[i] for i in range(3)]
        free(Omega)
        free(arrMinv_cpy)

    return output






#Inputs: a double (Q), a double (L), a double (M), a double from 0 to 1 (astar)
#Outputs: an array of 9 doubles (Minv), a double (emin)
#double CKerr_Emin(double Q, double L, double M, double astar)
def ckerr_emin(Q, L, M, astar):

    cdef double emax = -1
    emin = gr_wrapper.CKerr_Emin(<double> Q, <double> L, <double> M, <double> astar)
    return emin





#Inputs: a double (Q), a double (L), a double (M), a double from 0 to 1 (astar), an array of 1 double (arrJrmax)
#Outputs: a double (emax)
#double CKerr_Emax(double Q, double L, double M, double astar, double *Jrmax)
def ckerr_emax(Q, L, M, astar, arrJrmax):
    assert len(arrJrmax) == 1

    cdef double *arrJrmax_cpy = <double *> malloc(len(arrJrmax)*sizeof(double))
    for i in range(len(arrJrmax)):
        arrJrmax_cpy[i] = <double> arrJrmax[i]

    cdef double emax = -1
    try:
        emax = gr_wrapper.CKerr_Emax(<double> Q, <double> L, <double> M, <double> astar, arrJrmax_cpy)
    finally:
        free(arrJrmax_cpy)

    return emax




#Inputs: a double (Q), a double (L), a double (Jr), a double (M), a double from 0 to 1 (astar)
#Outputs: an array of 1 double (energy_value) (note that it returns Jtheta if unsuccesful)
#double CKerr_QLJr2E(double Q, double L, double Jr, double M, double astar, double *Jtheta)
def ckerr_qljr2e(Q, L, Jr, M, astar):

    cdef double *Jtheta = <double *> malloc (3 * sizeof(double))
    cdef double energy_value= -1
    try:
        energy_value = gr_wrapper.CKerr_QLJr2E(<double> Q, <double> L, <double> Jr, <double> M, <double> astar, Jtheta)
    finally:
        output = [Jtheta[i] for i in range(1)]
        free(Jtheta)

    return energy_value, output




#Inputs: an array of 3 doubles (arrJ), an array of 3 double (arrEQL), a double (M), a double from 0 to 1 (astar)
#Outputs: an array of 3 doubles (dEdJ), an array of 3 doubles (dQdJ), checker int (result)
#int CKerr_dEQdJ(double *J, double *EQL, double *dEdJ, double *dQdJ, double M, double astar)
def ckerr_deqdj(arrJ, arrEQL, M, astar):
    assert len(arrJ), len(arrEQL) == 3

    cdef double *arrJ_cpy = <double *> malloc(len(arrJ)*sizeof(double))
    for i in range(len(arrJ)):
        arrJ_cpy[i] = <double> arrJ[i]

    cdef double *arrEQL_cpy = <double *> malloc(len(arrEQL)*sizeof(double))
    for i in range(len(arrEQL)):
        arrEQL_cpy[i] = <double> arrEQL[i]

    cdef double *dEdJ = <double *> malloc (3 * sizeof(double))
    cdef double *dQdJ = <double *> malloc (3 * sizeof(double))
    try:
        result = gr_wrapper.CKerr_dEQdJ(arrJ_cpy, arrEQL_cpy, dEdJ, dQdJ, <double> M, <double> astar)
    finally:
        output1 = [dEdJ[i] for i in range(3)]
        output2 = [dQdJ[i] for i in range(3)]
        free(arrJ_cpy)
        free(arrEQL_cpy)
        free(dEdJ)
        free(dQdJ)

    return result, output1, output2
    



#Inputs: an array of 3 doubles (arrJ), a double (M), a double from 0 to 1 (astar)
#Outputs: an array of 9 doubles (Minv), checker int (result)
#int CKerr_MinverseND(double *J, double *Minv, double M, double astar)
def ckerr_minversend(arrJ, M, astar):
    assert len(arrJ) == 3

    cdef double *arrJ_cpy = <double *> malloc(len(arrJ)*sizeof(double))
    for i in range(len(arrJ)):
        arrJ_cpy[i] = <double> arrJ[i]

    cdef double *Minv = <double *> malloc (9 * sizeof(double))
    cdef int result = -1
    try:
        result = gr_wrapper.CKerr_MinverseND(arrJ_cpy, Minv, <double> M, <double> astar)
    finally:
        output = [Minv[i] for i in range(9)]
        free(Minv)
        free(arrJ_cpy)

    return result, output




#Inputs: an array of 3 doubles (arrJ), a double (M), a double from 0 to 1 (astar)
#Outputs: an array of 6 doubles (xu), checker int (result)
#int CKerr_TorusOrigin(double *J, double *xu, double M, double astar)
def ckerr_torusorigin(arrJ, M, astar):

    cdef double *arrJ_cpy = <double *> malloc(len(arrJ)*sizeof(double))
    for i in range(len(arrJ)):
        arrJ_cpy[i] = <double> arrJ[i]

    cdef double *xu = <double *> malloc (6 * sizeof(double))
    cdef int result = -1
    try:
        result = gr_wrapper.CKerr_TorusOrigin(arrJ_cpy, xu, <double> M, <double> astar)
    finally:
        output = [xu[i] for i in range(6)]
        free(xu)
        free(arrJ_cpy)

    return result, output





#Inputs: an array of 6 doubles (arrxu), a double (M), a double from 0 to 1 (astar)
#Outputs: an array of 6 doubles (xudot_H), an array of 5 doubles (xudot_Q)
#void CKerr_FlowVectors(double *xu, double *xudot_H, double *xudot_Q, double M, double astar)
def ckerr_flowvectors(arrxu, M, astar):
    assert len(arrxu) == 6

    cdef double *arrxu_cpy = <double *> malloc(len(arrxu)*sizeof(double))
    for i in range(len(arrxu)):
        arrxu_cpy[i] = <double> arrxu[i]

    cdef double *xudot_H = <double *> malloc (6 * sizeof(double))
    cdef double *xudot_Q = <double *> malloc (5 * sizeof(double))

    try:
        gr_wrapper.CKerr_FlowVectors(arrxu_cpy, xudot_H, xudot_Q, <double> M, <double> astar)
    finally:
        output1 = [xudot_H[i] for i in range(6)]
        output2 = [xudot_Q[i] for i in range(5)]
        free(xudot_H)
        free(xudot_Q)
        free(arrxu_cpy)

    return output1, output2



#Inputs: an array of 6 doubles (arrxu), an array of 9 double (arrMinv), a double (M), a double from 0 to 1 (astar)
#Outputs: an array of 6 doubles (xudot_r), an array of 6 doubles (xudot_theta)
#void CKerr_AngleFlowVectors(double *xu, double *xudot_r, double *xudot_theta, double *Minv, double M, double astar)
def ckerr_angleflowvectors(arrxu, arrMinv, M, astar):
    assert len(arrxu) == 6
    assert len(arrMinv) == 9

    cdef double *arrxu_cpy = <double *> malloc(len(arrxu)*sizeof(double))
    for i in range(len(arrxu)):
        arrxu_cpy[i] = <double> arrxu[i]

    cdef double *arrMinv_cpy = <double *> malloc(len(arrMinv)*sizeof(double))
    for i in range(len(arrMinv)):
        arrMinv_cpy[i] = <double> arrMinv[i]

    cdef double *xudot_r = <double *> malloc (6 * sizeof(double))
    cdef double *xudot_theta = <double *> malloc (6 * sizeof(double))

    try:
        gr_wrapper.CKerr_AngleFlowVectors(arrxu_cpy, xudot_r, xudot_theta, arrMinv_cpy, <double> M, <double> astar)
    finally:
        output1 = [xudot_r[i] for i in range(6)]
        output2 = [xudot_theta[i] for i in range(6)]
        free(xudot_r)
        free(xudot_theta)
        free(arrxu_cpy)

    return output1, output2







#Inputs: a double (M), a double from 0 to 1 (astar), a double (r)
#Outputs: an array of 7 doubles (info), checker int (result)
#int CKerr_getData_CircEq(double M, double astar, double r, double *info)
def ckerr_getdata_circeq(M, astar, r):

    cdef double *info = <double *> malloc (7 * sizeof(double))
    cdef int result = -1
    try:
        result = gr_wrapper.CKerr_getData_CircEq(<double> M, <double> astar, <double> r, info)
    finally:
        output = [info[i] for i in range(7)]
        free(info)
        
    return result, output



#Inputs: a double (M), a double from 0 to 1 (astar), a double (r0), a double (mm)
#Outputs: a double (lindbladresonance)
#double CKerr_FindLindbladResonance(double M, double astar, double r0, double mm)
def ckerr_findlindbladresonance(M, astar, r0, mm):

    cdef double lindbladresonance = -1
    lindbladresonance = gr_wrapper.CKerr_FindLindbladResonance(<double> M, <double> astar, <double> r0, <double> mm)
    return lindbladresonance




#Inputs: a double (Ic), a double (rp), a double (ra), a double (M), a double from 0 to 1 (astar)
#Outputs: an array of 3 doubles (EQL), checker int (result)
#int CKerr_FindEQL_IRR(double Ic, double rp, double ra, double *EQL, double M, double astar)
def ckerr_findeql_irr(Ic, rp, ra, M, astar):

    cdef double *EQL = <double *> malloc (3 * sizeof(double))
    cdef int result = -1
    try:
        result = gr_wrapper.CKerr_FindEQL_IRR(<double> Ic, <double> rp, <double> ra, EQL, <double> M, <double> astar)
    finally:
        output = [EQL[i] for i in range(3)]
        free(EQL)
        
    return result, output



#Inputs: a double (Ic), a double (r), a double (M), a double from 0 to 1 (astar)
#Outputs: an array of 3 doubles (EQL), checker int (result)
#int CKerr_FindEQL_IRCirc(double Ic, double r, double *EQL, double M, double astar)
def ckerr_findeql_ircirc(Ic, r, M, astar):

    cdef double *EQL = <double *> malloc (3 * sizeof(double))
    cdef int result = -1
    try:
        result = gr_wrapper.CKerr_FindEQL_IRCirc(<double> Ic, <double> r, EQL, <double> M, <double> astar)
    finally:
        output = [EQL[i] for i in range(3)]
        free(EQL)
        
    return result, output



#Inputs: a double >= 2 (fratio), a double (Ic), a double (M), a double from 0 to 1 (astar)
#Outputs: an array of 3 doubles (J), checker int (result)
# int CKerr_FindResCirc(int fratio, double Ic, double *J, double M, double astar)
def ckerr_findrescirc(fratio, Ic, M, astar):

    cdef double *J = <double *> malloc (3 * sizeof(double))
    cdef int result = -1
    try:
        result = gr_wrapper.CKerr_FindResCirc(<int> fratio, <double> Ic, J, <double> M, <double> astar)
    finally:
        output = [J[i] for i in range(3)]
        free(J)
        
    return result, output



#Inputs: a double (Ic), a double (M), a double from 0 to 1 (astar)
#Outputs: a double (isco_radius)
#double CKerr_FindISCO(double Ic, double M, double astar)
def ckerr_findisco(Ic, M, astar):

    cdef double isco_radius = -1
    isco_radius = gr_wrapper.CKerr_FindISCO(<double> Ic, <double> M, <double> astar)
    return isco_radius
















# All the function wrappers from kerrmode.c




#Inputs: a double (theta), a long (l), a long (s)
#Outputs: an array of 2*l+1 doubles (Y), an array of 2*l+1 doubles (dYdtheta) 
#void CKerr_Yslm(double theta, long l, long s, double *Y, double *dYdtheta);
def ckerr_yslm(theta, l, s):

    cdef int n = <int> l
    cdef double *Y = <double *> malloc ((2 * n + 1) * sizeof(double))
    cdef double *dYdtheta= <double *> malloc ((2 * n + 1) * sizeof(double))
    
    try:
        gr_wrapper.CKerr_Yslm(<double> theta, <long> l, <long> s, Y, dYdtheta)
    finally:
        output1 = [Y[i] for i in range(2 * l + 1)]
        output2 = [dYdtheta[i] for i in range(2 * l + 1)]
        free(Y)
        free(dYdtheta)

    return output1, output2







#Inputs: a long (s), a long (m), a long (nb), a double (chi), a double (errtol) 
#Outputs: an array of nb doubles (E), an array of nb doubles (b), a long (eigenvec_num)
#long CKerr_SbCoef(long s, long m, long nb, double *E, double *b, double chi, double errtol)
def ckerr_sbcoef(s, m, nb, chi, errtol):

    cdef double *E = <double *> malloc (nb * sizeof(double))
    cdef double *b= <double *> malloc (nb * sizeof(double))

    cdef long eigenvec_num = -1
    try:
        eigenvec_num = gr_wrapper.CKerr_SbCoef(<long> s, <long> m, <long> nb, E, b, <double> chi, <double> errtol)
    finally:
        output1 = [E[i] for i in range(nb)]
        output2 = [b[i] for i in range(nb)]
        free(E)
        free(b)

    return eigenvec_num, output1, output2




   



#Inputs: a long (s), a long (m), a long (nu), a double (chi), a double (errtol) 
#Outputs: an array of nu doubles (E), a long (nb), a double an numpy array of nb*nb doubles (output_array)
#long CKerr_SbCoefN(long s, long m, long nu, double *E, double **b, double chi, double errtol)
def ckerr_sbcoefn(s, m, nu, arrE, chi, errtol):

    cdef double *arrE_cpy = <double *> malloc(len(arrE)*sizeof(double))
    for i in range(len(arrE)):
        arrE_cpy[i] = <double> arrE[i]

    cdef double **my_b = <double **> malloc(sizeof(double *))
    try:
        nb = gr_wrapper.CKerr_SbCoefN(<long> s, <long> m, <long> nu, arrE_cpy, my_b, <double> chi, <double> errtol)
    finally:
        output_array = np.zeros((nb*nb))
        for i in range(nb*nb):
            output_array[i] = my_b[0][i]
        output = [arrE_cpy[i] for i in range(nu)]
        free(my_b)
        free(my_b[0])
        free(arrE_cpy)

    return output_array.reshape((nb,nb)), output, nb
  











#Inputs: a double (theta), a long (s), a long (m), a long (nu), a double (chi), a double (errtol), an array of nu*nb doubles (arrb)
#Outputs: an array of nu doubles (S), an array of nu doubles (dSdtheta)
#void CKerr_SpheroidalYSLM(double theta, long s, long m, long nu, long nb, double *b, double *S, double *dSdtheta)
def ckerr_spheroidalyslm(theta, s, m, nu, nb, arrb):

    cdef double *arrb_cpy = <double *> malloc(len(arrb)*sizeof(double))
    for i in range(len(arrb)):
        arrb_cpy[i] = <double> arrb[i]

    cdef double *S = <double *> malloc (nu * sizeof(double))
    cdef double *dSdtheta = <double *> malloc (nu * sizeof(double))
    try:
        gr_wrapper.CKerr_SpheroidalYSLM(<double> theta, <long> s, <long> m, <long> nu, <long> nb, arrb_cpy, S, dSdtheta)
    finally:
        output1 = [S[i] for i in range(nu)]
        output2 = [dSdtheta[i] for i in range(nu)]
        free(arrb_cpy)
        free(S)
        free(dSdtheta)

    return output1, output2





#Inputs: a double (r), a double (M), a double from 0 to 1 (astar)
#Outputs: a double (rstar)
#double CKerr_r2rstar(double r, double M, double astar) 
def ckerr_r2rstar(r, M, astar):

    cdef double rstar = -1
    rstar = gr_wrapper.CKerr_r2rstar(<double> r, <double> M, <double> astar)
    return rstar







#Inputs: a double (rstar), a double (M), a double from 0 to 1 (astar)
#Outputs: a double (r), an array of 1024 doubles (Delta)
#Why the value of 1024? It's the "CKERR_NSTEP_TORTOISE" defined in kerrmode.c
#double CKerr_rstar2r(double rstar, double M, double astar, double *Delta)
def ckerr_rstar2r(rstar, M, astar):

    cdef double *Delta = <double *> malloc (1024 * sizeof(double))
    cdef result = -1
    try:
        result = gr_wrapper.CKerr_rstar2r(<double> rstar, <double> M, <double> astar, Delta)
    finally:
        output = [Delta[i] for i in range(1024)]
        free(Delta)

    return result, output










#Inputs: a double (rstar), a double (M), a double from 0 to 1 (astar), a double (omega), a long (mm), a double (eigenlm)
#Outputs: an array of 1 double (ReV), an array of 1 double (ImV)
#void CKerr_Vr(double rstar, double M, double astar, double omega, long mm, double eigenlm, double *ReV, double *ImV)
def ckerr_vr(rstar, M, astar, omega, mm, eigenlm):

    cdef double *ReV = <double *> malloc (1 * sizeof(double))
    cdef double *ImV = <double *> malloc (1 * sizeof(double))
    try:
        gr_wrapper.CKerr_Vr(<double> rstar, <double> M, <double> astar, <double> omega, <long> mm, <double> eigenlm, ReV, ImV)
    finally:
        output1 = [ReV[i] for i in range(1)]
        output2 = [ImV[i] for i in range(1)]
        free(ReV)
        free(ImV)

    return output1, output2







#Inputs: a double (rf), a double (M), a double from 0 to 1 (astar), a double (omega), a long (mm), a double (eigenlm)
#Outputs: an array of 4 doubles (R1)
#void CKerr_IngoingR1(double rf, double M, double astar, double omega, long mm, double eigenlm, double *R1)
def ckerr_ingoingr1(rf, M, astar, omega, mm, eigenlm):

    cdef double *R1 = <double *> malloc (4 * sizeof(double))
    try:
        gr_wrapper.CKerr_IngoingR1(<double> rf, <double> M, <double> astar, <double> omega, <long> mm, <double> eigenlm, R1)
    finally:
        output = [R1[i] for i in range(1)]
        free(R1)

    return output







#Inputs: an array of 2 doubles (arr_r), an array of 4 doubles (arr_RR), a double (M), a double from 0 to 1 (astar), a double (omega), a long (mm), a double (eigenlm)
#Outputs: an array of 2 doubls (d2Rdr2)
#void CKerr_Radial2nd(double *r, double M, double astar, double omega, long mm, double eigenlm, double *RR, double *d2Rdr2)
def ckerr_radial2nd(arr_r, arr_RR, M, astar, omega, mm, eigenlm):
    assert len(arr_r) == 2
    assert len(arr_RR) == 4

    cdef double *arr_r_cpy = <double *> malloc(len(arr_r)*sizeof(double))
    for i in range(len(arr_r)):
        arr_r_cpy[i] = <double> arr_r[i]

    cdef double *arr_RR_cpy = <double *> malloc(len(arr_RR)*sizeof(double))
    for i in range(len(arr_RR)):
        arr_RR_cpy[i] = <double> arr_RR[i]

    cdef double *d2Rdr2 = <double *> malloc (2 * sizeof(double))
    try:
        gr_wrapper.CKerr_Radial2nd(arr_r_cpy, <double> M, <double> astar, <double> omega, <long> mm, <double> eigenlm, arr_RR_cpy, d2Rdr2)
    finally:
        output = [d2Rdr2[i] for i in range(2)]
        free(arr_r_cpy)
        free(arr_RR_cpy)
        free(d2Rdr2)

    return output






#Inputs: an array of 2 doubles (arr_r), an array of 2 doubles (arr_h), a double (M), a double from 0 to 1 (astar), a double (omega), a long (mm), a double (eigenlm), an array of 4 doubles (RR_old)
#Outputs: an array of 4 doubles (RR_new)
#void CKerr_RadialRKStep(double *r, double *h, double M, double astar, double omega, long mm, double eigenlm, double *RR_old, double *RR_new)
def ckerr_radialrkstep(arr_r, arr_h, M, astar, omega, mm, eigenlm, arr_RR_old):

    cdef double *arr_r_cpy = <double *> malloc(len(arr_r)*sizeof(double))
    for i in range(len(arr_r)):
        arr_r_cpy[i] = <double> arr_r[i]

    cdef double *arr_h_cpy = <double *> malloc(len(arr_h)*sizeof(double))
    for i in range(len(arr_h)):
        arr_h_cpy[i] = <double> arr_h[i]

    cdef double *arr_RR_old_cpy = <double *> malloc(len(arr_RR_old)*sizeof(double))
    for i in range(len(arr_RR_old)):
        arr_RR_old_cpy[i] = <double> arr_RR_old[i]

    cdef double *RR_new = <double *> malloc (4 * sizeof(double))
    try:
        gr_wrapper.CKerr_RadialRKStep(arr_r_cpy, arr_h_cpy, <double> M, <double> astar, <double> omega, <long> mm, <double> eigenlm, arr_RR_old_cpy, RR_new)
    finally:
        output = [RR_new[i] for i in range(4)]
        free(arr_r_cpy)
        free(arr_h_cpy)
        free(arr_RR_old_cpy)
        free(RR_new)

    return output









#Inputs: an array of 2 doubles (arr_r), an array of 2 doubles (arr_h), a double (M), a double from 0 to 1 (astar), a double (omega), a long (mm), a double (eigenlm), a long w/ >= 2 (N)
#Outputs: an array of 4 doubles (RR_new)
#void CKerr_RadialMidptStep(double *r, double *h, double M, double astar, double omega, long mm, double eigenlm, double *RR_old, double *RR_new, long N) {
def ckerr_radialmidptstep(arr_r, arr_h, M, astar, omega, mm, eigenlm, arr_RR_old, N):

    cdef double *arr_r_cpy = <double *> malloc(len(arr_r)*sizeof(double))
    for i in range(len(arr_r)):
        arr_r_cpy[i] = <double> arr_r[i]

    cdef double *arr_h_cpy = <double *> malloc(len(arr_h)*sizeof(double))
    for i in range(len(arr_h)):
        arr_h_cpy[i] = <double> arr_h[i]

    cdef double *arr_RR_old_cpy = <double *> malloc(len(arr_RR_old)*sizeof(double))
    for i in range(len(arr_RR_old)):
        arr_RR_old_cpy[i] = <double> arr_RR_old[i]

    cdef double *RR_new = <double *> malloc (4 * sizeof(double))
    try:
        gr_wrapper.CKerr_RadialMidptStep(arr_r_cpy, arr_h_cpy, <double> M, <double> astar, <double> omega, <long> mm, <double> eigenlm, arr_RR_old_cpy, RR_new, <long> N)
    finally:
        output = [RR_new[i] for i in range(4)]
        free(arr_r_cpy)
        free(arr_h_cpy)
        free(arr_RR_old_cpy)
        free(RR_new)

    return output








#Inputs: an array of 2 doubles (arr_r), an array of 4 doubles (arr_RR), a double (M), a double from 0 to 1 (astar), a double (omega), a long (mm), a double (eigenlm), an array of 4 doubles (arr_RR_old)
#Outputs: an array of 4 doubles (RR_new)
#void CKerr_RadialMidptExtrStep(double *r, double *h, double M, double astar, double omega, long mm, double eigenlm, double *RR_old, double *RR_new)
def ckerr_radialmidptextrstep(arr_r, arr_h, M, astar, omega, mm, eigenlm, arr_RR_old):

    cdef double *arr_r_cpy = <double *> malloc(len(arr_r)*sizeof(double))
    for i in range(len(arr_r)):
        arr_r_cpy[i] = <double> arr_r[i]

    cdef double *arr_h_cpy = <double *> malloc(len(arr_h)*sizeof(double))
    for i in range(len(arr_h)):
        arr_h_cpy[i] = <double> arr_h[i]

    cdef double *arr_RR_old_cpy = <double *> malloc(len(arr_RR_old)*sizeof(double))
    for i in range(len(arr_RR_old)):
        arr_RR_old_cpy[i] = <double> arr_RR_old[i]

    cdef double *RR_new = <double *> malloc (4 * sizeof(double))
    try:
        gr_wrapper.CKerr_RadialMidptExtrStep(arr_r_cpy, arr_h_cpy, <double> M, <double> astar, <double> omega, <long> mm, <double> eigenlm, arr_RR_old_cpy, RR_new)
    finally:
        output = [RR_new[i] for i in range(4)]
        free(arr_r_cpy)
        free(arr_h_cpy)
        free(arr_RR_old_cpy)
        free(RR_new)

    return output









#Inputs: a double (rf), a double (M), a double from 0 to 1 (astar), a double (omega), a long (mm), a double (eigenlm)
#Outputs: an array of 4 doubles (R3)
#void CKerr_OutgoingR3(double rf, double M, double astar, double omega, long mm, double eigenlm, double *R3)
def ckerr_outgoingr3(rf, M, astar, omega, mm, eigenlm):

    cdef double *R3 = <double *> malloc (4 * sizeof(double))
    try:
        gr_wrapper.CKerr_OutgoingR3(<double> rf, <double> M, <double> astar, <double> omega, <long> mm, <double> eigenlm, R3)
    finally:
        output = [R3[i] for i in range(4)]
        free(R3)

    return output










#Inputs: an array of 2 doubles (arr_r), an array of 2 doubles (arr_RR), a double (M), a double from 0 to 1 (astar), a double (omega), a long (mm), a double (eigenlm), an array of 4 doubles (arr_RR_old), a long (N)
#Outputs: an array of 4 doubles (RR_new)
#void CKerr_RadialStepMany(double *r_old, double *r_new, double M, double astar, double omega, long mm, double eigenlm, double *RR_old, double *RR_new, long Nstep)
def ckerr_radialstepmany(arr_r_old, arr_r_new, M, astar, omega, mm, eigenlm, arr_RR_old, Nstep):

    cdef double *arr_r_old_cpy = <double *> malloc(len(arr_r_old)*sizeof(double))
    for i in range(len(arr_r_old)):
        arr_r_old_cpy[i] = <double> arr_r_old[i]

    cdef double *arr_r_new_cpy = <double *> malloc(len(arr_r_new)*sizeof(double))
    for i in range(len(arr_r_new)):
        arr_r_new_cpy[i] = <double> arr_r_new[i]

    cdef double *arr_RR_old_cpy = <double *> malloc(len(arr_RR_old)*sizeof(double))
    for i in range(len(arr_RR_old)):
        arr_RR_old_cpy[i] = <double> arr_RR_old[i]

    cdef double *RR_new = <double *> malloc (4 * sizeof(double))
    try:
        gr_wrapper.CKerr_RadialStepMany(arr_r_old_cpy, arr_r_new_cpy, <double> M, <double> astar, <double> omega, <long> mm, <double> eigenlm, arr_RR_old_cpy, RR_new, <long> Nstep)
    finally:
        output = [RR_new[i] for i in range(4)]
        free(arr_r_old_cpy)
        free(arr_r_new_cpy)
        free(arr_RR_old_cpy)
        free(RR_new)

    return output







#Inputs: a double (M), a double from 0 to 1 (astar), a double (omega), a long (mm), a double (eigenlm)
#Outputs: an array of 16 doubles (cscat), an array of 4 doubles (aux)
#void CKerr_GWScatMatrix(double M, double astar, double omega, long mm, double eigenlm, double *cscat, double *aux) 
def ckerr_gwscatmatrix(M, astar, omega, mm, eigenlm):

    cdef double *cscat = <double *> malloc (16 * sizeof(double))
    cdef double *aux = <double *> malloc (4 * sizeof(double))
    try:
        gr_wrapper.CKerr_GWScatMatrix(<double> M, <double> astar, <double> omega, <long> mm, <double> eigenlm, cscat, aux)
    finally:
        output1 = [cscat[i] for i in range(16)]
        output2 = [aux[i] for i in range(4)]
        free(cscat)
        free(aux)

    return output1, output2











































































#Inputs: an array of 6 doubles (arr_xu), a double (M), a double from 0 to 1 (astar)
#Outputs: an array of 6 doubles (T)
#void CKerr_xu2intT(double *xu, double *T, double M, double astar)
def ckerr_xu2intt(arrxu, M, astar):
    assert len(arrxu) == 6

    cdef double *arrxu_cpy = <double *> malloc(len(arrxu)*sizeof(double))
    for i in range(len(arrxu)):
        arrxu_cpy[i] = <double> arrxu[i]

    cdef double *T = <double *> malloc (6 * sizeof(double))

    try:
        gr_wrapper.CKerr_xu2intT(arrxu_cpy, T, <double> M, <double> astar)
    finally:
        output = [T[i] for i in range(6)]
        free(T)
        free(arrxu_cpy)

    return output




#Inputs: an array of 9 doubles (arrMinv), an array of 6 (arrxuorig), a double (M), a double from 0 to 1 (astar), a long (qr), a long (qt), a long (mm), a long (Nr), a long (Nt), a long (nl)
#Outputs: an array of 4*nl doubles (Coef), an array of 1 double (omega), an array of nl doubles (E), checker int (result)
#int CKerr_RadialFunc(double *Minv, double *xuorig, double M, double astar, long qr, long qt, long mm, long Nr, long Nt, long nl, double *Coef, double *omega, double *E)
def ckerr_radialfunc(arrMinv, arrxuorig, M, astar, qr, qt, mm, Nr, Nt, nl):
    assert len(arrMinv) == 9
    assert len(arrxuorig) == 6

    cdef double *arrMinv_cpy = <double *> malloc(len(arrMinv)*sizeof(double))
    for i in range(len(arrMinv)):
        arrMinv_cpy[i] = <double> arrMinv[i]

    cdef double *arrxuorig_cpy = <double *> malloc(len(arrxuorig)*sizeof(double))
    for i in range(len(arrxuorig)):
        arrxuorig_cpy[i] = <double> arrxuorig[i]

    cdef int result = -1
    cdef double *Coef = <double *> malloc ((4*nl) * sizeof(double))
    cdef double *omega = <double *> malloc (1 * sizeof(double))
    cdef double *E = <double *> malloc (nl * sizeof(double))
    try:
        result = gr_wrapper.CKerr_RadialFunc(arrMinv_cpy, arrxuorig_cpy, <double> M, <double> astar, <long> qr, <long> qt, <long> mm, <long> Nr, <long> Nt, <long> nl, Coef, omega, E)
    finally:
        output1 = [Coef[i] for i in range(4*nl)]
        output2 = [omega[i] for i in range(1)]
        output3 = [E[i] for i in range(nl)]
        free(Coef)
        free(omega)
        free(E)
        free(arrMinv_cpy)
        free(arrxuorig_cpy)

    return result, output1, output2, output3





# Check!
#Inputs: a double (M), a double from 0 to 1 (astar), a double (omega), a long (mm), a double (eigenlm)
#Outputs: a double (alpha_coeff)
#double CKerr_HorizonFluxNorm(double M, double astar, double omega, long mm, double eigenlm);
def ckerr_horizonfluxnorm(M, astar, omega, mm, eigenlm):

    cdef double alpha_coeff = -1
    alpha_coeff  = gr_wrapper.CKerr_HorizonFluxNorm(<double> M, <double> astar, <double> omega, <long> mm, <double> eigenlm)
    return alpha_coeff 




#Inputs: a double (M), a double from 0 to 1 (astar), a long (qr), a long (qt), a long (mm), a double (nl), an array of 3 doubles (arrE), an array of 4*nl doubles (arrcoef)
#Outputs: an array of nl doubles (QuantaHor), an array of nl doubles (QuantaInf)
#void CKerr_QuantaEmitted(double M, double astar, long qr, long qt, long mm, long nl, double *E, double *Coef, double omega, double *QuantaInf, double *QuantaHor);
def ckerr_quantaemitted(M, astar, qr, qt, mm, nl, arrE, arrcoef, omega):

    cdef double *arrE_cpy = <double *> malloc(len(arrE)*sizeof(double))
    for i in range(len(arrE)):
        arrE_cpy[i] = <double> arrE[i]

    cdef double *arrcoef_cpy = <double *> malloc(len(arrcoef)*sizeof(double))
    for i in range(len(arrcoef)):
        arrcoef_cpy[i] = <double> arrcoef[i]

    cdef double *QuantaInf = <double *> malloc (nl * sizeof(double))
    cdef double *QuantaHor = <double *> malloc (nl * sizeof(double))
    try:
        gr_wrapper.CKerr_QuantaEmitted(<double> M, <double> astar, <long> qr, <long> qt, <long> mm, <long> nl, arrE_cpy, arrcoef_cpy, <double> omega, QuantaInf, QuantaHor)
    finally:
        output1 = [QuantaHor[i] for i in range(nl)]
        output2 = [QuantaInf[i] for i in range(nl)]
        free(QuantaHor)
        free(QuantaInf)
        free(arrE_cpy)
        free(arrcoef_cpy)

    return output1, output2








#Inputs: a double (M), a double from 0 to 1 (astar), a double (r0), a long (mm), a long (nl), a double (jr1)
#Outputs: an array of 1 double (r1), checker double (result)
#double CKerr_LindbladResonanceStrength(double M, double astar, double r0, long mm, long nl, double jr1, double *r1);
def ckerr_lindbladresonancestrength(M, astar, r0, mm, nl, jr1):

    cdef double *r1 = <double *> malloc (1 * sizeof(double))
    cdef double result = -1
    try:
        result = gr_wrapper.CKerr_LindbladResonanceStrength(<double> M, <double> astar, <double> r0, <long> mm, <long> nl, <double> jr1, r1)
    finally:
        output = [r1[i] for i in range(1)]
        free(r1)

    return result, output
         
    









#All the function wrappers from resonance_find.c



#Inputs: a double (ra), a double (rp), a double (I), a double from 0 to 1 (astar), a double (M)
#Outputs: an array of 3 doubles (EQL), checker double (result)
#double ra_rp_I2EQL(double ra, double *EQL, double rp, double I, double astar, double M)
def ra_rp_i2eql(ra, rp, I, astar, M):

    cdef double *EQL = <double *> malloc (3 * sizeof(double))
    cdef double result = -1
    try:
        result = gr_wrapper.ra_rp_I2EQL(<double> ra, EQL, <double> rp, <double> I, <double> M, <double> astar)
    finally:
        output = [EQL[i] for i in range(3)]
        free(EQL)

    return result, output





#Inputs: a double (radius), a double (M), a double (spin)
#Outputs: a double (omega_outer)
#double Omega_outer_direct(double radius, double M, double spin)
def omega_outer_direct(radius, M, spin):

    cdef double omega_outer = -1
    omega_outer = gr_wrapper.Omega_outer_direct(<double> radius, <double> M, <double> spin)
    return omega_outer




#Inputs: an int (n), an int (k), an int (m), a double (radius), a double (ra), a double (rp), a double I, a double from 0 to 1 (astar), a double (M)
#Outputs: a double (om_outer_circ)
#double ra_rp_I2Omega_OuterCirc(int n, int k, int m, double radius, double ra, double rp, double I, double astar, double M);
def ra_rp_i2omega_outercirc(n, k, m, radius, ra, rp, I, astar, M):

    cdef double om_outer_circ = -1
    om_outer_circ = gr_wrapper.ra_rp_I2Omega_OuterCirc(<int> n, <int> k, <int> m, <double> radius, <double> ra, <double> rp, <double> I, <double> astar, <double> M)
    return om_outer_circ



#Inputs: an int (n_inner), an int (k_inner), an int (m_inner), an int (n_outer), an int (k_outer), a double (ra_inner), a double (rp_inner), a double (I_inner), a double (ra_outer), a double (rp_outer), a double (I_outer), a double from 0 to 1 (astar), a double (M)
#Outputs: a double (om_res_cond)
#double ra_rp_I2Omega_generic(int n_inner, int k_inner, int m_inner, int n_outer, int k_outer, double ra_inner, double rp_inner, double I_inner, double ra_outer, double rp_outer, double I_outer, double astar, double M);
def ra_rp_i2omega_generic(n_inner, k_inner, m_inner, n_outer, k_outer, ra_inner, rp_inner, I_inner, ra_outer, rp_outer, I_outer, astar, M):

    cdef double om_generic = -1
    om_generic = gr_wrapper.ra_rp_I2Omega_generic(<int> n_inner, <int> k_inner, <int> m_inner, <int> n_outer, <int> k_outer, <double> ra_inner, <double> rp_inner, <double> I_inner, <double> ra_outer, <double> rp_outer, <double> I_outer, <double> astar, <double> M)
    return om_generic



#Inputs: an int (n), an int (k), an int (m), a double (radius), a double (guess1), a double (guess2), a double (rp), a double I, a double from 0 to 1 (astar), a double (M)
#Outputs: a double (res_apo)
#double find_resonance_apo_OuterCirc(int n, int k, int m, double radius, double guess1, double guess2, double rp, double I, double astar, double M)
def find_resonance_apo_outercirc(n, k, m, radius, guess1, guess2, rp, I, astar, M):

    cdef double res_apo = -1
    res_apo = gr_wrapper.find_resonance_apo_OuterCirc(<int> n, <int> k, <int> m, <double> radius, <double> guess1, <double> guess2, <double> rp, <double> I, <double> astar, <double> M)
    return res_apo







#All the function wrappers from J_dot.c



#Inputs: an int (nl), an int (nmax), an int (kmax), an int (m), a double (apo), a double (rp), a double (radius_outer), a double (I), a double (M), a double from 0 to 1 (astar)
#Outputs: a double (self_force), an array of 3 doubles (J_dot_sf)
#int J_dot_selfforce(int nl, int nmax, int kmax, int mmax, double apo, double rp, double radius_outer, double I, double M, double astar, double *J_dot_sf)
def j_dot_selfforce(nl, nmax, kmax, mmax, apo, rp, radius_outer, I, M, astar):

    cdef double *J_dot_sf = <double *> malloc (3 * sizeof(double))
    cdef int self_force = -1
    try:
        self_force = gr_wrapper.J_dot_selfforce(<int> nl, <int> nmax, <int> kmax, <int> mmax, <double> apo, <double> rp, <double> radius_outer, <double> I, <double> M, <double> astar, J_dot_sf)
    finally:
        output = [J_dot_sf[i] for i in range(3)]
        free(J_dot_sf)

    return self_force, output




#Inputs: an int (nl), an int (N_res), an int (n_res_inner), an int (n_res_outer), an int (k_res_inner), an int (k_res_outer), an int (m_res_inner), an int (m_res_outer), a double (apo), a double (rp), a double (radius_outer), a double (I), a double (M), a double from 0 to 1 (astar), a double (theta_res_F), an array of 3 doubles (J_dot_td)
#Outputs: a double (j_dot_t)
#int J_dot_tidal(int nl, int N_res, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double apo, double rp, double radius_outer, double I, double M, double astar, double theta_res_F, double *J_dot_td)
def j_dot_tidal(nl, N_res, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, apo, rp, radius_outer, I, M, astar, theta_res_F):

    cdef double *J_dot_td = <double *> malloc (3 * sizeof(double))
    cdef int j_dot_t = -1
    try:
        j_dot_t = gr_wrapper.J_dot_tidal(<int> nl, <int> N_res, <int> n_res_inner, <int> n_res_outer, <int> k_res_inner, <int> k_res_outer, <int> m_res_inner, <int> m_res_outer, <double> apo, <double> rp, <double> radius_outer, <double> I, <double> M, <double> astar, <double> theta_res_F, J_dot_td)
    finally:
        output = [J_dot_td[i] for i in range(3)]
        free(J_dot_td)

    return j_dot_t, output





#Inputs: a double (mu_outer), a double (r_outer), a double (apo), a double (pari), a double (incline), a double (Theta_res)
#Outputs: a double (torque_res)
#double J_dot_phi_Kepler(double mu_outer, double r_outer, double apo, double peri, double incline, double Theta_res)
def j_dot_phi_kepler(mu_outer, r_outer, apo, peri, incline, Theta_res):
    
    cdef double torque_res = -1
    torque_res = gr_wrapper.J_dot_phi_Kepler(<double> mu_outer, <double> r_outer, <double> apo, <double> peri, <double> incline, <double> Theta_res)
    return torque_res













#All the function wrappers from Gamma.c




#Inputs: an int (nl), an int (n_res_inner), an int (k_res_inner), an int (m_res_inner), an int (n_res_outer), an int (k_res_outer), a double (ra_inner), a double (rp_inner), a double (I_inner), a double (ra_outer), a double (rp_outer), a double (I_outer), a double (astar), a double (M), a double (radius_outer), a double (delta_t)
#Outputs: an array of 2 doubles (Gamma)
#int omega_dot(int nl, int n_res_inner, int k_res_inner, int m_res_inner, int n_res_outer, int k_res_outer, double ra_inner, double rp_inner, double I_inner, double ra_outer, double rp_outer, double I_outer, double astar, double M, double radius_outer, double delta_t, double *Gamma)
def omegadot(nl, n_res_inner, k_res_inner, m_res_inner, n_res_outer, k_res_outer, ra_inner, rp_inner, I_inner, ra_outer, rp_outer, I_outer, astar, M, radius_outer, delta_t):

    cdef double *Gamma = <double *> malloc (2 * sizeof(double))
    cdef int result = -1
    try:
        result = gr_wrapper.omega_dot(<int> nl, <int> n_res_inner, <int> k_res_inner, <int> m_res_inner, <int> n_res_outer, <int> k_res_outer, <double> ra_inner, <double> rp_inner, <double> I_inner, <double> ra_outer, <double> rp_outer, <double> I_outer, <double> astar, <double> M, <double> radius_outer, <double> delta_t, Gamma)
    finally:
        output = [Gamma[i] for i in range(2)]
        free(Gamma)

    return result, output

    








#All the function wrappers from Delta_J.c


#double Delta_J_tidal(int nl, int n_res_inner, int n_res_outer, int k_res_inner, int k_res_outer, int m_res_inner, int m_res_outer, double apo, double rp, double radius_outer, double I, double M, double astar, double theta_res_F, double *Delta_J_r_tidal, double *Delta_J_theta_tidal, double *Delta_J_phi_tidal)
#def delta_j_tidal(nl, n_res_inner, n_res_outer, k_res_inner, k_res_outer, m_res_inner, m_res_outer, apo, rp, radius_outer, I, M, astar, theta_res_F):

#    cdef double *Delta_J_r_tidal = <double *> malloc (2 * sizeof(double))
#    cdef double *Delta_J_theta_tidal = <double *> malloc (2 * sizeof(double))
#    cdef double *Delta_J_phi_tidal = <double *> malloc (2 * sizeof(double))
#    cdef double djt = -1
#    try:
#        djt = gr_wrapper.Delta_J_tidal(<int> nl, <int> n_res_inner, <int> n_res_outer, <int> k_res_inner, <int> k_res_outer, <int> m_res_inner, <int> m_res_outer, <double> apo, <double> rp, <double> radius_outer, <double> I, <double> M, <double> astar, <double> theta_res_F, Delta_J_r_tidal, Delta_J_theta_tidal, Delta_J_phi_tidal)
#    finally:
#        output1 = [Delta_J_r_tidal[i] for i in range(2)]
#        output2 = [Delta_J_theta_tidal[i] for i in range(2)]
#        output3 = [Delta_J_phi_tidal[i] for i in range(2)]
#        free(Delta_J_r_tidal)
#        free(Delta_J_theta_tidal)
#        free(Delta_J_phi_tidal)

#    return djt, output1, output2, output3
