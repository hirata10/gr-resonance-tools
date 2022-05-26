from libc.stdlib cimport malloc, free
cimport gr_wrapper



# All the function wrappers from kerrtraj.c

#Check!
# int CKerr_J2EQL(double *J, double *EQL, double M, double astar);
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




#Check!
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





#Check!
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






#Check everything below here!


#Check!
#double CKerr_Emin(double Q, double L, double M, double astar)
def ckerr_emin(Q, L, M, astar):

    cdef double result = -1
    result = gr_wrapper.CKerr_Emin(<double> Q, <double> L, <double> M, <double> astar)
    return result



#Check!
#double CKerr_Emax(double Q, double L, double M, double astar, double *Jrmax)
def ckerr_emax(Q, L, M, astar, arrJmax):
    assert len(arrJmax) == 1

    cdef double *arrJmax_cpy = <double *> malloc(len(arrJmax)*sizeof(double))
    for i in range(len(arrJmax)):
        arrJmax_cpy[i] = <double> arrJmax[i]

    cdef double result = -1
    try:
        result = gr_wrapper.CKerr_Emax(<double> Q, <double> L, <double> M, <double> astar, arrJmax_cpy)
    finally:
        free(arrJmax_cpy)

    return result


#Check!
#double CKerr_QLJr2E(double Q, double L, double Jr, double M, double astar, double *Jtheta)
def ckerr_qljr2e(Q, L, Jr, M, astar):

    cdef double *Jtheta = <double *> malloc (3 * sizeof(double))

    cdef double result = -1
    try:
        result = gr_wrapper.CKerr_QLJr2E(<double> Q, <double> L, <double> Jr, <double> M, <double> astar, Jtheta)
    finally:
        output = [Jtheta[i] for i in range(1)]
        free(Jtheta)

    return result, output


#Check!
#Computes the derivatives of the energy and Carter constant with respect to each of the actions.  Returns 1 if successful, 0 if no such orbit.
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
        #Note how L does not get one!
        output1 = [dEdJ[i] for i in range(3)]
        output2 = [dQdJ[i] for i in range(3)]
        free(arrJ_cpy)
        free(arrEQL_cpy)
        free(dEdJ)
        free(dQdJ)

    return result, output1, output2
    
#Check!
#/* Computes the M-inverse matrix and returns it in row-order (0,0) (0,1)
# * (0,2) ... (2,2).  Uses 9 slots.  This routine uses the numerical
# * differentiation method.
# * Returns 1 if successful, 0 if no such orbit.
# */
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


#Check!
#/* Returns the phase-space coordinates (r,theta,phi,u_r,u_theta,u_phi) of the point
# * on the torus with actions J and at zero angles psi=(0,0,0).
# * Returns 1 if successful and 0 if the torus does not exist.
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


#Check!
#/* Returns the flow vector for the Hamiltonian and Carter constant, I dH and I dQ,
# * in 6D phase space.
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


#Check!
#/* Returns the flow vector in the psi_r and psi_theta directions (I dJr and I dJtheta)
# * in 6D phase space.  Requires the M-inverse matrix Minv[9], which can be obtained using
# * CKerr_Minverse, in order to convert H and Q flow-vectors into angle flows.
# */
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






#Check!
# Returns the angular momentum of a circular equatorial orbit with the specified
# radius.  Formulas from Chandrasekhar ch.7.
# double CKerr_getL_CircEq(double M, double astar, double r) 
# Returns dynamical data on circular equatorial orbit with the specified
# radius.  Formulas from Chandrasekhar ch.7 and Paper I.
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


#Check!
#/* Returns the location of the order-mm Lindblad resonance around the
# * specified black hole.  positive mm finds ILRs; negative mm finds the
# * OLR with order |mm|.
# *
# * Only integer mm gives Lindblad resonances, but we preserve the option to
# * insert fractional mm to search for higher-order eccentricity resonances --
# * e.g. the 3:1 resonance involving two factors of the test particle eccentricity
# * can be found by setting mm = 1.5.
# *
# * However may be ill-behaved if |mm|<1.
# */
#double CKerr_FindLindbladResonance(double M, double astar, double r0, double mm)
def ckerr_findlindbladresonance(M, astar, r0, mm):

    cdef double result = -1
    result = gr_wrapper.CKerr_FindLindbladResonance(<double> M, <double> astar, <double> r0, <double> mm)
    return result



#Check!
# Obtains constants (E,Q,L) for a trajectory of specified Carter inclination,
# pericenter, and apocenter.  Returns 1 if successful, 0 if no solution.
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


#Check!
# /* Obtains constants (E,Q,L) for a trajectory of specified Carter inclination,
# * and a circular orbit.  Returns 1 if successful, 0 if no solution.
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


#Check!
#/* Finds the orbit at which the ratio of azimuthal to vertical frequencies is fratio,
# * for (Carter) inclination Ic.  Also returns the actions J[3].
# * Currently requires fratio >= 2, astar > 0.
# * Returns 1 if the orbit exists, 0 otherwise.
# int CKerr_FindResCirc(int fratio, double Ic, double *J, double M, double astar)
# Giving issues? Converting double pointer to python object
def ckerr_findrescirc(fratio, Ic, M, astar):
    cdef double *J = <double *> malloc (3 * sizeof(double))
    cdef int result = -1
    try:
        result = gr_wrapper.CKerr_FindResCirc(<int> fratio, <double> Ic, J, <double> M, <double> astar)
    finally:
        output = [J[i] for i in range(3)]
        free(J)
        
    return result, output

#Check!
#Finds the ISCO radius for a general inclination orbit, as a function of (Carter) inclination Ic.
#double CKerr_FindISCO(double Ic, double M, double astar)
def ckerr_findisco(Ic, M, astar):

    cdef double result = -1
    result = gr_wrapper.CKerr_FindISCO(<double> Ic, <double> M, <double> astar)
    return result






