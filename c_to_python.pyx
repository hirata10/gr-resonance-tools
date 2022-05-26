from libc.stdlib cimport malloc, free
cimport c_to_python

# int CKerr_J2EQL(double *J, double *EQL, double M, double astar);
def ckerr_j2eql(arrJ, M, astar):
    assert len(arrJ) == 3

    cdef double *arrJ_cpy = <double *> malloc(len(arrJ)*sizeof(double))
    for i in range(len(arrJ)):
        arrJ_cpy[i] = <double> arrJ[i]

    cdef int result = -1 
    cdef double *EQL = <double *> malloc (3 * sizeof(double))
    try:
        result = c_to_python.CKerr_J2EQL(arrJ_cpy, EQL, <double> M, <double> astar)
    finally:
        output = [EQL[i] for i in range(3)]
        free(EQL)
        free(arrJ_cpy)

    return result, output





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
        result = c_to_python.CKerr_EQL2J(arrEQL_cpy, J, <double> M, <double> astar, ancilliary)
    finally:
        output1 = [J[i] for i in range(3)]
        output2 = [ancilliary[i] for i in range(3)]
        free(J)
        free(ancilliary)
        free(arrEQL_cpy)

    return result, output1, output2






#int CKerr_Minverse(double *J, double *Minv, double M, double astar);
def ckerr_minverse(arrJ, M, astar):
    assert len(arrJ) == 3

    cdef double *arrJ_cpy = <double *> malloc(len(arrJ)*sizeof(double))
    for i in range(len(arrJ)):
        arrJ_cpy[i] = <double> arrJ[i]

    cdef int result = -1 
    cdef double *Minv = <double *> malloc (9 * sizeof(double))
    try:
        result = c_to_python.CKerr_Minverse(arrJ_cpy, Minv, <double> M, <double> astar)
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
        c_to_python.CKerr_Minv2Omega(arrMinv_cpy, Omega)
    finally:
        output = [Omega[i] for i in range(3)]
        free(Omega)
        free(arrMinv_cpy)

    return output










