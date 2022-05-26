cdef extern from "CKerr.h":
    int CKerr_J2EQL(double *J, double *EQL, double M, double astar)
    int CKerr_EQL2J(double *EQL, double *J, double M, double astar, double *ancillary)
    int CKerr_Minverse(double *J, double *Minv, double M, double astar)
    void CKerr_Minv2Omega(double *Minv, double *Omega)








