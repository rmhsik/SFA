import numpy as np
import ctypes as ct

def SFA(Ip, efield, t, nt, nthreads=1):
    sfa_c = ct.CDLL("build/libsfa.so")
    ND_POINTER_1_double = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C") 
    ND_POINTER_1_complex = np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1, flags="C")

    sfa_c.SFA.argtypes = [ct.c_double, ND_POINTER_1_double, ND_POINTER_1_double, ct.c_int32, ct.c_int32, ND_POINTER_1_complex]
    sfa_c.SFA.restype = ct.c_void_p
    dipole_array = np.zeros(efield.shape, dtype=np.complex128)
    sfa_c.SFA(Ip, efield.astype(np.float64), t.astype(np.float64), nt, nthreads, dipole_array)
    del sfa_c
    return dipole_array
