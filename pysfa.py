import numpy as np
import ctypes as ct
import os

def pySFA(Ip, Z, n_prin, efield, t, nthreads=1):
    nt = t.shape[0]
    if efield.shape[0] != nt:
        raise Exception("efield array length must be equal to t array length")
        
    sfa_c = ct.CDLL(os.getcwd() + "/pySFA/build/libsfa.so")
    ND_POINTER_1_double = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C") 
    ND_POINTER_1_complex = np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1, flags="C")

    sfa_c.SFA.argtypes = [ct.c_double, ct.c_double, ct.c_double, ND_POINTER_1_double, ND_POINTER_1_double, ct.c_int32, ct.c_int32, ND_POINTER_1_double]
    sfa_c.SFA.restype = ct.c_void_p
    dipole_array = np.zeros(efield.shape, dtype=np.float64)
    sfa_c.SFA(Ip, Z, n_prin, efield.astype(np.float64), t.astype(np.float64), nt, nthreads, dipole_array)
    del sfa_c
    return dipole_array

def pySFA3D(Ip, Z, n_prin, efield, t, nthreads=1):
    nt = t.shape[0]
    if efield.shape[0] != 3:
        raise Exception("efield axis 0 length must be 3")
    if efield.shape[1] != nt:
        raise Exception("efield axis 1 length must be equal to t array length")
        
    sfa_c = ct.CDLL(os.getcwd() + "/pySFA/build/libsfa.so")
    ND_POINTER_1_double = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C") 
    ND_POINTER_1_complex = np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1, flags="C")

    sfa_c.SFA3D.argtypes = [ct.c_double, ct.c_double, ct.c_double, ND_POINTER_1_double, ND_POINTER_1_double, ND_POINTER_1_double, ND_POINTER_1_double, ct.c_int32, ct.c_int32, ND_POINTER_1_double, ND_POINTER_1_double, ND_POINTER_1_double]
    sfa_c.SFA3D.restype = ct.c_void_p
    dipole_x = np.zeros(efield[0].shape, dtype=np.float64)
    dipole_y = np.zeros(efield[1].shape, dtype=np.float64)
    dipole_z = np.zeros(efield[2].shape, dtype=np.float64)
    sfa_c.SFA3D(Ip, Z, n_prin, efield[0].astype(np.float64), efield[1].astype(np.float64), efield[2].astype(np.float64),
             t.astype(np.float64), nt, nthreads,
             dipole_x, dipole_y, dipole_z)

    del sfa_c
    return np.array([dipole_x, dipole_y, dipole_z])
