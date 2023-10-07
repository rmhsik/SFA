#ifndef SFA_H
#define SFA_H

#include <complex>
#define cdouble std::complex<double>
#define I std::complex<double>(0.0,1.0)

extern "C" void SFA(double Ip, double *efield, double *t, const int nt, const int nthreads, cdouble *dipole);
extern "C" void calc_vecpot(double *efield, double *t, const int nt, double *afield);
extern "C" double calc_momentum(double *afield, const int idx_ion, const int idx_rec, double *t, const int nt);
extern "C" double calc_action(double *afield, double Ip, const int idx_ion, const int idx_rec, double *t, const int nt);

#endif
