#ifndef UTILS_H
#define UTILS_H

#include <complex>
#define cdouble std::complex<double>
#define I std::complex<double>(0.0,1.0)

double max_efield(double *efield, const int n);
void calc_vecpot(double *efield, double *t, const int nt, double *afield);
double calc_momentum(double *afield, const int idx_ion, const int idx_rec, double *t);
double calc_action(double *afield, double momentum, double Ip, const int idx_ion, const int idx_rec, double *t);
cdouble calc_matrix_element(double ps, double Z, double n, double E0, double Ip);
void calc_momentum3D(double *afield_x, double *afield_y, double *afield_z, const int idx_ion, const int idx_rec, double*t, double *momentum);
double calc_action3D(double *afield_x, double *afield_y, double *afield_z, double *momentum, double Ip, const int idx_ion, const int idx_rec, double *t);
void calc_matrix_element3D(double *ps, double afield_x, double afield_y, double afield_z, double Z, double n, double E0, double Ip, cdouble *matrix_element);
#endif
