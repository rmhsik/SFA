#include "pysfa.h"
#include "utils.h"
#include <iostream>
#include <cmath>
#include <omp.h>

cdouble fastPow1_5(cdouble a){
    cdouble sqrt_a = sqrt(a);
    return sqrt_a*a;
}

extern "C"{
    void pySFA(double Ip, double Z, double n_prin, double *efield, double *t, const int nt, const int nthreads, double *dipole){
        double dt = t[1] -t[0];
        double *afield;
        double maxE0 = 0.0;
        afield = new double[nt];
        for(int i=0; i<nt;i++){
            afield[i] = 0.0;
        }
        omp_set_num_threads(nthreads);
        maxE0 = max_efield(efield, nt);
        calc_vecpot(efield, t, nt, afield);
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<nt; i++){
            cdouble dipole_i = 0.0;
            cdouble tmp_value = 0.0;
            double action   = 0.0;
            //double momentum = 0.0;
            double *momentum_array;
            double *action_array;
            cdouble *matrix_dipole_ion;
            cdouble *matrix_dipole_rec;
            momentum_array = new double[i];
            action_array = new double[i];
            matrix_dipole_ion = new cdouble[i];
            matrix_dipole_rec = new cdouble[i];
            for(int j=0; j<i; j++){
                momentum_array[j] = 0.0;
                action_array[j] = 0.0;
                matrix_dipole_ion[j] = cdouble(0.0,0.0);
                matrix_dipole_rec[j] = cdouble(0.0,0.0);
            }
            calc_momentum_array(afield, momentum_array, i, t);
            calc_action_array(afield, momentum_array, action_array, Ip, i, t);
            calc_matrix_element_array(momentum_array, afield, matrix_dipole_ion,
                                      matrix_dipole_rec, i, Z, n_prin, maxE0, Ip);
            for(int j=0; j<i; j++){
                tmp_value = 0.0;
                action = action_array[j];
                tmp_value = fastPow1_5(2.0*M_PI/( 0.00001 + I*(t[i]-t[j])) );
                tmp_value *= efield[j];
                tmp_value *= exp(-I*action);
                dipole_i += tmp_value*matrix_dipole_ion[j]*std::conj(matrix_dipole_rec[j]);
            }
            delete [] momentum_array;
            delete [] action_array;
            delete [] matrix_dipole_ion;
            delete [] matrix_dipole_rec;
            dipole_i *= I*dt;
            dipole[i] = dipole_i.real();
        }
        delete [] afield;
    }

    void pySFA3D(double Ip, double Z, double n_prin, double *efield_x, double *efield_y, double *efield_z, double *t, const int nt, const int nthreads, double *dipole_x, double *dipole_y, double *dipole_z){
        double dt = t[1] - t[0];
        double *afield_x;
        double *afield_y;
        double *afield_z;
        double maxE0 = 0.0;
        double maxE0x = 0.0;
        double maxE0y = 0.0;
        double maxE0z = 0.0;
        afield_x = new double[nt];
        afield_y = new double[nt];
        afield_z = new double[nt];
        for(int i=0; i<nt; i++){
            afield_x[i] = 0.0;
            afield_y[i] = 0.0;
            afield_z[i] = 0.0;
        }
        maxE0x = max_efield(efield_x, nt);
        maxE0y = max_efield(efield_y, nt);
        maxE0z = max_efield(efield_z, nt);
        maxE0 = maxE0z;
        if (maxE0x>maxE0){
            maxE0 = maxE0x;
        }
        if (maxE0y>maxE0){
            maxE0 = maxE0y;
        }
        omp_set_num_threads(nthreads);
        calc_vecpot(efield_x, t, nt, afield_x);
        calc_vecpot(efield_y, t, nt, afield_y);
        calc_vecpot(efield_z, t, nt, afield_z);
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<nt; i++){
            cdouble tmp_dipole_x = 0.0;
            cdouble tmp_dipole_y = 0.0;
            cdouble tmp_dipole_z = 0.0;
            cdouble tmp_value_x = 0.0;
            cdouble tmp_value_y = 0.0;
            cdouble tmp_value_z = 0.0;
            cdouble prefactor = 0.0;
            double action = 0.0;
            double *momentum = new double[3];
            double *momentum_array = new double[3*i];
            double *action_array = new double[i];
            cdouble *matrix_dipole_ion = new cdouble[3];
            cdouble *matrix_dipole_rec = new cdouble[3];
            for(int j = 0; j<i; j++){
                action_array[j] = 0.0;
                momentum_array[j*3 + 0] = 0.0;
                momentum_array[j*3 + 1] = 0.0;
                momentum_array[j*3 + 2] = 0.0;
            }
            calc_momentum3D_array(afield_x, afield_y, afield_z, i, t, momentum_array);
            calc_action3D_array(afield_x, afield_y, afield_z, momentum_array, action_array, Ip, i, t);
            for(int j=0; j<i; j++){
                tmp_value_x = 0.0;
                tmp_value_y = 0.0;
                tmp_value_z = 0.0;
                std::copy(&momentum_array[j*3], &momentum_array[j*3+3], momentum);
                action = action_array[j];
                calc_matrix_element3D(momentum, afield_x[j], 
                                                afield_y[j], 
                                                afield_z[j], Z, n_prin, maxE0, Ip, matrix_dipole_ion);
                calc_matrix_element3D(momentum, afield_x[i],
                                                afield_y[i],
                                                afield_z[i], Z, n_prin, maxE0, Ip, matrix_dipole_rec);
                prefactor = fastPow1_5(2.0*M_PI/(0.00001 + I*(t[i]-t[j])));
                tmp_value_x = prefactor;
                tmp_value_y = prefactor;
                tmp_value_z = prefactor;
                tmp_value_x *= exp(-I*action);
                tmp_value_y *= exp(-I*action);
                tmp_value_z *= exp(-I*action);
                tmp_dipole_x += tmp_value_x*(matrix_dipole_ion[0]*efield_x[j] +
                                             matrix_dipole_ion[1]*efield_y[j] +
                                             matrix_dipole_ion[2]*efield_z[j]) *
                                             std::conj(matrix_dipole_rec[0]);
                tmp_dipole_y += tmp_value_y*(matrix_dipole_ion[0]*efield_x[j] +
                                             matrix_dipole_ion[1]*efield_y[j] +
                                             matrix_dipole_ion[2]*efield_z[j]) *
                                             std::conj(matrix_dipole_rec[1]);
                tmp_dipole_z += tmp_value_z*(matrix_dipole_ion[0]*efield_x[j] +
                                             matrix_dipole_ion[1]*efield_y[j] +
                                             matrix_dipole_ion[2]*efield_z[j]) *
                                             std::conj(matrix_dipole_rec[2]);
            }
            delete[] momentum;
            delete[] momentum_array;
            delete[] action_array;
            delete[] matrix_dipole_ion;
            delete[] matrix_dipole_rec;
            tmp_dipole_x *= I*dt;
            tmp_dipole_y *= I*dt;
            tmp_dipole_z *= I*dt;
            dipole_x[i] = tmp_dipole_x.real();
            dipole_y[i] = tmp_dipole_y.real();
            dipole_z[i] = tmp_dipole_z.real();
        }
        delete[] afield_x;
        delete[] afield_y;
        delete[] afield_z;
    }
}
