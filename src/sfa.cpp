#include "sfa.h"
#include <iostream>
#include <cmath>
#include <omp.h>

extern "C"{
    void SFA(double Ip, double *efield, double *t, const int nt, const int nthreads, cdouble *dipole){
        double dt = t[1] -t[0];
        double *afield;
        afield = new double[nt];
        for(int i=0; i<nt;i++){
            afield[i] = 0.0;
        }
        omp_set_num_threads(nthreads);
        calc_vecpot(efield, t, nt, afield);
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<nt; i++){
            cdouble dipole_i = 0.0;
            cdouble tmp_value = 0.0;
            double action   = 0.0;
            for(int j=0; j<i; j++){
                tmp_value = 0.0;
                action = calc_action(afield, Ip, j, i, t, nt);
                tmp_value  = pow(M_PI/( 0.1 + I*(t[i]-t[j]) ), 1.5);
                tmp_value *= efield[j];
                tmp_value *= exp(-I*action);
                dipole_i += tmp_value;
            }
            dipole_i *= I*dt;
            dipole[i] = dipole_i;
        }
        delete [] afield;
    }

    void calc_vecpot(double *efield, double *t, const int nt, double *afield){
        double dt = t[1]-t[0];
        for(int i=0; i<nt; i++){
            for(int j=0;j<i;j++){
                afield[i] += -efield[j];
            }
            afield[i] *=dt;
        }
    }

    double calc_momentum(double *afield, const int idx_ion, const int idx_rec, double *t, const int nt){
        double dt = t[1]-t[0];
        double canonical_momentum = 0.0;
        double tau_ion = t[idx_ion];
        double tau_rec = t[idx_rec];
        double prefactor = -1.0/(tau_rec-tau_ion+0.001);

        for(int i=idx_ion; i<idx_rec;i++){
            canonical_momentum += afield[i];
        }
        canonical_momentum *= prefactor*dt;
        return canonical_momentum;
    }

    double calc_action(double *afield, double Ip, const int idx_ion, const int idx_rec, double *t, const int nt){
        double dt = t[1]-t[0];
        double action = 0.0;
        double momentum = calc_momentum(afield, idx_ion, idx_rec, t, nt);
        for(int i=idx_ion; i<idx_rec; i++){
            action += 0.5*(momentum + afield[i])*(momentum + afield[i]) +Ip;
        }
        action *= dt;
        return action;
    }
}
