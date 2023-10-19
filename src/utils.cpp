#include "utils.h"
#include <cmath>
#include <iostream>

double max_efield(double *efield, const int nt){
    double largestE = 0.0;

    for(int i=0; i<nt; i++){
        if(std::abs(efield[i])>largestE){
            largestE = std::abs(efield[i]);
        }
    }
    return largestE;
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

double calc_momentum(double *afield, const int idx_ion, const int idx_rec, double *t){
    double dt = t[1]-t[0];
    double canonical_momentum = 0.0;//1E-6;
    double tau_ion = t[idx_ion];
    double tau_rec = t[idx_rec];
    double prefactor = -1.0/(tau_rec-tau_ion);

    for(int i=idx_ion; i<idx_rec;i++){
        canonical_momentum += afield[i];
    }
    canonical_momentum *= prefactor*dt;
    if(std::abs(canonical_momentum)<1E-5){
        canonical_momentum=1E-5;
    }   
    return canonical_momentum;
}

void calc_momentum_array(double *afield, double *momentum, const int idx_rec, double *t){
    double dt = t[1]-t[0];
    double tau_ion = 0.0;
    double tau_rec = t[idx_rec];
    double prefactor = 0.0;
    for(int i=0; i<idx_rec; i++){
        momentum[0] += afield[i];
    }
    if(idx_rec>0){
        for(int i=1; i<idx_rec; i++){
            momentum[i] = momentum[i-1]-afield[i-1]; 
        }
    }
    for(int i=0; i<idx_rec; i++){
        tau_ion = t[i];
        prefactor = -1.0/(tau_rec-tau_ion);
        momentum[i] *= prefactor*dt;  
    }
}

void calc_momentum3D(double *afield_x, double *afield_y, double* afield_z, const int idx_ion, const int idx_rec, double *t, double *momentum){
    double dt = t[1] - t[0];
    double canonical_momentum_x = 0.0;
    double canonical_momentum_y = 0.0;
    double canonical_momentum_z = 0.0;
    double tau_ion = t[idx_ion];
    double tau_rec = t[idx_rec];
    double prefactor = -1.0/(tau_rec-tau_ion+0.000001);
    for(int i=idx_ion; i<idx_rec; i++){
        canonical_momentum_x += afield_x[i];
        canonical_momentum_y += afield_y[i];
        canonical_momentum_z += afield_z[i];
    }
    momentum[0] = canonical_momentum_x*prefactor*dt;
    momentum[1] = canonical_momentum_y*prefactor*dt;
    momentum[2] = canonical_momentum_z*prefactor*dt;
    return;
}

void calc_momentum3D_array(double *afield_x, double *afield_y, double* afield_z, const int idx_rec, double *t, double *momentum){
    double dt = t[1] - t[0];
    double tau_ion = 0.0;
    double tau_rec = t[idx_rec];
    double prefactor = 0.0;
    for(int i=0; i<idx_rec; i++){
        momentum[0*3 + 0] += afield_x[i];
        momentum[0*3 + 1] += afield_y[i];
        momentum[0*3 + 2] += afield_z[i];
    }
    if(idx_rec>0){
        for(int i=1; i<idx_rec; i++){
            momentum[i*3 + 0] = momentum[(i-1)*3+0] - afield_x[i-1];
            momentum[i*3 + 1] = momentum[(i-1)*3+1] - afield_y[i-1];
            momentum[i*3 + 2] = momentum[(i-1)*3+2] - afield_z[i-1];
        }
    }
    for(int i=0; i<idx_rec; i++){
        tau_ion = t[i];
        prefactor = -1.0/(tau_rec-tau_ion);
        momentum[i*3 + 0] *= prefactor*dt;
        momentum[i*3 + 1] *= prefactor*dt;
        momentum[i*3 + 2] *= prefactor*dt;
    }
}

    double calc_action(double *afield, double momentum, double Ip, const int idx_ion, const int idx_rec, double *t){
        double dt = t[1]-t[0];
        double action = 0.0;
        // double momentum = calc_momentum(afield, idx_ion, idx_rec, t);
        for(int i=idx_ion; i<idx_rec; i++){
            action += 0.5*(momentum + afield[i])*(momentum + afield[i]) +Ip;
        }
        action *= dt;
        return action;
    }
void calc_action_array(double *afield, double *momentum, double *action, double Ip, const int idx_rec, double *t){
    double dt = t[1]-t[0];
    for(int i=0; i<idx_rec; i++){
        action[0] += (0.5*(momentum[i] + afield[i])*(momentum[i] + afield[i]) + Ip)*dt;
    }
    
    // double momentum = calc_momentum(afield, idx_ion, idx_rec, t);
    for(int i=1; i<idx_rec; i++){
        action[i] = action[i-1] - (0.5*(momentum[i-1] + afield[i-1])*(momentum[i-1] + afield[i-1]) +Ip)*dt;
    }
}   
    double calc_action3D(double *afield_x, double *afield_y, double *afield_z, double *momentum, double Ip, const int idx_ion, const int idx_rec, double *t){
        double dt = t[1] - t[0];
        double action = 0.0;
        //double *momentum = calc_momentum3D(afield_x, afield_y, afield_z, idx_ion, idx_rec, t);
        for (int i = idx_ion; i<idx_rec; i++){
            action += 0.5*(momentum[0]*momentum[0] + afield_x[i]*afield_x[i] + 2.0*momentum[0]*afield_x[i]);
            action += 0.5*(momentum[1]*momentum[1] + afield_y[i]*afield_y[i] + 2.0*momentum[1]*afield_y[i]);
            action += 0.5*(momentum[2]*momentum[2] + afield_z[i]*afield_z[i] + 2.0*momentum[2]*afield_z[i]);
            action += Ip;
        }
        action *= dt;
        return action;
    }
void calc_action3D_array(double *afield_x, double *afield_y, double *afield_z, double *momentum, double *action, double Ip, const int idx_rec, double *t){
        double dt = t[1] - t[0];
        //double *momentum = calc_momentum3D(afield_x, afield_y, afield_z, idx_ion, idx_rec, t);
        for (int i = 0; i<idx_rec; i++){
            action[0] += (0.5*(momentum[i*3+0]+afield_x[i])*(momentum[i*3+0]+afield_x[i]))*dt;
            action[0] += (0.5*(momentum[i*3+1]+afield_y[i])*(momentum[i*3+1]+afield_y[i]))*dt;
            action[0] += (0.5*(momentum[i*3+2]+afield_z[i])*(momentum[i*3+2]+afield_z[i]) +Ip)*dt;
        }
        for(int i=1; i<idx_rec; i++){
            action[i] = action[i-1] 
                        -(0.5*(momentum[(i-1)*3+0]+afield_x[i-1])*(momentum[(i-1)*3+0]+afield_x[i-1]))*dt
                        -(0.5*(momentum[(i-1)*3+1]+afield_y[i-1])*(momentum[(i-1)*3+1]+afield_y[i-1]))*dt
                        -(0.5*(momentum[(i-1)*3+2]+afield_z[i-1])*(momentum[(i-1)*3+2]+afield_z[i-1])+Ip)*dt;
        }
    }

    cdouble calc_matrix_element(double ps, double Z, double n, double E0, double Ip){
        double a0 = 1.0; double q0 = 1.0; 
        double ps_mod = std::sqrt(ps*ps);
        double ps_mod2 = ps_mod*ps_mod;
        double atan_term = 1.0;//atan2(ps_mod*a0,Z);
        double Z_a0 = Z/a0;
        cdouble tmp_value = 0.0;
        if (ps_mod<1E-5){
            return tmp_value;
        }
        tmp_value = I*0.450158158*pow(Z_a0,1.5);
        tmp_value *= ps/pow(ps_mod,3)*pow(ps_mod2 + (Z_a0)*(Z_a0),-(3.0-n)/2.0);
        tmp_value *= 1.0;//tgamma(2.0-n);
        tmp_value *= ps_mod*(2-n)*cos((3-n)*atan_term) 
                     - sqrt(ps_mod2+(Z_a0)*(Z_a0))*sin((2.0-n)*atan_term);
        //tmp_value = -I/M_PI*pow(2/a0,1.5)*ps/(pow(ps_mod*ps_mod+(1.0/a0)*(1/a0),0));
        //tmp_value *= pow(4.0*Ip/(q0*E0),n);
        return tmp_value;
    }

    void calc_matrix_element_array(double *momentum, double *afield, cdouble *matrix_element_ion, cdouble *matrix_element_rec, const int idx_rec, double Z, double n, double E0, double Ip){
        const double a0 = 1.0; const double q0 = 1.0; 
        const double Z_a0 = Z/a0;
        const double gamma_term = std::tgamma(2.0-n);
        const double sqrt2_pi = std::sqrt(2.0)/M_PI;
        const cdouble prefactor = I*sqrt2_pi*std::pow(Z_a0,1.5);
        double ps = 0.0;
        double ps_mod = 0.0;
        double ps_mod2 = 0.0;
        double atan_term = 0.0;
        for (int i=0; i<idx_rec; i++){
            ps = momentum[i] + afield[i];
            ps_mod = std::sqrt(ps*ps);
            ps_mod2 = ps_mod*ps_mod;
            if(ps_mod<1E-5){
                matrix_element_ion[i] = 0.0;
                continue;
            }
            else{
                atan_term = atan2(ps_mod*a0,Z);
                matrix_element_ion[i]  = prefactor * ps/(ps_mod*ps_mod*ps_mod);
                matrix_element_ion[i] *= pow(ps_mod2 + Z_a0*Z_a0, -(3.0-n)/2.0);
                matrix_element_ion[i] *= gamma_term;
                matrix_element_ion[i] *= ps_mod*(2.0-n)*cos((3.0-n)*atan_term)
                                        -sqrt(ps_mod2+Z_a0*Z_a0)*sin((2.0-n)*atan_term);
            }
            ps = momentum[i] + afield[idx_rec];
            ps_mod = std::sqrt(ps*ps);
            ps_mod2 = ps_mod*ps_mod;
            if(ps_mod<1E-5){
                matrix_element_rec[i] = 0.0;
            }
            else{
                atan_term = atan2(ps_mod*a0,Z);
                matrix_element_rec[i]  = prefactor * ps/(ps_mod*ps_mod*ps_mod);
                matrix_element_rec[i] *= pow(ps_mod2 + Z_a0*Z_a0, -(3.0-n)/2.0);
                matrix_element_rec[i] *= gamma_term;
                matrix_element_rec[i] *= ps_mod*(2.0-n)*cos((3.0-n)*atan_term)
                                        -sqrt(ps_mod2+Z_a0*Z_a0)*sin((2.0-n)*atan_term);
            }

        }
    }



    void calc_matrix_element3D(double *ps, double afield_x, double afield_y, double afield_z, double Z, double n, double E0, double Ip, cdouble *matrix_element){
        double a0 = 1.0; double q0 = 1.0;
        double ps_mod = sqrt(  (ps[0]+afield_x)*(ps[0]+afield_x) 
                             + (ps[1]+afield_y)*(ps[1]+afield_y) 
                             + (ps[2]+afield_z)*(ps[2]+afield_z));
        matrix_element[0] = 0.0;
        matrix_element[1] = 0.0;
        matrix_element[2] = 0.0;
        if(ps_mod < 1E-5){
            return;    
        }
        matrix_element[0] = I*sqrt(2.0)/M_PI*pow(Z/a0,1.5);
        matrix_element[0] *= (ps[0]+afield_x)/pow(ps_mod,3)*pow(ps_mod*ps_mod + (Z/a0)*(Z/a0),-(3.0-n)/2.0);
        matrix_element[0] *= tgamma(2.0-n);
        matrix_element[0] *= ps_mod*(2-n)*cos((3-n)*atan(ps_mod*a0/Z)) 
                     - sqrt(ps_mod*ps_mod+(Z/a0)*(Z/a0))*sin((2.0-n)*atan(ps_mod*a0/Z));
        matrix_element[0] *= pow(4.0*Ip/(q0*E0),n);

        matrix_element[1] = I*sqrt(2.0)/M_PI*pow(Z/a0,1.5);
        matrix_element[1] *= (ps[1]+afield_y)/pow(ps_mod,3)*pow(ps_mod*ps_mod + (Z/a0)*(Z/a0),-(3.0-n)/2.0);
        matrix_element[1] *= tgamma(2.0-n);
        matrix_element[1] *= ps_mod*(2-n)*cos((3-n)*atan(ps_mod*a0/Z)) 
                     - sqrt(ps_mod*ps_mod+(Z/a0)*(Z/a0))*sin((2.0-n)*atan(ps_mod*a0/Z));
        matrix_element[1] *= pow(4.0*Ip/(q0*E0),n);

        matrix_element[2] = I*sqrt(2.0)/M_PI*pow(Z/a0,1.5);
        matrix_element[2] *= (ps[2]+afield_z)/pow(ps_mod,3)*pow(ps_mod*ps_mod + (Z/a0)*(Z/a0),-(3.0-n)/2.0);
        matrix_element[2] *= tgamma(2.0-n);
        matrix_element[2] *= ps_mod*(2-n)*cos((3-n)*atan(ps_mod*a0/Z)) 
                     - sqrt(ps_mod*ps_mod+(Z/a0)*(Z/a0))*sin((2.0-n)*atan(ps_mod*a0/Z));
        matrix_element[2] *= pow(4.0*Ip/(q0*E0),n);
        return;
    }

