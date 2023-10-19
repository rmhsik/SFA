#include "propag/harmonic_pysfa.h"
#include "propag/config.h"
#include "pysfa.h"
#include <algorithm>
#include <fstream>
#include <chrono>
#include <filesystem>

void createDirectories_pysfa(std::string name){
    try{
        bool status = std::filesystem::create_directories(name);
        if (!status)
            throw(status);
    }
    catch(bool status){
        SPDLOG_ERROR("[HarmonicTDSESolver] Temporal directories could not be created.");
    }
}


void writeField(double *field_x, double *field_y, double *field_z, const int nt, std::string tmppath){
    std::string path_field_x = tmppath + "/field_x.dat";
    std::string path_field_y = tmppath + "/field_y.dat";
    std::string path_field_z = tmppath + "/field_z.dat";
    std::ofstream ffieldx(path_field_x);
    if(ffieldx.is_open()){
        for(int i=0; i<nt; i++){
            ffieldx<<std::setprecision(12)<<std::fixed<<field_x[i]<<std::endl;
        }
        ffieldx.close();
    }
    std::ofstream ffieldy(path_field_y);
    if(ffieldy.is_open()){
        for(int i=0; i<nt; i++){
            ffieldy<<std::setprecision(12)<<std::fixed<<field_y[i]<<std::endl;
        }
        ffieldy.close();
    }
    std::ofstream ffieldz(path_field_z);
    if(ffieldz.is_open()){
        for(int i=0; i<nt; i++){
            ffieldz<<std::setprecision(12)<<std::fixed<<field_z[i]<<std::endl;
        }
        ffieldz.close();
    }
}


void writeDipole(double *dipole_x, double *dipole_y, double *dipole_z, const int nt, std::string tmppath){
    std::string path_dipole_x = tmppath + "/dipole_x.dat";
    std::string path_dipole_y = tmppath + "/dipole_y.dat";
    std::string path_dipole_z = tmppath + "/dipole_z.dat";
    std::ofstream fdipolex(path_dipole_x);
    if(fdipolex.is_open()){
        for(int i=0; i<nt; i++){
            fdipolex<<std::setprecision(12)<<std::fixed<<dipole_x[i]<<std::endl;
        }
        fdipolex.close();
    }
    std::ofstream fdipoley(path_dipole_y);
    if(fdipoley.is_open()){
        for(int i=0; i<nt; i++){
            fdipoley<<std::setprecision(12)<<std::fixed<<dipole_y[i]<<std::endl;
        }
        fdipoley.close();
    }
    std::ofstream fdipolez(path_dipole_z);
    if(fdipolez.is_open()){
        for(int i=0; i<nt; i++){
            fdipolez<<std::setprecision(12)<<std::fixed<<dipole_z[i]<<std::endl;
        }
        fdipolez.close();
    }
}


HarmonicPySFA::HarmonicPySFA(std::shared_ptr<Settings> settings) : Harmonic(settings){
    Ip = mSettings->mData.pysfa_Ip;
    Z  = mSettings->mData.pysfa_Z;
    n_prin = mSettings->mData.pysfa_nprin;
    nthreads = mSettings->mData.numThreads;
    nt = mSettings->NDT;
    dt = mSettings->DT;
    t = new double[nt];
    createTemporalArray();
}

HarmonicPySFA::~HarmonicPySFA(){
    delete [] t;
}

void HarmonicPySFA::calculateAcceleration(Field *field){
    std::string timestamp = std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(
                                          std::chrono::system_clock::now().time_since_epoch()).count());
    std::string tmppath = "tmp/results_"+timestamp;
    createDirectories_pysfa(tmppath);
    double *dipole_x, *dipole_y, *dipole_z;
    double *efieldx, *efieldy, *efieldz;
    dipole_x = new double[nt];
    dipole_y = new double[nt];
    dipole_z = new double[nt];
    efieldx  = new double[nt];
    efieldy  = new double[nt];
    efieldz  = new double[nt];
    for(int i=0; i<nt; i++){
        dipole_x[i] = 0.0;
        dipole_y[i] = 0.0;
        dipole_z[i] = 0.0;
        efieldx[i]  = field->mFieldX[i].real();
        efieldy[i]  = field->mFieldY[i].real();
        efieldz[i]  = field->mFieldZ[i].real();
    } 

    pySFA3D(Ip, Z, n_prin, efieldx, 
                           efieldy, 
                           efieldz, 
                           t, nt, nthreads,
                           dipole_x,
                           dipole_y,
                           dipole_z);
    writeField(efieldx, efieldy, efieldz, nt, tmppath);
    writeDipole(dipole_x, dipole_y, dipole_z, nt, tmppath);
    std::copy(dipole_x, dipole_x+nt, mAccelX);
    std::copy(dipole_y, dipole_y+nt, mAccelY);
    std::copy(dipole_z, dipole_z+nt, mAccelZ);

    delete [] dipole_x;
    delete [] dipole_y;
    delete [] dipole_z;
    delete [] efieldx;
    delete [] efieldy;
    delete [] efieldz;
}

void HarmonicPySFA::createTemporalArray(){
    for(int i=0; i<nt; i++){
        t[i] = i*dt;
    }
}
