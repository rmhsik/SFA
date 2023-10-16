#ifndef HARMONIC_PYSFA_H
#define HARMONIC_PYSFA_H

#include "propag/harmonic.h"
#include "pysfa.h"

class HarmonicPySFA : public Harmonic{
    public:
        HarmonicPySFA(std::shared_ptr<Settings> settings);
        virtual ~HarmonicPySFA();
        virtual void calculateAcceleration(Field *field);

    protected:
        double Ip;
        double Z;
        double n_prin;
        double *t;
        double dt;
        int nt;
        int nthreads;

    private:
        void createTemporalArray();
};
#endif 
