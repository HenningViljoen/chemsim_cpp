#include "global.h"

int calcSimIterations()
{
    double delta = (global::SampleT == 0) ? 0.01 : 0;
    return (int)(global::SimTime / (global::SampleT + delta)); //Nr of iterations of the simulation.
}

int calcSimVectorLength()
{
    double delta = (global::SampleT == 0) ? 0.01 : 0;
    return (int)(global::SimTime / (global::SimVectorUpdateT + delta)); //Nr of iterations of the simulation.
}

double *initsimtimevector(double *simtimevector)
{
    if (simtimevector != nullptr)
    {
        delete simtimevector;
    }
    else
    {
        simtimevector = new double[global::SimVectorLength];
    }

    simtimevector[0] = 0;
    for (int i = 1; i < global::SimVectorLength; i++)
    {
        simtimevector[i] = simtimevector[i - 1] + global::SimVectorUpdateT;
    }
    return simtimevector;
}


int global::TimerInterval = 1; // micro seconds
double global::SpeedUpFactor = 200; //50, CT and normal sim: 200; //factor  Heat exchangers: 30000
double global::SampleT = global::TimerInterval / 1000.0 * global::SpeedUpFactor; // seconds - at this point SampleT
double global::SimVectorUpdateT = 1.0;
int global::SimVectorUpdatePeriod = (int)(round(global::SimVectorUpdateT / global::SampleT)); //Nr. samples between saving in vect.
double global::TrendUpdateT = 1.0; //seconds.  The simulation period of updating trends in simulation.
                                                 //; HX: 30s; normal sim with cooling tower: 1s.
int global::TrendUpdateIterPeriod =
        (int)(round(global::TrendUpdateT / global::SampleT)); //Nr. of samples between update trend.;
double global::SimTime = 3600.0*4; // Valve stepping: 3600.0*1; Normal full model: 3600.0*4 ; CT alone:  64830; //seconds 3600*1;//3600*24;  135: ES004 53340 ; for 172EP004: 34440;for CT fitting: 12 hours/3hours.
int global::SimIterations = calcSimIterations(); //Nr of iterations of the simulation.
int global::SimVectorLength = calcSimVectorLength();
double global::PlantDataSampleT = 30.0; //30 seconds for CT;  The sample frequency of the IP.21 data used for fitting.
double *global::simtimevector = initsimtimevector(global::simtimevector);

global::global()
{

}

global::~global()
{
    if (simtimevector != nullptr)
    {
        delete simtimevector;
    }
}


double global::calcSampleT()
{
    return TimerInterval / 1000.0 * SpeedUpFactor; // seconds - at this point SampleT
}







