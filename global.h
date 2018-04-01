#ifndef GLOBAL_H
#define GLOBAL_H

#include <math.h>

enum DrawModeEntity {EditMode, ValveMode};

enum objecttypes { FTReactor, GasPipe, LiquidPipe, Pump, Tank, Valve, Tee, Mixer, StreamObjectType, PIDController, HX,
        HeatExchangerSimple, SteamGenerator, Flange, NMPC, CoolingTower, CoolingTowerSimple, CoolingTowerHeatExchangerSimple, DistillationColumn, Signal,
        Selector, ControlMVSignalSplitter };

enum liquidpipeflowreference { PipeEntrance, PipeEnd };

enum materialphase { Solid, Liquid, Gas };

enum piddirection { Direct = -1, Reverse = 1 };

enum baseclasstypeenum { UnitOp, StreamBaseClassType, Block };

enum calculationmethod { DetermineFlow, DeterminePressure };

enum nmpcalgorithm { UnconstrainedLineSearch, InteriorPoint1, ActiveSet1, GeneticAlgorithm1, ParticleSwarmOptimisation1 };

enum datasourceforvar { Simulation, Exceldata };

enum typesofselector { LowSelector, HighSelector };

enum components { Naphtha = 0, Air = 1 };

const DrawModeEntity DefaultDrawModeEntity = EditMode;



class global
{
public:
    //Timing constants ---------------------------------------------------------------------------------------------------------------------------
    static int TimerInterval; // micro seconds
    static double SpeedUpFactor; //50, CT and normal sim: 200; //factor  Heat exchangers: 30000
    static double SampleT;
    static double SimVectorUpdateT; //CT: 1.0; //seconds; HX: 30s; normal sim with cooling tower: 1s.
    static int SimVectorUpdatePeriod; //Nr. samples between saving in vect.
    static double TrendUpdateT; //seconds.  The simulation period of updating trends in simulation.
                                                     //; HX: 30s; normal sim with cooling tower: 1s.
    static int TrendUpdateIterPeriod;
    static double SimTime; // Valve stepping: 3600.0*1; Normal full model: 3600.0*4 ; CT alone:  64830; //seconds 3600*1;//3600*24;  135: ES004 53340 ; for 172EP004: 34440;for CT fitting: 12 hours/3hours.
    static int SimIterations;
    static int SimVectorLength;
    static double PlantDataSampleT; //30 seconds for CT;  The sample frequency of the IP.21 data used for fitting.
    static double *simtimevector;

    global();
    ~global();

    double calcSampleT();
    //void initsimtimevector();


};


#endif // GLOBAL_H
