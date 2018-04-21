#include <math.h>
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

std::vector<component> &generate_fluid_pacage()
{
    //public molecule(string anabreviation, string aname, double amolarmass, double adynamicviscosity, double adensity, double aTc = 500, double aPc = 50*100000,
    //double aomega = 0.3, double aCpA = -4.224, double aCpB = 0.3063, double aCpC = -1.586e-04, double aCpD = 3.215e-08)

    std::vector<component> fluidpackage;

    fluidpackage.push_back(component(new molecule("Naphtha", "GTL Naphtha", 0.157, 0.00164, 661.4959, 273 + 495, 1.2411 * pow(10, 7),
                    -1 - log10(110000 /1.2411 * pow(10, 7)) , 150.5, 0.6, 0, 0), 0));

    fluidpackage.push_back(component(new molecule("Air", "Air", 0.02897, 1.983 * pow(10, -5), 1.225, 132.41, 3.72 * pow(10, 6),
                    0.0335,
                    0.8*31.15  + 0.2*28.11, 0.8*(-0.01357)  + 0.2*(-3.7) * pow(10, -6), 0.8*2.68*pow(10,-5)  + 0.2*1.746 * pow(10, -5),
                0.8 * (-1.168) * pow(10, -8) + 0.2 * (-1.065) * pow(10, -8)), 0));
                //Density: 1.977 kg/m3 (gas at 1 atm and 0 °C)

    fluidpackage.push_back(component(new molecule("CO2", "Carbon Dioxide", 0.018, 0.07 * 0.001, 1.977, 304.25, 7.39 * pow(10, 6), 0.228,
                    19.8, 0.07344, -5.602E-05, 1.715E-08), 0));
                //Density: 1.977 kg/m3 (gas at 1 atm and 0 °C)

    fluidpackage.push_back(component(new molecule("CO", "Carbon Monoxide", 0.02801, 0.0001662 * 0.001, 1.145), 0));
                //Density: 1.145 kg/m3 at 25 °C, 1 atm

    fluidpackage.push_back( component(new molecule("H2", "Hydrogen", 0.0020158, 8.76 * pow(10, -6), 0.08988), 0));
                //Density: 0.08988 g/L = 0.08988 kg/m3 (0 °C, 101.325 kPa)

    fluidpackage.push_back( component(new molecule("He", "Helium", global::HeMolarMass, 0, 0.1786, 5.1953, 5.1953E6, -0.390,
                    20.8, 0, 0, 0), 0));
                //essentially no viscosity.

    fluidpackage.push_back( component(new molecule("CH4", "Methane", 0.01604, 0.0001027 * 0.001, 0.6556), 0));
                //Density: 0.6556 g L−1 = 0.6556 kg/m3

    fluidpackage.push_back( component(new molecule("CH4O", "Methanol", 0.03204, 5.9E-04, 791.8, 513, 80.9 * 100000, 0.556,
                    21.15, 0.07092, 2.587E-05, -2.852E-08), 0));
                //Density: 0.6556 g L−1 = 0.6556 kg/m3

    fluidpackage.push_back( component(new molecule("N", "Nitrogen", 0.028, 0.018 * 0.001, 1.251,126.192, 3.3958*pow(10,6), 0.04,
                    31.15, -0.01357, 2.68*pow(10,-5), -1.168*pow(10,-8)), 0));
                //Density: 1.251 g/L = 1.251 kg/m3
    fluidpackage.push_back( component(new molecule("O2", "Oxygen", 0.016, 2.04 * pow(10, -5), 1.429, 154.581, 5.043 * pow(10, 6),
                    0.022, 28.11, -3.7 * pow(10, -6), 1.746 * pow(10, -5), -1.065 * pow(10, -8)), 0));

    fluidpackage.push_back( component(new molecule("H2O", "Water", 0.0180153, global::WaterDensity, 1, 647.096, 22060000, 0.344,
                    7.243e01, 1.039e-2, -1.497e-6, 0), 1.0));
                //Density: 1000 kg/m3
                //Dynamic viscosity for water at 20 Deg C

    fluidpackage.push_back( component(new molecule("C2H6", "Ethane", 0.03007, 0, 1.3562 * 100), 0)); //Just assume no viscosity for the moment.
    fluidpackage.push_back( component(new molecule("C3H8", "Propane", 0.03007, 0, 1.3562 * 100, 369.8, 42.5 * 100000, 0.153,
                    -4.224, 0.3063, -1.586e-04, 3.215e-08), 0)); //Just assume no viscosity for the moment.
    fluidpackage.push_back( component(new molecule("C4H10", "Butane", 58.12 / 1000, 0, 2.48 * 100, 425.2, 38.0 * 100000, 0.199,
                    9.487, 0.3313, -1.108e-04, -2.822e-09), 0)); //Just assume no viscosity for the moment.
    fluidpackage.push_back( component(new molecule("C5H12", "Pentane", 72.15 / 1000, 240 / 1000000,
                    0.626 * 1000), 0));
    fluidpackage.push_back( component(new molecule("C6H14", "2-Methylpentane", 86.18 / 1000, 0,
                    653), 0));
    fluidpackage.push_back( component(new molecule("C7H16", "Heptane", 100.20 / 1000, 386 / 1000000,
                    679.5), 0));
    fluidpackage.push_back( component(new molecule("C8H18", "Octane", 114.23 / 1000, 542 / 1000000,
                    0.703 * 1000), 0));
    fluidpackage.push_back( component(new molecule("C9H20", "Nonane", 128.26 / 1000, 0.711 / 1000, 718), 0));
    fluidpackage.push_back( component(new molecule("C10H22", "Decane", 142.28 / 1000, 0.920 / 1000, 730, 617.8, 21.1 * 100000), 0));
    fluidpackage.push_back( component(new molecule("C11H24", "Undecane", 156.30826 / 1000, 0.920 / 1000, 740.2), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C12H26", "Dodecane", 170.33 / 1000, 1.35 / 1000, 780.8), 0));
    fluidpackage.push_back( component(new molecule("C13H28", "Tridecane", 184.36 / 1000, 0, 756), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C14H30", "Tetradecane", 198.39 / 1000, 2.18 / 1000, 756), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C15H32", "Pentadecane", 212.41 / 1000, 0, 769), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C16H34", "Hexadecane", 226.44 / 1000, 3.34 / 1000, 770), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C17H36", "Heptadecane", 240.47 / 1000, 0, 777), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C18H38", "Octadecane", 254.494 / 1000, 0, 777), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C19H40", "Nonadecane", 268.5209 / 1000, 0, 786), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C20H42", "Icosane", 282.55 / 1000, 0, 786), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C21H44", "Heneicosane", 296.6 / 1000, 0, 792), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C22H46", "Docosane", 310.61 / 1000, 0, 778), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C23H48", "Tricosane ", 324.63 / 1000, 0, 797), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C24H50", "Tetracosane", 338.66 / 1000, 0, 797), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C25H52", "Pentacosane", 352.69 / 1000, 0, 801), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C26H54", "Hexacosane", 366.71 / 1000, 0, 778), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C27H56", "Heptacosane", 380.74 / 1000, 0, 780), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C28H58", "Octacosane", 394.77 / 1000, 0, 807), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C29H60", "Nonacosane", 408.80 / 1000, 0, 808), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C30H62", "Triacontane", 422.82 / 1000, 0, 810), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C31H64", "Hentriacontane", 436.85 / 1000, 0, 781), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C32H66", "Dotriacontane", 450.88 / 1000, 0, 812), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C33H68", "Tritriacontane", 464.90 / 1000, 0, 811), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C34H70", "Tetratriacontane ", 478.93 / 1000, 0, 812), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C35H72", "Pentatriacontane ", 492.96 / 1000, 0, 813), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C36H74", "Hexatriacontane", 506.98 / 1000, 0, 814), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C37H76", "Heptatriacontane", 520.99 / 1000, 0, 815), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C38H78", "Octatriacontane", 535.03 / 1000, 0, 816), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C39H80", "Nonatriacontane", 549.05 / 1000, 0, 817), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C40H82", "Tetracontane", 563.08 / 1000, 0, 817), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C41H84", "Hentetracontane", 577.11 / 1000, 0, 818), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C42H86", "Dotetracontane", 591.13 / 1000, 0, 819), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C43H88", "Triatetracontane", 605.15 / 1000, 0, 820), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C44H90", "Tetratetracontane", 619.18 / 1000, 0, 820), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C45H92", "Pentatetracontane", 633.21 / 1000, 0, 821), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C46H94", "Hexatetracontane", 647.23 / 1000, 0, 822), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C47H96", "Heptatetracontane", 661.26 / 1000, 0, 822), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C48H98", "Octatetracontane", 675.29 / 1000, 0, 823), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C49H100", "Nonatetracontane", 689.32 / 1000, 0, 823), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C50H102", "Pentacontane", 703.34 / 1000, 0, 824), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C51H104", "Henpentacontane", 717.37 / 1000, 0, 824), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C52H106", "Dopentacontane", 731.39 / 1000, 0, 825), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C53H108", "Tripentacontane", 745.42 / 1000, 0, 825), 0)); // Just assume zero viscosity for now.
    fluidpackage.push_back( component(new molecule("C54H110", "Tetrapentacontane", 759.45 / 1000, 0, 826), 0)); // Just assume zero viscosity for now.

                //public molecule(string anabreviation, string aname, double amolarmass, double adynamicviscosity (Pa·s), double adensity, double adefaultmolefraction)
    return fluidpackage;
}

std::vector<std::string> global::objecttypes_strings =
    { "FTReactor", "GasPipe", "LiquidPipe", "Pump", "Tank", "Valve", "Tee", "Mixer", "StreamObjectType", "PIDController", "HX",
       "HeatExchangerSimple", "SteamGenerator", "Flange", "NMPC", "CoolingTower", "CoolingTowerSimple", "CoolingTowerHeatExchangerSimple",
       "DistillationColumn", "Signal", "Selector", "ControlMVSignalSplitter" };

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

//Scientific constants -----------------------SOME OF THESE WILL LATER BE ABLE TO BE MOVED TO THE FLUID PACKAGE-----------------------------------------------------------------------------
double global::g = 9.81; //m/s^2
double global::R = 8.314; // J/Kelvin/molecule
double global::Ps = 100000; //Pa
double global::Ts = 273.15; //Calvin
        //Water constants
double global::WaterDensity = 1000; //kg/m^3
double global::DeltaHWater = 2257*1000; //Joule / kg
        //Air constants
double global::AirDensity = 1.225; //kg/m^3

//Differential constants ----------------------------------------------------------------------------------------------------------------------
double global::limithm = 0.00001; //0.00001; //h in the limit to zero as a multiplier for derivative calculation.
double global::limitjacnmpc = 0.00001;  //0.00001  h in the limit to zero as a multiplier for derivative calculation for nmpc.
double global::limitjacnmpcadd = 0.0001;  //h in the limit to zero as an added infinitesimal constant for derivative calculation for nmpc.
int global::RungaKuta = 4; //Number of runga kuta iteration calculations for integration of the diff equations.

        //Calculation constants
double global::Epsilon = 0.00000001; //small number that is added to a denominator in order to divide by it safely.
double global::ConvergeDiffFrac = 0.001; //Fraction difference or less that will be treated as convergance.

//Screen constants
double global::GScale = 1600.0 / 100.0; //pixels per m
int global::OriginX = 0; //pixels; 10 in csharp
int global::OriginY = 0; //pixels; 60 in csharp
double global::DefaultLocationX = 0;
double global::DefaultLocationY = 0;
double global::MinDistanceFromPoint = 0.5; //m : Minimum Distance from each point for it to be selected

double global::MinDistanceFromGasPipe = 0.5; //m : Minimum distance from each gas pipe for it to be selected

        //Simulation-wide constants
double global::RelativeHumidity = 80; //%.  Average RH for Doha through the year.  Might need to make this an instantaneous figure later
                                                    //    on.  wwww.qatar.climatemps.com/humidity.php
double global::AmbientTemperature = 25 + 273; //K.  To be converted to degC for Stull equation.  This is the average daily mean
                                                            //through out the year.  To be changed later based on when the simulation is run.
                                                            //Currently from en.wikipedia.org/wiki/Doha#Climate
double global::HeCircuitFlowT0 = 85.0; //kg/s
double global::SGPGasInlet = 65 * 100000; //Pa
double global::SGPGasOutlet = 58.5 * 100000; //Pa
double global::PBPGasInlet = 6.9e6; //Pa
double global::PBPGasOutlet = 6.5e6; //Pa
double global::PBTGasInlet = 259.1 + 273; //Kelvin.
double global::PBTGasOutlet = 700 + 273; //Kelvin.

//Material properties
double global::AirMolarMass = 0.02897; //kg/mol
double global::CO2MolarMass = 0.018; //kg/mol
double global::H2OMolarMass = 0.0180153; //kg/mol
double global::HeMolarMass = 4.002602e-3; //kg/mol

        //in and out point constants
double global::InOutPointWidth = 0.2; //m
double global::InOutPointHeight = 0.08; //m

        //baseprocessclass class default properties
double global::baseprocessclassInitMass = 10; //kg
double global::baseprocessclassInitMassFlow = 0; //kg/h
double global::baseprocessclassInitPressure = 1 * Ps; //Pa    210*Ps
double global::baseprocessclassInitTemperature = Ts + 25; //Kelvin.  25 degC   Ts + 170
double global::baseprocessclassInitVolume = 1; //m3

//chromosome class constants ------------------------------------------------------------------------------------------------------------------
int global::DefaultMaxValueforMVs = 100; //the top end of the frac/perc range of each MV.
int global::MaxBinaryLengthChromosome = 7;

        //complex class constants ----------------------------------------------------------------------------------------------------------------------
double global::ZeroImaginary = 0.0001;

//coolingtower class constants (CT) ---------This class is initially modelled on the Swedish paper, Marques later added-------------------------------------------------------------------------

int global::CTRK4ArraySize = 5;

double global::CTHeight = 14.1;  //m
double global::CTWidth = 14.4; //m
double global::CTLength = 14.4; //m
int global::CTDefaultNrStages = 10;  //The number of discretisations of the water, interface and air streams

double global::CTDefaultFanDP = 211.9; //Pa.  Calculated on model Excel sheet.
double global::CTDefaultFanTotalEfficiency = 0.866; //Fraction.  From CT datasheet.
double global::CTFanSpeedT0 = 120.1 / 60.0; //revolutions per second; from data sheet for fan.
double global::CTFanPowerT0 = 137000; //W; as per CT datasheet.
double global::CTFanShutdownSpeed = 0.1; //rps.  Speed at shutdown to keep simulating the air that does exchange heat with the water when the fan is turned off.

        //Below are the model paramaters for the second order power transient for the fans and pumps
double global::RotatingPercOS = 50; //% percentage overshoot
double global::RotatingTsettle = 15; //seconds; settling time.
double global::RotatingZeta = -log(global::RotatingPercOS/100)/(sqrt(pow(M_PI,2) + pow(log(global::RotatingPercOS/100),2)));
        //public static double RotatingOmegaN = 4 / (RotatingZeta * RotatingTsettle);
double global::RotatingOmegaN = -log(0.02 * sqrt(1 - pow(global::RotatingZeta, 2))) / (global::RotatingZeta * global::RotatingTsettle);
        //public static double RotatingZeta = 0.4;
        //public static double RotatingOmegaN = 1 / 8;

double global::Rotatingb0 = pow(global::RotatingOmegaN, 2);
double global::Rotatinga0 = pow(global::RotatingOmegaN, 2);
double global::Rotatinga1 = 2 * global::RotatingZeta * global::RotatingOmegaN;
double global::Rotatinga2 = 1;

double global::CTTotalInterfaceArea = global::CTWidth * global::CTHeight;
double global::CTTotalHorizontalArea = global::CTWidth * global::CTLength;
double global::CTTotalVolume = global::CTWidth * global::CTLength * global::CTHeight;
double global::CTFillVolume = 1161; //m^3;  as per CT data sheet.
double global::CTDefaultSegmentVolume = global::CTWidth * global::CTLength * global::CTHeight / global::CTDefaultNrStages;

double global::CTWaterVolumeFraction = 0.1;
double global::CTPackingVolumeFraction = global::CTFillVolume / global::CTTotalVolume;
double global::CTAirVolumeFraction = 1 - global::CTWaterVolumeFraction - global::CTPackingVolumeFraction;

double global::CTDropletRadius = 0.001; //m; Assumption  0.001
double global::CTDropletVolume = 4 / 3 * M_PI * pow(global::CTDropletRadius, 3);

double global::CTDropletSurfaceArea = 4 * M_PI * pow(global::CTDropletRadius, 2);

double global::CTLewisFactor = 1; //This is from Lewis' work.
double global::CTCpAir = 1013; //J/kgK; At 400K.

double global::CTTransferCoefCoef = 1 / global::CTDefaultSegmentVolume; //Multiplier to scale all the
double global::CTDefaultMassTransferCoefficientAir = 0.000657; // kg/(s*m^2) 0.0001; CTTransferCoefCoef*2.71E-14; //; fitted value . 9.2688E-07
double global::CTDefaultHeatTransferCoefficientWater = 64.395; //W/(m^2*K) . 1.0; fitted value .CTTransferCoefCoef*14.814
double global::CTDefaultHeatTransferCoefficientAir = 0.6658;  //W/(m^2*K) . 1.0; fitted value .CTTransferCoefCoef*5.729;

int global::CTNIn = 2;  //One flow in is water (strm1), the other is air (strm2).
int global::CTNOut = 2;

double global::AAntoineWater = 8.07131; //Antoine equation coefficients for water vapour pressure.  This is for the equation yielding mmHg
double global::BAntoineWater = 1730.63;
double global::CAntoineWater = 233.426;
double global::AbsHumidityConst = 2.16679 / 1000.0; //kg*K/J
double global::ConvertmmHgtoPa = 133.3223; //Pa per mmHg

double global::WaterSatPressC1 = -7.85951;
double global::WaterSatPressC2 = 1.844;
double global::WaterSatPressC3 = -11.786;
double global::WaterSatPressC4 = 22.68;
double global::WaterSatPressC5 = -15.96;
double global::WaterSatPressC6 = 1.801;

double global::BuckC1 = 611.21;
double global::BuckC2 = 18.678;
double global::BuckC3 = 234.5;
double global::BuckC4 = 257.14;

double global::CTTuningFactor = 0.9; //factor to throttle the amount of cooling for tuning the total model purposes (to be removed later).

double global::CTHeightDraw = 5.0; //meter
double global::CTWidthDraw = 5.0; //meter
double global::CTInPointsFraction[] = { 0.1, 0.9 }; //input 1: Cooling water return
double global::CTOutPointsFraction[] = { 0.10, 0.9 }; //Output 1: Cooling water supply

double global::CTTemperatureTau = 15 * 60; //seconds.  From Muller Craig article.

        //Strm1 is normally the warm stream, and Strm2 the cold stream.  So for the CT strm1 will then be the water that is coming in, and strm2 the air.
double global::CTDefaultU = 497; //from datasheet    330 * 1000000 / 3600; //W/(m^2*K);  Taken from the Muller/Craig article and converted to SI units.
double global::CTDefaultA = 287; //m^2 ; From the Muller/Craig article this figure would have been 100.

double global::CTMassFlowStrm0T0 = 6169960 / 3600.0; //kg/s , based on Flows to Equipment sheet.
double global::CTMolFlowStrm0T0 = global::CTMassFlowStrm0T0 / global::H2OMolarMass; //MOLAR MASS here is for the water flow then back from the plant.

double global::CTMassFlowStrm1T0 = 2340671/3600.0; //kg/s , from fitted data per cooling tower.
double global::CTMolFlowStrm1T0 = global::CTMassFlowStrm1T0 / global::AirMolarMass; //MOLAR MASS TO BE CHANGED HERE TO BE GENERIC

double global::CTPStrm0Inlet = 2.0 * Ps + Ps; //Pa  2.0 barg
double global::CTPStrm0Outlet = global::Ps;
double global::CTPStrm1Inlet = global::Ps;
double global::CTPStrm1Outlet = global::CTPStrm1Inlet - 0.1 * global::Ps; //THE FAN MODEL WILL NEED TO BE ADDED HERE LATER TO MAKE THIS MORE ACCURATE.

double global::CTTStrm0Inlet = 273 + 45; //Kelvin (water)
double global::CTTStrm0Outlet = 273 + 35; //Kelvin (water)
double global::CTTStrm1Inlet = global::AmbientTemperature; //Kelvin (air)
double global::CTTStrm1Outlet = 273 + 39.6; //Kelvin (air)

double global::CTTInterfaceT0 = 0.5 * (global::CTTStrm0Inlet + global::CTTStrm1Inlet);

double global::CTStrm0ValveOpeningDefault = 1.0; //fraction
double global::CTStrm0Cv = global::CTMassFlowStrm0T0 / global::CTStrm0ValveOpeningDefault /
            sqrt((global::CTPStrm0Inlet - global::CTPStrm0Outlet) / global::WaterDensity);

double global::CTPMaxFactorIncreaseperSampleT = 2.0;

double global::CTStrm1FlowCoefficient = global::CTMassFlowStrm1T0 / sqrt((global::CTPStrm1Inlet - global::CTPStrm1Outlet) / global::AirDensity);
double global::CTStrm0TempTau = global::CTTemperatureTau; //seconds.
double global::CTStrm1TempTau = global::CTTemperatureTau; //seconds.
double global::CTStrm0FlowTau = 60; //seconds.  Based on Muller-Craig.
double global::CTStrm1FlowTau = 60; //seconds.  Based on Muller-Craig.

//coolingtowerheatexchangersimple class constants (CTHES) -------------------------------------------------------------------------------------

int global::CTHESNIn = 2;  //One flow in is water (strm1), the other is air (strm2).
int global::CTHESNOut = 2;

double global::CTHESTuningFactor = 0.9; //factor to throttle the amount of cooling for tuning the total model purposes (to be removed later).

double global::CTHESHeight = 5.0; //meter
double global::CTHESWidth = 5.0; //meter
double global::CTHESInPointsFraction[] = { 0.1, 0.9 }; //input 1: Cooling water return
double global::CTHESOutPointsFraction[] = { 0.10, 0.9}; //Output 1: Cooling water supply

double global::CTHESTemperatureTau = 15 * 60; //seconds.  From Muller Craig article.

        //Strm1 is normally the warm stream, and Strm2 the cold stream.  So for the CT strm1 will then be the water that is coming in, and strm2 the air.
double global::CTHESDefaultU = 497; //from datasheet    330 * 1000000 / 3600; //W/(m^2*K);  Taken from the Muller/Craig article and converted to SI units.
double global::CTHESDefaultA = 287; //m^2 ; From the Muller/Craig article this figure would have been 100.

double global::CTHESMassFlowStrm1T0 = 6169960 / 3600.0; //kg/s , based on Flows to Equipment sheet.
double global::CTHESMolFlowStrm1T0 = global::CTHESMassFlowStrm1T0 / global::H2OMolarMass; //MOLAR MASS here is for the water flow then back from the plant.
double global::CTHESMassFlowStrm2T0 = global::CTHESMassFlowStrm1T0; //kg/s , FOR NOW. TO BE CHANGED LATER.
double global::CTHESMolFlowStrm2T0 = global::CTHESMassFlowStrm2T0 / global::AirMolarMass; //MOLAR MASS TO BE CHANGED HERE TO BE GENERIC

double global::CTHESPStrm1Inlet = 3.7 * global::Ps + global::Ps; //Pa
double global::CTHESPStrm1Outlet = global::Ps;
double global::CTHESPStrm2Inlet = global::Ps;
double global::CTHESPStrm2Outlet = global::CTHESPStrm2Inlet - 0.1*global::Ps; //THE FAN MODEL WILL NEED TO BE ADDED HERE LATER TO MAKE THIS MORE ACCURATE.

double global::CTHESTStrm1Outlet = 373 + 35; //Kelvin (water)
double global::CTHESTStrm2Outlet = 373 + 39.6; //Kelvin (air)

double global::CTHESStrm1FlowCoefficient = global::CTHESMassFlowStrm1T0 / sqrt((global::CTHESPStrm1Inlet - global::CTHESPStrm1Outlet) / global::WaterDensity);

double global::CTHESStrm2FlowCoefficient = global::CTHESMassFlowStrm2T0 / sqrt((global::CTHESPStrm2Inlet - global::CTHESPStrm2Outlet) / global::AirDensity);
double global::CTHESStrm1TempTau = global::CTHESTemperatureTau; //seconds.
double global::CTHESStrm2TempTau = global::CTHESTemperatureTau; //seconds.
double global::CTHESStrm1FlowTau = 60; //seconds.  Based on Muller-Craig.
double global::CTHESStrm2FlowTau = 60; //seconds.  Based on Muller-Craig.

//coolingtowersimple class constants (CTS) ---------------------------------------------------------------------------------------------------
double global::CTSApproach = 4.5;  //K or degC (different in Temps so either unit applies).  From Muller&Craig article.
int global::CoolingTowerSimpleNIn = 1;  //For now this will just be the main process flow in and out for simplicity,  can be changed again
        //later to add more complexity.
int global::CoolingTowerSimpleNOut = 1;
double global::CTSFlowMakeUp = 144.0 * 1000 / 3600; //kg/s  Taken from ORYX GTL plant data for the cooling tower make-up.
double global::CTSVaporisationFraction = 0.00153; //fraction per degree C
double global::CTSTuningFactor = 0.1; //factor to throttle the amount of cooling for tuning the total model purposes (to be removed later).

double global::CTSHeight = 5.0; //meter
double global::CTSWidth = 5.0; //meter
double global::CTSInPointsFraction[] = { 0.1 }; //input 1: Cooling water return
double global::CTSOutPointsFraction[] = { 0.50 }; //Output 1: Cooling water supply

double global::CTSTemperatureTau = 5 * 60; //seconds.  From Muller Craig article.

//distillation column default properties -------------------------------------------------------------------------------------------------------
double global::DistillationColumnRadius = 2.4; //m
double global::DistillationColumnHeight = 22.1; //m
int global::NTrays = 1;
int global::DistillationColumnNIn = 1;
int global::DistillationColumnNOut = 2;
double global::InitialDCTrayVolume = 0.61 * pow(0.6, 2) * M_PI; //Some dimentions from p 471
        //Chemical Process Equipment - Selection and Design
double global::InitialDCTrayU = 100; //Joule.  VROOM guess.
double global::InitialDCTrayn = 100000; //Moles.  VROOM guess.

double global::DistillationColumnInPointsFraction[] = { 0.50 }; //input 1: Gas feed, Input 2: Catalyst slurry.
double global::DistillationColumnOutPointsFraction[] = { 0.05, 0.95 }; //Output 1: Products, Output 2: Catalyst take-out.

        //embeddedtrend class constants
double global::EmbeddedTrendWidth = 10; //meters
double global::EmbeddedTrendHeight = 7.5; //meters

        //flange class constants -----------------------------------------------------------------------------------------------------------------------
double global::FlangeLength = 0.2; //m

//ftreactor class constants
double global::FTReactorRadius = 4.8; //m
double global::FTReactorHeight = 44.2; //m
int global::FTReactorNIn = 2;
int global::FTReactorNOut = 2;
double global::FTReactorInPointsFraction[] = { 0.60, 0.95 }; //input 1: Gas feed, Input 2: Catalyst slurry.
double global::FTReactorOutPointsFraction[] = { 0.50, 0.95 }; //Output 1: Products, Output 2: Catalyst take-out.

double global::FTReactorMaxVolume = 10; //m3 - Arbritrary number at this stage.
double global::FTReactorInitInventory = 50.0; //%
int global::OrigLoading = 13;
double global::CatDecayRate = -0.4 / (60 * 60 * 24); //% productivity change/tonne/second
double global::FreshCatAct = 9 / (60 * 60 * 24); //Tonnes product per second per tonne catalyst : Productivity
double global::RegenCatAct = FreshCatAct * 0.90;  //% : productivity of regenerated catalyst - It looses 10% productivity
double global::NrRegen = 4.0;
double global::LostInRegen = 5; //% : percentage weight lost during the regen process that can never be used again since it is fines - basically spent
        // catalyst.
double global::CatTake = 20 * 1000 * 12 / 365 / 24 / 3600; //kg/second out of the reactor (spent catalyst + catalyst to be regenned).  Does not include losses in regen.
double global::RegenIn = global::CatTake * global::NrRegen / (1 + global::NrRegen) * (100 - global::LostInRegen) / 100;  //tonnes per month
double global::RegenLosses = global::CatTake * global::NrRegen / (1 + global::NrRegen) * (global::LostInRegen) / 100;  //tonnes per month
double global::CatSpent = global::CatTake - global::RegenIn - global::RegenLosses;
double global::FreshCatIn = global::CatTake - global::RegenIn; //tonnes per month
double global::FreshCatLoadingConst[] = {106,
                                          80,
                                          40,
                                          30,
                                          20,
                                          20,
                                          20,
                                          20,
                                          20,
                                          20,
                                          20,
                                          20,
                                          20,
                                          20};
//gaspipe class constants
double global::PipeDefaultLength = 100; //m
double global::PipeDefaultDiameter = 0.5; //m
double global::PipeDefaultPressure = Ps; //Pa
double global::PipeDefaultTemperature = global::Ts + 25; //Kelvin
        //public static double PipeDefaultMoles = 100; //This will need to be changed.
double global::PipeDefaultFiLocation = 0.5;

        //liquidpipe class constants

std::string global::liquidpipeflowreferencestrings[] = {
            "Pipe entrance",
            "Pipe end"};

//heatexchanger class ---------------------------------------------------------------------------------------------------
double global::HEThermalPower = 5*1000000.0; //W  200.0 * 1000000   - > 200 MW.
int global::HeatExchangerNIn = 2;
int global::HeatExchangerNOut = 2;
double global::HeatExchangerInPointsFraction[] = { 0.05, 0.95 }; //input 1: Hot Gas feed, Input 2: Water
double global::HeatExchangerOutPointsFraction[] = { 0.05, 0.95 }; //Output 1: Cooled Gas, Output 2: Steam.
double global::HeatExchangerRadius = 1; //m
double global::HeatExchangerWidth = 6; //m

int global::NStrm2Coils = 455;
int global::HENSegments = 3;
int global::HENNodes = global::HENSegments + 1; //The nodes are the boundaries, and the segements are what is between the nodes.

        //This is based on the Modelling design Excel file and Areva design.

double global::HETStrm1Inlet = 188 + 273;//699.9 + 273; //K
double global::HETStrm1Outlet = 45 + 273;//244.7 + 273; //245 + 273; //K
double global::HETStrm2Inlet = 35 + 273; //35 + 273;//170.0 + 273; //K
double global::HETStrm2Outlet = 45 + 273;//325 + 273;//530 + 273; //K

double global::HEPStrm1Inlet = 10.3*100000; //65 * 100000; //Pa
double global::HEPStrm1Outlet = global::HEPStrm1Inlet - 0.97 * 100000;//58.5 * 100000; //Pa
double global::HEPStrm2Inlet = 1*100000;//210 * 100000; //Pa
double global::HEPStrm2Outlet = global::HEPStrm2Inlet - 0.46*100000;//HEPStrm2Inlet - 20 * 100000; //Pa

double global::HEPStrm1Delta = (global::HEPStrm1Inlet - global::HEPStrm1Outlet) / (global::HENSegments); //Not NSegments plus 1 since
        //outflow[0] is now going to
        //just be an extension of the final
        //segment
double global::HEPStrm2Delta = (global::HEPStrm2Inlet - global::HEPStrm2Outlet) / (global::HENSegments);

double global::HETStrm1T0[] = { global::HETStrm1Outlet, 0.5 * (global::HETStrm1Outlet + global::HETStrm1Inlet), global::HETStrm1Inlet };
double global::HETStrm2T0[] = { global::HETStrm2Inlet, 0.5 * (global::HETStrm2Inlet + global::HETStrm2Outlet), global::HETStrm2Outlet };
double global::HEPStrm1T0[] = { global::HEPStrm1Outlet, global::HEPStrm1Inlet - 2*global::HEPStrm1Delta, global::HEPStrm1Inlet - global::HEPStrm1Delta};
double global::HEPStrm2T0[] = { global::HEPStrm2Inlet - global::HEPStrm2Delta, global::HEPStrm2Inlet - 2*global::HEPStrm2Delta,
                                global::HEPStrm2Inlet - 3*global::HEPStrm2Delta };

double global::HEMassFlowStrm1T0 = 99600.0 / 3600.0 / global::NStrm2Coils; //kg/s Per tube.
double global::HEMassFlowArrayStrm1T0[] = { global::HEMassFlowStrm1T0, global::HEMassFlowStrm1T0, global::HEMassFlowStrm1T0,
                                            global::HEMassFlowStrm1T0 }; //kg/s  From AREVA design.
double global::MolFlowStrm1T0 = global::HEMassFlowStrm1T0 / global::CO2MolarMass; //MOLAR MASS TO BE CHANGED HERE TO BE GENERIC

double global::HEMassFlowStrm2T0 = 477000.0 / 3600.0 / global::NStrm2Coils; //kg/s  Needs to be slighly higher
double global::MolFlowStrm2T0 = global::HEMassFlowStrm2T0 / global::H2OMolarMass; //MOLAR MASS TO BE CHANGED HERE TO BE GENERIC

        //heatexchanger class : Metal differential equation constants in particular
double global::HEM = 0.1;  //2.42305008; //kg/m; The mass of steam tube per unit length.  From Excel sheet where the properties

        //heatexchanger class : pressure drop / energy drop due to friction constants
        //public static double HEAddFriction = 10.0;
double global::HEStrm1AddFriction = 0.1; //1
double global::HEStrm2AddFriction = 1; //1,   10.0;
double global::HEStrm1DeltaPK[] = { 1, (global::HEPStrm1T0[1] - global::HEPStrm1T0[0]) / pow(global::MolFlowStrm1T0, 2.0),
                                                    (global::HEPStrm1T0[2] - global::HEPStrm1T0[1]) / pow(global::MolFlowStrm1T0, 2.0),
                                                    (global::HEPStrm1Inlet - global::HEPStrm1T0[2]) / pow(global::MolFlowStrm1T0, 2.0)};

double global::HEStrm2DeltaPK[] = {   (global::HEPStrm2Inlet - global::HEPStrm2T0[0]) / pow(global::MolFlowStrm2T0, 2.0),
                                                    (global::HEPStrm2T0[0] - global::HEPStrm2T0[1]) / pow(global::MolFlowStrm2T0, 2.0),
                                                    (global::HEPStrm2T0[1] - global::HEPStrm2T0[2]) / pow(global::MolFlowStrm2T0, 2.0),
                                                    1};
        //public static double[] HEStrm2DeltaPK = HEPStrm2Delta / Math.Pow(MolFlowStrm2T0, 2.0);

        //heatexchanger class : heat exchange constants
double global::HEHeatExchangeSurfaceArea = 329; //m2
double global::HEOutsideDiameterTube = 19.05 / 1000.0; //m
double global::HETubeWallThickness = 2.11 / 1000.0; //m
double global::HEInsideDiameterTube = global::HEOutsideDiameterTube - 2 * global::HETubeWallThickness; //m
double global::HENrPassesThroughShell = 2.0;
double global::HETubeCircOutside = M_PI * global::HEOutsideDiameterTube;  //m; Tube Circumferance on the outside.
double global::HETubeCircInside = M_PI * global::HEInsideDiameterTube;  //m; Tube Circumferance on the inside.
double global::HETubeCircAve = 0.5 * (HETubeCircOutside + HETubeCircInside);
double global::HEAveLengthPerTube = 6.1; //6.1; //102.6490726; //m; Lenth per tube in AREVA design as per sheet.
double global::HEAveLengthPerSegment = global::HEAveLengthPerTube / global::HENSegments;
double global::HEAStrm2 = M_PI * pow(global::HEInsideDiameterTube / 2.0, 2.0); //m2; Cross sectional area of the pipes for the steam/water
double global::HEStrm2TubeVolume = global::HEAStrm2* global::HEAveLengthPerTube; //0.042648188; //m^3
double global::HEStrm2SegmentVolume = global::HEStrm2TubeVolume / global::HENSegments; //m^3
double global::HEShellVolume = 1.040 * 6.096 - global::HEStrm2TubeVolume * global::NStrm2Coils;
double global::HEStrm1TubeVolume = global::HEShellVolume / global::NStrm2Coils; //3.641659803; //m^3
double global::HEAStrm1 = global::HEStrm1TubeVolume / global::HEAveLengthPerTube;  //m2;  From Excel sheet from AREVA design.
double global::HEEffGasTubeCircInside = 2 * sqrt(global::HEAStrm1 / M_PI) * M_PI;//2 * Math.Sqrt(HEAStrm1 / Math.PI) * Math.PI;
double global::HEStrm1SegmentVolume = global::HEStrm1TubeVolume / global::HENSegments; //m^3
double global::HERi = 0.0001;  //0.001;  // K∙s³/kg = K m^2 / W .  Thermal conductivity's reciprocal times thickness of wall.

        //public static double HEAg = 0.035476792;  //m2;  From Excel sheet from AREVA design.
double global::HEInPointsFraction[] = { 0.05, 0.95 }; //input 1: Hot Gas feed, Input 2: Water
double global::HEOutPointsFraction[] = { 0.05, 0.95 }; //Output 1: Cooled Gas, Output 2: Steam.
        //public static double HETubeDiameter = 0.023; //m; Diameter of the tubes in the steam generator.  From Excel sheet
        //public static double HEAs = Math.PI * Math.Pow(HETubeDiameter / 2.0, 2.0); //m2; Cross sectional area of the pipes for the steam/water

        //thermal resistivity of the metal times the tube thickness.
        //This kgm and Ri, needs to be backed up with some more science.  Why is it so
        //dificult to get these values and to fix things up properly?
double global::HECsi = 0.1 ; //0.8;  //Dimensionless. These values need to be backed up with some more research.  Why is it so difficult to
        //get proper values for these?
double global::HECgi = 0.1; //0.8;  //Dimensionless. These values need to be backed up with some more research.  Why is it so difficult to
        //get proper values for these?
double global::HEHeatTransferArea = global::HETubeCircAve * global::HEAveLengthPerTube / global::HENSegments;
double global::HEKgm[] = {(global::HEHeatTransferArea * (0.5 * (global::HETStrm1T0[1] - global::HETStrm2T0[1])) /
            (global::HEThermalPower*1 / global::NStrm2Coils) - global::HERi*0.5) /
            pow(global::HEMassFlowStrm1T0,-global::HECgi),
            (global::HEHeatTransferArea * (0.5 * (global::HETStrm1T0[1] - global::HETStrm2T0[1])) /
            (global::HEThermalPower*1 / global::NStrm2Coils) - global::HERi*0.5) /
            pow(global::HEMassFlowStrm1T0,-global::HECgi),
            (global::HEHeatTransferArea * (0.5 * (global::HETStrm1T0[1] - global::HETStrm2T0[1])) /
            (global::HEThermalPower*1 / global::NStrm2Coils) - global::HERi*0.5) /
            pow(global::HEMassFlowStrm1T0,-global::HECgi)
            };
double global::HEKms[] = {(global::HEHeatTransferArea * (0.5 * (global::HETStrm1T0[1] - global::HETStrm2T0[1])) /
            (global::HEThermalPower*1 / global::NStrm2Coils) - global::HERi*0.5) /
            pow(global::HEMassFlowStrm2T0,-global::HECsi),
            (global::HEHeatTransferArea * (0.5 * (global::HETStrm1T0[1] - global::HETStrm2T0[1])) /
            (global::HEThermalPower*1 / global::NStrm2Coils) - global::HERi*0.5) /
            pow(global::HEMassFlowStrm2T0,-global::HECsi),
            (global::HEHeatTransferArea * (0.5 * (global::HETStrm1T0[1] - global::HETStrm2T0[1])) /
            (global::HEThermalPower*1 / global::NStrm2Coils) - global::HERi*0.5) /
            pow(global::HEMassFlowStrm2T0,-global::HECsi)};
        //public static double[] HEKgm = new double[] {(HEHeatTransferArea * (0.5 * (HETStrm1T0[0] - HETStrm2T0[0])) /
        //    (HEThermalPower*1 / NStrm2Coils / HENSegments) - HERi*0.5) /
        //    Math.Exp(-HEMassFlowStrm1T0*HECgi),
        //    (HEHeatTransferArea * (0.5 * (HETStrm1T0[1] - HETStrm2T0[1])) /
        //    (HEThermalPower*1 / NStrm2Coils / HENSegments) - HERi*0.5) /
        //    Math.Exp(-HEMassFlowStrm1T0*HECgi),
        //    (HEHeatTransferArea * (0.5 * (HETStrm1T0[2] - HETStrm2T0[2])) /
        //    (HEThermalPower*1 / NStrm2Coils / HENSegments) - HERi*0.5) /
        //    Math.Exp(-HEMassFlowStrm1T0*HECgi)};
        //public static double[] HEKms = new double[] {(HEHeatTransferArea * 0.5 * (HETStrm1T0[0] - HETStrm2T0[0]) /
        //    (HEThermalPower*1 / NStrm2Coils / HENSegments) - HERi*0.5) /
        //    Math.Exp(-HEMassFlowStrm2T0*HECsi),
        //    (HEHeatTransferArea * 0.5 * (HETStrm1T0[1] - HETStrm2T0[1]) /
        //    (HEThermalPower*1 / NStrm2Coils / HENSegments) - HERi*0.5) /
        //    Math.Exp(-HEMassFlowStrm2T0*HECsi),
        //    (HEHeatTransferArea * 0.5 * (HETStrm1T0[2] - HETStrm2T0[2]) /
        //    (HEThermalPower*1 / NStrm2Coils / HENSegments) - HERi*0.5) /
        //    Math.Exp(-HEMassFlowStrm2T0*HECsi)};

//heatexchangersimple (HES) class constants  -----------------------------------------------------------------------------
        //Strm1 is normally the warm stream, and Strm2 the cold stream
double global::HeatExchangerSimpleDefaultU = 500; // 497 from datasheet    330 * 1000000 / 3600; //W/(m^2*K);  Taken from the Muller/Craig article and converted to SI units.
double global::HeatExchangerSimpleDefaultA = 441; // 287 m^2 ; From the Muller/Craig article this figure would have been 100.

double global::HESMassFlowStrm1T0 = 325000.0 / 3600.0; //kg/s , biggest CW exchanger in CW circuit.
double global::HESMolFlowStrm1T0 = global::HESMassFlowStrm1T0 / global::CO2MolarMass; //MOLAR MASS TO BE CHANGED HERE TO BE GENERIC
double global::HESMassFlowStrm2T0 = 1366538.0 / 3600.0; //kg/s , for the whole strm2.
double global::HESMolFlowStrm2T0 = global::HESMassFlowStrm2T0 / global::H2OMolarMass; //MOLAR MASS TO BE CHANGED HERE TO BE GENERIC

double global::HESPStrm1Inlet = 3.7*global::Ps + global::Ps; //Pa
double global::HESPStrm1Outlet = global::HESPStrm1Inlet - 0.5*global::Ps;
double global::HESPStrm2Inlet = 3.5*global::Ps + global::Ps;
double global::HESPStrm2Outlet = global::HESPStrm2Inlet - 2 * global::Ps;

double global::HESStrm1FlowCoefficient = global::HESMassFlowStrm1T0/sqrt((global::HESPStrm1Inlet - global::HESPStrm1Outlet)/
                                                                         global::WaterDensity); //SPECIFIC GRAVITY OF THE STREAM TO BE ADDED LATER TO THIS CALC
double global::HESStrm2FlowCoefficient = 10.21076364; //This is fitted in the model - average one for all exchangers in model.
double global::HESStrm1TempTau = 6 * 60; //seconds.  Based on Muller-Craig.
double global::HESStrm2TempTau = 6 * 60; //seconds.  Based on Muller-Craig.
double global::HESStrm1FlowTau = 60; //seconds.  Based on Muller-Craig.
double global::HESStrm2FlowTau = 60; //seconds.  Based on Muller-Craig.

int global::HESNSegments = 2; //For the simple heat exchanger, the inflow and outflow streams on the two sides will be modelled as the only
                                            //2 segments.
int global::HESNStrm2Coils = 1; //heatexchangersimple class will be modelled with one big coild for strm2 only.

        //heatexchangersimple class : heat exchange constants
        //public static double HEHeatExchangeSurfaceArea = 329; //m2
        //public static double HEOutsideDiameterTube = 19.05 / 1000.0; //m
        //public static double HETubeWallThickness = 2.11 / 1000.0; //m
        //public static double HEInsideDiameterTube = HEOutsideDiameterTube - 2 * HETubeWallThickness; //m
        //public static double HENrPassesThroughShell = 2.0;
        //public static double HETubeCircOutside = Math.PI * HEOutsideDiameterTube;  //m; Tube Circumferance on the outside.
        //public static double HETubeCircInside = Math.PI * HEInsideDiameterTube;  //m; Tube Circumferance on the inside.
        //public static double HETubeCircAve = 0.5 * (HETubeCircOutside + HETubeCircInside);
        //public static double HEAveLengthPerTube = 6.1; //6.1; //102.6490726; //m; Lenth per tube in AREVA design as per sheet.
        //public static double HEAveLengthPerSegment = HEAveLengthPerTube / HENSegments;
        //public static double HEAStrm2 = Math.PI * Math.Pow(HEInsideDiameterTube / 2.0, 2.0); //m2; Cross sectional area of the pipes for the steam/water
        //public static double HEStrm2TubeVolume = HEAStrm2 * HEAveLengthPerTube; //0.042648188; //m^3
double global::HESStrm2Volume = global::HEStrm2TubeVolume * global::NStrm2Coils;
double global::HESStrm2SegmentVolume = global::HESStrm2Volume / global::HENSegments; //m^3
double global::HESShellVolume = 1.040 * 6.096 - global::HESStrm2Volume;
double global::HESStrm1TubeVolume = global::HEShellVolume / global::NStrm2Coils; //3.641659803; //m^3
        //public static double HEAStrm1 = HEStrm1TubeVolume / HEAveLengthPerTube;  //m2;  From Excel sheet from AREVA design.
        //public static double HEEffGasTubeCircInside = 2 * Math.Sqrt(HEAStrm1 / Math.PI) * Math.PI;//2 * Math.Sqrt(HEAStrm1 / Math.PI) * Math.PI;
double global::HESStrm1SegmentVolume = global::HESShellVolume / global::HESNSegments; //m^3

//material class constants ---------------------------------------------------------------------------------------------------------------------
double global::Udefault = 19500 * 100000; //Joule
double global::Vdefault = 18.01528 * 100000 / 1000 / 1000; //m^3
double global::fdefault = 0.5; //Vapour molar fraction. just a value in order to get some convergiance.
std::vector<component> global::fluidpackage = generate_fluid_pacage();
materialphase global::MaterialInitPhase = Liquid;
int global::NMaterialIterations = 10;
double global::ZNotDefined = -999999.0;
double global::epsilonadd = 0.0001; //for derivatives where the denominator is added to the variable being differentiated with respect to, and not
        //multiplied.
double global::epsilonfrac = 1.001;
double global::Inf = 9999999999;

//mixer class constants ------------------------------------------------------------------------------------------------------------------------
double global::MixerLength = 1; //m
double global::MixerDistanceBetweenBranches = 1; //m
int global::MixerDefaultNIn = 2;
double global::MixerBranchThickness = 0.1; //m
double global::MixerInitRadiusDefault = global::MixerDefaultNIn * (global::MixerDistanceBetweenBranches + global::MixerBranchThickness);

//molecule class constants
double global::InitTc = 500;              //Kelvin
double global::InitPc = 50 * 100000;      //Pa;
double global::Initomega = 0.3;           //Default Acentric factor.

//nmpc class constants ----------------------------------------------------------------------------------------------------------------------
double global::NMPCWidth = 2; //m
double global::NMPCHeight = 2; //m
int global::DefaultN = 9000; //3000; //Default Optimisation horison 80
int global::DefaultInitialDelay = 0;
int global::DefaultRunInterval = 300; //Assuming TSample is 10 sec, so then the interval would be a multiple of
                                                    //that.
double global::Defaultalphak = 1.0; //0.1; //0.1   0.001; //How much of line search delta is implemented.
nmpcalgorithm global::DefaultNMPCAlgorithm = ParticleSwarmOptimisation1; //nmpcalgorithm.ParticleSwarmOptimisation1; //nmpcalgorithm.UnconstrainedLineSearch;
double global::DefaultNMPCSigma = 0.8; //Multiplier of mubarrier for each iteration of the nmpc algorithm.

        //Interior Point 1 - ConstrainedLineSearch algorithm constants
double global::DefaultMuBarrier = 0.00000000000000001; //initial value / guess , for initialisation.
double global::Defaultsvalue = 0.6; //since the constraints at this point are the valve openings in ChemSim, a fraction half will be
                                                        //the initial value of the valve openings, and thus the values will be 0.5 away from zero.
double global::DefaultsvalueMultiplier = 0.9; //1.0;
        //public static double DefaultSigma = 0.1; //The fraction multiplied with mubarrier at the end of each iteration.
double global::DefaulttauIP = 0.95; //Tau constant for Interior Point methods.
double global::DefaultIPErrorTol = 0.0001; //If the max of the norm is below this figure, then the algorithm will stop, and we are close enough to a solution.
double global::CholeskyDelta = 0.01; //small delta used in Hessian modification.
double global::CholeskyBeta = 1000000; //Big value that will be used to try and cut down D matrix values if they become too large in Hessian modification.
double global::MVMaxMovePerSampleTT0 = 0.1; //Fraction of MV range.
double global::NMPCIPWeightPreTerm = 1;
double global::NMPCIPWeightTerminal = 10;

        //Genetic Algorithm 1 constants (mostly default values for variables in the nmpc class).
int global::DefaultNrChromosomes = 12; //The total number of total solutions that will be kept in memory each iteration.
int global::DefaultNrSurvivingChromosomes = 9; //Nr of chromosomes that will be passed to the next iteration and not replaced by new random ones.
int global::DefaultNrPairingChromosomes = 8; //The nr of chromosomes of the total population that will be pairing and producing children.
double global::DefaultProbabilityOfMutation = 0.1; //The probability that a child will be mutated in one bit.
int global::DefaultNrIterations = 10; //The number of iterations until the best GA solution will be passed to the update method.
int global::DefaultCrossOverPoint = 3; //Bit index nr (starting from zero) from right to left in the binary representation, where
                                                     //cross over and mutation will start.

        //PSO constants
int global::DefaultNrContinuousParticles = 20; //The total number of total solutions that will be kept in memory each iteration.
        //public static int DefaultNrParticles = 20; //The total number of total solutions that will be kept in memory each iteration.
int global::DefaultNrBooleanParticles = global::DefaultNrContinuousParticles; //The total number of total solutions that will be kept in memory each iteration.
int global::PSOMVBoundaryBuffer = 10; //Distance from boundary that particles are put at random when they cross the boundary.
double global::PSOMaxBooleanSpeed = 1.0; //Max probability paramater for sigmoid function for boolean PSO.

//pidcontroller class constants -------------------------------------------------------------------------------------------------------------
int global::Direct = -1;
int global::Reverse = 1;
double global::PIDControllerInitRadius = 0.4; //m
double global::PIDControllerInitK = 1;
double global::PIDControllerInitI = 100;
double global::PIDControllerInitD = 0;
double global::PIDControllerInitMinOP = 0; //Engineering units.
double global::PIDControllerInitMaxOP = 1; //Engineering units.
double global::PIDControllerInitMinPV = 0; //Engineering units.
double global::PIDControllerInitMaxPV = 1; //Engineering units.

//pump class default properties -------------------------------------------------------------------------------------------------------------
double global::PumpInitMaxDeltaPressure = 6.7*2*100000; //Pa
double global::PumpInitMinDeltaPressure = 0; //Pa
double global::PumpInitMaxActualFlow = 8700000*2/2.0 / 3600.0 / 1000; //m3/s  Dividing by 2 since we will now have 2 pumps in parallel.
                                                                                    //assume a density of 1000 as well.
double global::PumpMinActualFlow = 0.01 * global::PumpInitMaxActualFlow; //This is for when the pumps is off its curve due to too high DP.
double global::PumpCurveYAxis = global::WaterDensity * g * 70; //Pa.  Making this more than the data sheet for now to make sure the pump
                                                                     //can survive well at that level.
double global::PumpCurvef1 = 8500/3600.0; //m3/s;  Actual flow.  From pump data sheet.
double global::PumpCurvep1 = global::WaterDensity * g * 50; //Pa
double global::PumpCurvef2 = 15000 / 3600.0; //m3/s; Actual flow.
double global::PumpCurveSpeedT0 = 740 / 60.0; //rev per second.  From pump data sheet.
double global::PumpSpeedTau = 60; //seconds.

double global::PumpInitActualVolumeFlow = 0; //m3/s
double global::PumpInitOn = 1; //0 for off, 1 for on
double global::PumpInitRadius = 0.4; //m
double global::PumpInitOutletLength = 0.5; //m
double global::PumpInitOutletRadius = 0.05; //m

//stream class default properties -------------------------------------------------------------------------------------------------------------
int global::SignalNrPropDisplay = 1;

        //steamgenerator class constants
double global::SteamGeneratorRadius = 2.4; //m
double global::SteamGeneratorHeight = 12.1; //m
int global::SteamGeneratorNSegments = 3;
int global::SteamGeneratorNNodes = global::SteamGeneratorNSegments + 1; //The nodes are the boundaries, and the segements are what is between the nodes.
double global::SteamGeneratorWaterTubeVolume = 0.042648188; //m^3
double global::SteamGeneratorWaterSegmentVolume = global::SteamGeneratorWaterTubeVolume / global::SteamGeneratorNSegments; //m^3
double global::SteamGeneratorGasTubeVolume = 3.641659803; //m^3
double global::SteamGeneratorGasSegmentVolume = global::SteamGeneratorGasTubeVolume / global::SteamGeneratorNSegments; //m^3
int global::SteamGeneratorNIn = 2;
int global::SteamGeneratorNOut = 2;
int global::NSteamCoils = 220; //From AREVA design.
double global::ThermalPower = 200000000.0; //W  200.0 * 1000000   - > 200 MW.
double global::TubeCircOutside = 0.092991143;  //m; Tube Circumferance on the outside.
        //This is based on the Modelling design Excel file and Areva design.
double global::TubeCircInside = 0.072256631;  //m; Tube Circumferance on the inside.
double global::TubeCircAve = 0.5 * (global::TubeCircOutside + global::TubeCircInside);
double global::AveLengthPerTube = 102.6490726; //m; Lenth per tube in AREVA design as per sheet.
double global::AveLengthPerSegment = global::AveLengthPerTube / global::SteamGeneratorNSegments;
double global::SteamGeneratorInPointsFraction[] = { 0.05, 0.95 }; //input 1: Hot Gas feed, Input 2: Water
double global::SteamGeneratorOutPointsFraction[] = { 0.05, 0.95 }; //Output 1: Cooled Gas, Output 2: Steam.

double global::Ws0T0 = 77.0 / global::NSteamCoils; //kg/s  From AREVA design.

double global::SGTGasInlet = 699.9 + 273; //K
double global::SGTGasOutlet = 244.7 + 273; //245 + 273; //K
double global::SGTWaterInlet = 170.0 + 273; //K
double global::SGTWaterOutlet = 325 + 273;//530 + 273; //K
double global::WgasT0 = global::HeCircuitFlowT0 / global::NSteamCoils; //kg/s
double global::MolFlowGasT0 = global::WgasT0 / global::HeMolarMass;
double global::WwaterT0 = 77.0 / global::NSteamCoils; //kg/s
double global::MolFlowWaterT0 = global::WwaterT0 / global::H2OMolarMass;

double global::SGPGasDelta = (global::SGPGasInlet - global::SGPGasOutlet) / (global::SteamGeneratorNSegments); //Not NSegments plus 1 since
        //outflow[0] is now going to
        //just be an extension of the final
        //segment
double global::SGPWaterInlet = 210 * 100000; //Pa
double global::SGPWaterOutlet = global::SGPWaterInlet - 20 * 100000; //Pa
double global::SGPWaterDelta = (global::SGPWaterInlet - global::SGPWaterOutlet) / (global::SteamGeneratorNSegments);

double global::TGasT0[] = { global::SGTGasOutlet, 0.5 * (global::SGTGasOutlet + global::SGTGasInlet), global::SGTGasInlet };
double global::TWaterT0[] = { global::SGTWaterInlet, 0.5 * (global::SGTWaterInlet + global::SGTWaterOutlet), global::SGTWaterOutlet };
double global::PGasT0[] = { global::SGPGasOutlet, global::SGPGasInlet - 2*global::SGPGasDelta,
            global::SGPGasInlet - global::SGPGasDelta};
double global::PWaterT0[] = { global::SGPWaterInlet - global::SGPWaterDelta,
            global::SGPWaterInlet - 2*global::SGPWaterDelta, global::SGPWaterInlet - 3*global::SGPWaterDelta };

        //steamgenerator class : Heat exchange constants
double global::Ri = 0.001;   // K∙s³/kg = K m^2 / W .  Thermal conductivity's reciprocal times thickness of wall.
        //thermal resistivity of the metal times the tube thickness.
        //This kgm and Ri, needs to be backed up with some more science.  Why is it so
        //dificult to get these values and to fix things up properly?
double global::Csi = 0.8;  //Dimensionless. These values need to be backed up with some more research.  Why is it so difficult to
        //get proper values for these?
double global::Cgi = 0.8;  //Dimensionless. These values need to be backed up with some more research.  Why is it so difficult to
        //get proper values for these?
double global::HeatTransferArea = global::TubeCircAve * global::AveLengthPerTube / global::SteamGeneratorNSegments;
        //public static double Kgm =
        //    (HeatTransferArea * (0.5 * (0.5 * (SGTGasInlet + SGTGasOutlet) - 0.5 * (SGTWaterOutlet + SGTWaterInlet))) /
        //    (ThermalPower / NSteamCoils / SteamGeneratorNSegments) - Ri * 0.5) /
        //    Math.Exp(-WgasT0 * Cgi); //THIS IS TO BE PUT BACK AS THE WAY OF INITIALISING LATER ON.
        //public static double Kgm = (HeatTransferArea * (0.5 * (0.5 * (TGasInlet + TGasOutlet) - 0.5 * (TWaterOutlet + TWaterInlet))) /
        //    (ThermalPower / NSteamCoils / SteamGeneratorNSegments)) /
        //    Math.Pow(WgasT0, -Cgi);

double global::Kgm[] = {(global::HeatTransferArea * (0.5 * (global::TGasT0[0] - global::TWaterT0[0])) /
            (global::ThermalPower*1 / global::NSteamCoils / global::SteamGeneratorNSegments) - Ri*0.5) /
            exp(-global::WgasT0*Cgi),
            (global::HeatTransferArea * (0.5 * (global::TGasT0[1] - global::TWaterT0[1])) /
            (global::ThermalPower*1 / global::NSteamCoils / global::SteamGeneratorNSegments) - global::Ri*0.5) /
            exp(-global::WgasT0*Cgi),
            (global::HeatTransferArea * (0.5 * (global::TGasT0[2] - global::TWaterT0[2])) /
            (global::ThermalPower*1 / global::NSteamCoils / global::SteamGeneratorNSegments) - Ri*0.5) /
            exp(-global::WgasT0*Cgi)};

        //public static double Kms =
        //    (HeatTransferArea * (0.5 * (0.5 * (SGTGasInlet + SGTGasOutlet) - 0.5 * (SGTWaterOutlet + SGTWaterInlet))) /
        //    (ThermalPower / NSteamCoils / SteamGeneratorNSegments) - Ri * 0.5) /
        //    Math.Exp(-WwaterT0 * Cgi);  //THIS IS TO BE PUT BACK AS THE WAY OF INITIALISING LATER ON.
double global::Kms[] = {(global::HeatTransferArea * 0.5 * (global::TGasT0[0] - global::TWaterT0[0]) /
            (global::ThermalPower*1 / global::NSteamCoils / global::SteamGeneratorNSegments) - global::Ri*0.5) /
            exp(-global::WwaterT0*Cgi),
            (global::HeatTransferArea * 0.5 * (global::TGasT0[1] - global::TWaterT0[1]) /
            (global::ThermalPower*1 / global::NSteamCoils / global::SteamGeneratorNSegments) - global::Ri*0.5) /
            exp(-global::WwaterT0*global::Cgi),
            (global::HeatTransferArea * 0.5 * (global::TGasT0[2] - global::TWaterT0[2]) /
            (global::ThermalPower*1 / global::NSteamCoils / global::SteamGeneratorNSegments) - Ri*0.5) /
            exp(-global::WwaterT0*global::Cgi)};

double global::WgT0[] = { 85.0 / global::NSteamCoils, 85.0 / global::NSteamCoils, 85.0 / global::NSteamCoils, 85.0 / global::NSteamCoils }; //kg/s  From AREVA design.
double global::Ag = 0.035476792;  //m2;  From Excel sheet from AREVA design.
double global::EffGasTubeCircInside = 2 * sqrt(global::Ag / M_PI) * M_PI;
double global::TubeDiameter = 0.023; //m; Diameter of the tubes in the steam generator.  From Excel sheet
double global::As = M_PI * pow(global::TubeDiameter / 2.0, 2.0); //m2; Cross sectional area of the pipes for the steam/water

        //steamgenerator class : Metal differential equation constants in particular
double global::M = 2.42305008; //kg/m; The mass of steam tube per unit length.  From Excel sheet where the properties

        //steamgenerator class : pressure drop / energy drop due to friction constants
double global::SteamGenAddFriction = 10.0;
double global::SteamGenGasAddFriction = 1;
double global::SGWaterDeltaPK = global::SGPWaterDelta / pow(global::MolFlowWaterT0, 2.0);
double global::SGGasDeltaPK = global::SGPGasDelta / pow(global::MolFlowGasT0, 2.0);

        //steamgenerator class : MOMENTUM EQUATION constants
double global::DynViscA = 0.2259; //Coefficient in front of Dynamic Viscocity equation as per Excel sheet.
double global::DynViscB = -0.018; //Exponent coefficient of Dynamic Viscocity equation as per Excel sheet.

        //public static double[] Kfric = new double[] {(PsaveT0[0] - PsaveT0[1])/Math.Pow(Ws0T0,2),
        //                                             (PsaveT0[1] - PsaveT0[2])/Math.Pow(Ws0T0,2),
        //                                             (PsaveT0[2] - PsaveT0[3])/Math.Pow(Ws0T0,2)};

//stream class default properties -------------------------------------------------------------------------------------------------------------
double global::MinDistanceFromStream = 15; //pixels  0.5; //m Minimum distance from each stream for it to be selected
double global::StreamArrowAngle = 30.0 / 180.0 * M_PI; //radians
double global::StreamArrowLength = 0.5; //m
double global::StreamMaxMassFlow = 100000000; //kg/s (100,000 tps)
double global::StreamMinMassFlow = -global::StreamMaxMassFlow;
int global::StreamNrPropDisplay = 3;

//tank class default properties ----------------------------------------------------------------------------------------------------------------
double global::TankInitRadius = 22; //12.88280927; //m; from Cooling tower sump design info.
double global::TankInitHeight = 13.2; //13.2; //m; //from cooling tower sump design.
double global::TankRadiusDraw = 6.0; //meter
double global::TankHeightDraw = 2.0; //meter
double global::TankInitMaxVolume = M_PI * pow(global::TankInitRadius, 2) * global::TankInitHeight; //m3
double global::TankInitFracInventory = 0.5; //Fraction
double global::TankInitInOutletDistanceFraction = 0.95; //fraction from the bottom or top that the inlet or outlet of the tank will be situated.
double global::TankMinFracInventory = 0.02; //fraction

//tee class constants --------------------------------------------------------------------------------------------------------------------------
double global::TeeLength = 1; //m
double global::TeeDistanceBetweenBranches = 1; //m
int global::TeeDefaultNOut = 2;
double global::TeeBranchThickness = 0.1; //m
double global::TeeInitRadiusDefault = global::TeeDefaultNOut * (global::TeeDistanceBetweenBranches + global::TeeBranchThickness);

        //trend class constants -----------------------------------------------------------------------------------------------------------------------
float global::XAxisMargin = 20; //pixels
float global::YAxisMargin = 20; //pixels
double global::YIncrement = 0.1; //fraction
double global::YmaxMaxFactorAbove = 2.0; //factor

//valve class constants -----------------------------------------------------------------------------------------------------------------------
double global::ValveDefaultActualFlow = 800 / 3600.0; //From average HX flow in unit CW circuit.
double global::ValveDefaultDP = 1.2 * Ps;
double global::ValveDefaultOpening = 0.5; //Default valve opening.

double global::ValveEqualPercR = 40; //Dimensionless constant for valve equalpercentage.
double global::ValveDefaultCv = global::ValveDefaultActualFlow /
        (pow(global::ValveEqualPercR, global::ValveDefaultOpening - 1) * sqrt(global::ValveDefaultDP)); //m^3/s/(Pa^0.5)
double global::ValveHydraulicTau = 10.0; //seconds.  Time constant of valve hydraulics.

double global::ValveLength = 0.5; //m
double global::ValveWidth = 0.4; //m


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







