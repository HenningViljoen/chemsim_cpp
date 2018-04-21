#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <vector>
#include "component.h"

enum DrawModeEntity {EditMode, ValveMode};

enum objecttypes { FTReactor, GasPipe, LiquidPipe, Pump, Tank, Valve, Tee, Mixer, StreamObjectType, PIDController, HX,
        HeatExchangerSimple, SteamGenerator, Flange, NMPC, CoolingTower, CoolingTowerSimple, CoolingTowerHeatExchangerSimple,
                   DistillationColumn, Signal, Selector, ControlMVSignalSplitter };



enum liquidpipeflowreference { PipeEntrance, PipeEnd };

enum materialphase { Solid, Liquid, Gas };

enum piddirection { PIDDirectionDirect = -1, PIDDirectionReverse = 1 };

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
    static std::vector<std::string> objecttypes_strings;

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

    //Scientific constants -----------------------SOME OF THESE WILL LATER BE ABLE TO BE MOVED TO THE FLUID PACKAGE-----------------------------------------------------------------------------
    static double g;
    static double R;
    static double Ps;
    static double Ts;
            //Water constants
    static double WaterDensity;
    static double DeltaHWater;
            //Air constants
    static double AirDensity;

    //Differential constants ----------------------------------------------------------------------------------------------------------------------
    static double limithm;
    static double limitjacnmpc;
    static double limitjacnmpcadd;
    static int RungaKuta;

            //Calculation constants
    static double Epsilon;
    static double ConvergeDiffFrac;

    //Screen constants
    static double GScale;
    static int OriginX;
    static int OriginY;
    static double DefaultLocationX;
    static double DefaultLocationY;
    static double MinDistanceFromPoint;

    static double MinDistanceFromGasPipe;

            //Simulation-wide constants
    static double RelativeHumidity;
    static double AmbientTemperature;
    static double HeCircuitFlowT0;
    static double SGPGasInlet;
    static double SGPGasOutlet;
    static double PBPGasInlet;
    static double PBPGasOutlet;
    static double PBTGasInlet;
    static double PBTGasOutlet;

    //Material properties
    static double AirMolarMass;
    static double CO2MolarMass;
    static double H2OMolarMass;
    static double HeMolarMass;


            //in and out point constants
    static double InOutPointWidth;
    static double InOutPointHeight;

            //baseprocessclass class default properties
    static double baseprocessclassInitMass;
    static double baseprocessclassInitMassFlow;
    static double baseprocessclassInitPressure;
    static double baseprocessclassInitTemperature;
    static double baseprocessclassInitVolume;

    //chromosome class constants ------------------------------------------------------------------------------------------------------------------
    static int DefaultMaxValueforMVs;
    static int MaxBinaryLengthChromosome;

            //complex class constants ----------------------------------------------------------------------------------------------------------------------
    static double ZeroImaginary;

    //coolingtower class constants (CT) ---------This class is initially modelled on the Swedish paper, Marques later added-------------------------------------------------------------------------
    static int CTRK4ArraySize;

    static double CTHeight;
    static double CTWidth;
    static double CTLength;
    static int CTDefaultNrStages;

    static double CTDefaultFanDP;
    static double CTDefaultFanTotalEfficiency;
    static double CTFanSpeedT0;
    static double CTFanPowerT0;
    static double CTFanShutdownSpeed;

            //Below are the model paramaters for the second order power transient for the fans and pumps
    static double RotatingPercOS;
    static double RotatingTsettle;
    static double RotatingZeta;
    static double RotatingOmegaN;

    static double Rotatingb0;
    static double Rotatinga0;
    static double Rotatinga1;
    static double Rotatinga2;

    static double CTTotalInterfaceArea;
    static double CTTotalHorizontalArea;
    static double CTTotalVolume;
    static double CTFillVolume;
    static double CTDefaultSegmentVolume;

    static double CTWaterVolumeFraction;
    static double CTPackingVolumeFraction;
    static double CTAirVolumeFraction;

    static double CTDropletRadius;
    static double CTDropletVolume;

    static double CTDropletSurfaceArea;

    static double CTLewisFactor;
    static double CTCpAir;

    static double CTTransferCoefCoef;
    static double CTDefaultMassTransferCoefficientAir;
    static double CTDefaultHeatTransferCoefficientWater;
    static double CTDefaultHeatTransferCoefficientAir;

    static int CTNIn;
    static int CTNOut;

    static double AAntoineWater;
    static double BAntoineWater;
    static double CAntoineWater;
    static double AbsHumidityConst;
    static double ConvertmmHgtoPa;

    static double WaterSatPressC1;
    static double WaterSatPressC2;
    static double WaterSatPressC3;
    static double WaterSatPressC4;
    static double WaterSatPressC5;
    static double WaterSatPressC6;

    static double BuckC1;
    static double BuckC2;
    static double BuckC3;
    static double BuckC4;

    static double CTTuningFactor;

    static double CTHeightDraw;
    static double CTWidthDraw;
    static double CTInPointsFraction[];
    static double CTOutPointsFraction[];

    static double CTTemperatureTau;

            //Strm1 is normally the warm stream, and Strm2 the cold stream.  So for the CT strm1 will then be the water that is coming in, and strm2 the air.
    static double CTDefaultU;
    static double CTDefaultA;

    static double CTMassFlowStrm0T0;
    static double CTMolFlowStrm0T0;

    static double CTMassFlowStrm1T0;
    static double CTMolFlowStrm1T0 ;

    static double CTPStrm0Inlet;
    static double CTPStrm0Outlet;
    static double CTPStrm1Inlet;
    static double CTPStrm1Outlet;

    static double CTTStrm0Inlet;
    static double CTTStrm0Outlet;
    static double CTTStrm1Inlet;
    static double CTTStrm1Outlet;

    static double CTTInterfaceT0;

    static double CTStrm0ValveOpeningDefault;
    static double CTStrm0Cv;

    static double CTPMaxFactorIncreaseperSampleT;

    static double CTStrm1FlowCoefficient;
    static double CTStrm0TempTau;
    static double CTStrm1TempTau;
    static double CTStrm0FlowTau;
    static double CTStrm1FlowTau;

    //coolingtowerheatexchangersimple class constants (CTHES) -------------------------------------------------------------------------------------

    static int CTHESNIn;
    static int CTHESNOut;

    static double CTHESTuningFactor;
    static double CTHESHeight;
    static double CTHESWidth;
    static double CTHESInPointsFraction[];
    static double CTHESOutPointsFraction[];

    static double CTHESTemperatureTau;

            //Strm1 is normally the warm stream, and Strm2 the cold stream.  So for the CT strm1 will then be the water that is coming in, and strm2 the air.
    static double CTHESDefaultU;
    static double CTHESDefaultA;

    static double CTHESMassFlowStrm1T0;
    static double CTHESMolFlowStrm1T0;
    static double CTHESMassFlowStrm2T0;
    static double CTHESMolFlowStrm2T0;

    static double CTHESPStrm1Inlet;
    static double CTHESPStrm1Outlet;
    static double CTHESPStrm2Inlet;
    static double CTHESPStrm2Outlet;

    static double CTHESTStrm1Outlet;
    static double CTHESTStrm2Outlet;

    static double CTHESStrm1FlowCoefficient;

    static double CTHESStrm2FlowCoefficient;
    static double CTHESStrm1TempTau;
    static double CTHESStrm2TempTau;
    static double CTHESStrm1FlowTau;
    static double CTHESStrm2FlowTau;

    //coolingtowersimple class constants (CTS) ---------------------------------------------------------------------------------------------------
    static double CTSApproach;
    static int CoolingTowerSimpleNIn;
    static int CoolingTowerSimpleNOut;
    static double CTSFlowMakeUp;
    static double CTSVaporisationFraction;
    static double CTSTuningFactor;

    static double CTSHeight;
    static double CTSWidth;
    static double CTSInPointsFraction[];
    static double CTSOutPointsFraction[];

    static double CTSTemperatureTau;

    //distillation column default properties -------------------------------------------------------------------------------------------------------
    static double DistillationColumnRadius;
    static double DistillationColumnHeight;
    static int NTrays;
    static int DistillationColumnNIn;
    static int DistillationColumnNOut;
    static double InitialDCTrayVolume;
            //Chemical Process Equipment - Selection and Design
    static double InitialDCTrayU;
    static double InitialDCTrayn;

    static double DistillationColumnInPointsFraction[];
    static double DistillationColumnOutPointsFraction[];

            //embeddedtrend class constants
    static double EmbeddedTrendWidth;
    static double EmbeddedTrendHeight;

            //flange class constants -----------------------------------------------------------------------------------------------------------------------
    static double FlangeLength;

    //ftreactor class constants
    static double FTReactorRadius;
    static double FTReactorHeight;
    static int FTReactorNIn;
    static int FTReactorNOut;
    static double FTReactorInPointsFraction[];
    static double FTReactorOutPointsFraction[];

    static double FTReactorMaxVolume;
    static double FTReactorInitInventory;
    static int OrigLoading;
    static double CatDecayRate;
    static double FreshCatAct;
    static double RegenCatAct;
    static double NrRegen;
    static double LostInRegen;
    static double CatTake;
    static double RegenIn;
    static double RegenLosses;
    static double CatSpent;
    static double FreshCatIn;
    static double FreshCatLoadingConst[];

    //gaspipe class constants
    static double PipeDefaultLength;
    static double PipeDefaultDiameter;
    static double PipeDefaultPressure;
    static double PipeDefaultTemperature;
            //public static double PipeDefaultMoles = 100; //This will need to be changed.
    static double PipeDefaultFiLocation;

            //liquidpipe class constants

    static std::string liquidpipeflowreferencestrings[];

    //heatexchanger class ---------------------------------------------------------------------------------------------------
    static double HEThermalPower;
    static int HeatExchangerNIn;
    static int HeatExchangerNOut;
    static double HeatExchangerInPointsFraction[];
    static double HeatExchangerOutPointsFraction[];
    static double HeatExchangerRadius;
    static double HeatExchangerWidth;

    static int NStrm2Coils;
    static int HENSegments;
    static int HENNodes;

            //This is based on the Modelling design Excel file and Areva design.

    static double HETStrm1Inlet;
    static double HETStrm1Outlet;
    static double HETStrm2Inlet;
    static double HETStrm2Outlet;

    static double HEPStrm1Inlet;
    static double HEPStrm1Outlet;
    static double HEPStrm2Inlet;
    static double HEPStrm2Outlet;

    static double HEPStrm1Delta;
            //outflow[0] is now going to
            //just be an extension of the final
            //segment
    static double HEPStrm2Delta;

    static double HETStrm1T0[];
    static double HETStrm2T0[];
    static double HEPStrm1T0[];
    static double HEPStrm2T0[];

    static double HEMassFlowStrm1T0;
    static double HEMassFlowArrayStrm1T0[];
    static double MolFlowStrm1T0;

    static double HEMassFlowStrm2T0;
    static double MolFlowStrm2T0;

            //heatexchanger class : Metal differential equation constants in particular
    static double HEM;

            //heatexchanger class : pressure drop / energy drop due to friction constants
            //public static double HEAddFriction = 10.0;
    static double HEStrm1AddFriction;
    static double HEStrm2AddFriction;
    static double HEStrm1DeltaPK[];

    static double HEStrm2DeltaPK[];
            //public static double[] HEStrm2DeltaPK = HEPStrm2Delta / Math.Pow(MolFlowStrm2T0, 2.0);

            //heatexchanger class : heat exchange constants
    static double HEHeatExchangeSurfaceArea;
    static double HEOutsideDiameterTube;
    static double HETubeWallThickness;
    static double HEInsideDiameterTube;
    static double HENrPassesThroughShell;
    static double HETubeCircOutside;
    static double HETubeCircInside;
    static double HETubeCircAve;
    static double HEAveLengthPerTube;
    static double HEAveLengthPerSegment;
    static double HEAStrm2;
    static double HEStrm2TubeVolume;
    static double HEStrm2SegmentVolume;
    static double HEShellVolume;
    static double HEStrm1TubeVolume;
    static double HEAStrm1;
    static double HEEffGasTubeCircInside;
    static double HEStrm1SegmentVolume;
    static double HERi;

            //public static double HEAg = 0.035476792;  //m2;  From Excel sheet from AREVA design.
    static double HEInPointsFraction[];
    static double HEOutPointsFraction[];
            //public static double HETubeDiameter = 0.023; //m; Diameter of the tubes in the steam generator.  From Excel sheet
            //public static double HEAs = Math.PI * Math.Pow(HETubeDiameter / 2.0, 2.0); //m2; Cross sectional area of the pipes for the steam/water

            //thermal resistivity of the metal times the tube thickness.
            //This kgm and Ri, needs to be backed up with some more science.  Why is it so
            //dificult to get these values and to fix things up properly?
    static double HECsi;
    static double HECgi;
    static double HEHeatTransferArea;
    static double HEKgm[];
    static double HEKms[];
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
    static double HeatExchangerSimpleDefaultU;
    static double HeatExchangerSimpleDefaultA;

    static double HESMassFlowStrm1T0;
    static double HESMolFlowStrm1T0;
    static double HESMassFlowStrm2T0;
    static double HESMolFlowStrm2T0;

    static double HESPStrm1Inlet;
    static double HESPStrm1Outlet;
    static double HESPStrm2Inlet;
    static double HESPStrm2Outlet;

    static double HESStrm1FlowCoefficient;
    static double HESStrm2FlowCoefficient;
    static double HESStrm1TempTau;
    static double HESStrm2TempTau;
    static double HESStrm1FlowTau;
    static double HESStrm2FlowTau;

    static int HESNSegments;
    static int HESNStrm2Coils;

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
    static double HESStrm2Volume;
    static double HESStrm2SegmentVolume;
    static double HESShellVolume;
    static double HESStrm1TubeVolume;
            //public static double HEAStrm1 = HEStrm1TubeVolume / HEAveLengthPerTube;  //m2;  From Excel sheet from AREVA design.
            //public static double HEEffGasTubeCircInside = 2 * Math.Sqrt(HEAStrm1 / Math.PI) * Math.PI;//2 * Math.Sqrt(HEAStrm1 / Math.PI) * Math.PI;
    static double HESStrm1SegmentVolume;

    //material class constants ---------------------------------------------------------------------------------------------------------------------
    static double Udefault;
    static double Vdefault;
    static double fdefault;
    static std::vector<component> fluidpackage;
    static materialphase MaterialInitPhase;
    static int NMaterialIterations;
    static double ZNotDefined;
    static double epsilonadd;
    static double epsilonfrac;
    static double Inf;

    //mixer class constants ------------------------------------------------------------------------------------------------------------------------
    static double MixerLength;
    static double MixerDistanceBetweenBranches;
    static int MixerDefaultNIn;
    static double MixerBranchThickness;
    static double MixerInitRadiusDefault;

    //molecule class constants
    static double InitTc;
    static double InitPc;
    static double Initomega;

    //nmpc class constants ----------------------------------------------------------------------------------------------------------------------
    static double NMPCWidth;
    static double NMPCHeight;
    static int DefaultN;
    static int DefaultInitialDelay;
    static int DefaultRunInterval;
                                                        //that.
    static double Defaultalphak;
    static nmpcalgorithm DefaultNMPCAlgorithm;
    static double DefaultNMPCSigma;

            //Interior Point 1 - ConstrainedLineSearch algorithm constants
    static double DefaultMuBarrier;
    static double Defaultsvalue;
                                                            //the initial value of the valve openings, and thus the values will be 0.5 away from zero.
    static double DefaultsvalueMultiplier;
            //public static double DefaultSigma = 0.1; //The fraction multiplied with mubarrier at the end of each iteration.
    static double DefaulttauIP;
    static double DefaultIPErrorTol;
    static double CholeskyDelta;
    static double CholeskyBeta;
    static double MVMaxMovePerSampleTT0;
    static double NMPCIPWeightPreTerm;
    static double NMPCIPWeightTerminal;

            //Genetic Algorithm 1 constants (mostly default values for variables in the nmpc class).
    static int DefaultNrChromosomes;
    static int DefaultNrSurvivingChromosomes;
    static int DefaultNrPairingChromosomes;
    static double DefaultProbabilityOfMutation;
    static int DefaultNrIterations;
    static int DefaultCrossOverPoint;

            //PSO constants
    static int DefaultNrContinuousParticles;
    static int DefaultNrBooleanParticles;
    static int PSOMVBoundaryBuffer;
    static double PSOMaxBooleanSpeed;

    //pidcontroller class constants -------------------------------------------------------------------------------------------------------------
    static int Direct;
    static int Reverse;;
    static double PIDControllerInitRadius;
    static double PIDControllerInitK;
    static double PIDControllerInitI;
    static double PIDControllerInitD;
    static double PIDControllerInitMinOP;
    static double PIDControllerInitMaxOP;
    static double PIDControllerInitMinPV;
    static double PIDControllerInitMaxPV;

    //pump class default properties -------------------------------------------------------------------------------------------------------------
    static double PumpInitMaxDeltaPressure;
    static double PumpInitMinDeltaPressure;
    static double PumpInitMaxActualFlow;
    static double PumpMinActualFlow;
    static double PumpCurveYAxis;
    static double PumpCurvef1;
    static double PumpCurvep1;
    static double PumpCurvef2;
    static double PumpCurveSpeedT0;
    static double PumpSpeedTau;

    static double PumpInitActualVolumeFlow;
    static double PumpInitOn;
    static double PumpInitRadius;
    static double PumpInitOutletLength;
    static double PumpInitOutletRadius;

    //stream class default properties -------------------------------------------------------------------------------------------------------------
    static int SignalNrPropDisplay;

            //steamgenerator class constants
    static double SteamGeneratorRadius;
    static double SteamGeneratorHeight;
    static int SteamGeneratorNSegments;
    static int SteamGeneratorNNodes;
    static double SteamGeneratorWaterTubeVolume;
    static double SteamGeneratorWaterSegmentVolume;
    static double SteamGeneratorGasTubeVolume;
    static double SteamGeneratorGasSegmentVolume;
    static int SteamGeneratorNIn;
    static int SteamGeneratorNOut;
    static int NSteamCoils;
    static double ThermalPower;
    static double TubeCircOutside;
            //This is based on the Modelling design Excel file and Areva design.
    static double TubeCircInside;
    static double TubeCircAve;
    static double AveLengthPerTube;
    static double AveLengthPerSegment;
    static double SteamGeneratorInPointsFraction[];
    static double SteamGeneratorOutPointsFraction[];

    static double Ws0T0;

    static double SGTGasInlet;
    static double SGTGasOutlet;
    static double SGTWaterInlet;
    static double SGTWaterOutlet;
    static double WgasT0;
    static double MolFlowGasT0;
    static double WwaterT0;
    static double MolFlowWaterT0;

    static double SGPGasDelta;
            //outflow[0] is now going to
            //just be an extension of the final
            //segment
    static double SGPWaterInlet;
    static double SGPWaterOutlet;
    static double SGPWaterDelta;

    static double TGasT0[];
    static double TWaterT0[];
    static double PGasT0[];
    static double PWaterT0[];

            //steamgenerator class : Heat exchange constants
    static double Ri;
    static double Csi;
    static double Cgi;
    static double HeatTransferArea;

    static double Kgm[];

            //public static double Kms =
            //    (HeatTransferArea * (0.5 * (0.5 * (SGTGasInlet + SGTGasOutlet) - 0.5 * (SGTWaterOutlet + SGTWaterInlet))) /
            //    (ThermalPower / NSteamCoils / SteamGeneratorNSegments) - Ri * 0.5) /
            //    Math.Exp(-WwaterT0 * Cgi);  //THIS IS TO BE PUT BACK AS THE WAY OF INITIALISING LATER ON.
    static double Kms[];

    static double WgT0[];
    static double Ag;
    static double EffGasTubeCircInside;
    static double TubeDiameter;
    static double As;

            //steamgenerator class : Metal differential equation constants in particular
    static double M;

            //steamgenerator class : pressure drop / energy drop due to friction constants
    static double SteamGenAddFriction;
    static double SteamGenGasAddFriction;
    static double SGWaterDeltaPK;
    static double SGGasDeltaPK;

            //steamgenerator class : MOMENTUM EQUATION constants
    static double DynViscA;
    static double DynViscB;

    //stream class default properties -------------------------------------------------------------------------------------------------------------

    static double MinDistanceFromStream;
    static double StreamArrowAngle;
    static double StreamArrowLength;
    static double StreamMaxMassFlow;
    static double StreamMinMassFlow;
    static int StreamNrPropDisplay;

    //tank class default properties ----------------------------------------------------------------------------------------------------------------
    static double TankInitRadius;
    static double TankInitHeight;
    static double TankRadiusDraw;
    static double TankHeightDraw;
    static double TankInitMaxVolume;
    static double TankInitFracInventory;
    static double TankInitInOutletDistanceFraction;
    static double TankMinFracInventory;

    //tee class constants --------------------------------------------------------------------------------------------------------------------------
    static double TeeLength;
    static double TeeDistanceBetweenBranches;
    static int TeeDefaultNOut;
    static double TeeBranchThickness;
    static double TeeInitRadiusDefault;

            //trend class constants -----------------------------------------------------------------------------------------------------------------------
    static float XAxisMargin;
    static float YAxisMargin;
    static double YIncrement;
    static double YmaxMaxFactorAbove;

    //valve class constants -----------------------------------------------------------------------------------------------------------------------
    static double ValveDefaultActualFlow;
    static double ValveDefaultDP;
    static double ValveDefaultOpening;

    static double ValveEqualPercR;
    static double ValveDefaultCv;
    static double ValveHydraulicTau;

    static double ValveLength;
    static double ValveWidth;


    global();
    ~global();

    //void initsimtimevector();


};


#endif // GLOBAL_H
