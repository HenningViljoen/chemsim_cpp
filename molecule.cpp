#include "molecule.h"
#include <math.h>


molecule::molecule(std::string anabreviation, std::string aname, double amolarmass, double adynamicviscosity, double adensity,
                   double aTc, double aPc,
            double aomega, double aCpA, double aCpB, double aCpC, double aCpD)
        //Propane is the default Cp coefficients that have been inserted here.  The Cp arguments here are heat capacity coefficients for the molecule.
{
    abreviation = anabreviation;
    name = aname;
    molarmass = amolarmass;
    dynamicviscosity = adynamicviscosity;
    density = adensity;

    Tc = aTc;
    Pc = aPc;
    omega = aomega;
    CpA = aCpA; //coeficient in equation to calculate Cp is a function of T
    CpB = aCpB; //coeficient in equation to calculate Cp is a function of T
    CpC = aCpC; //coeficient in equation to calculate Cp is a function of T
    CpD = aCpD; //coeficient in equation to calculate Cp is a function of T
}

double molecule::calcCp(double T) //given the Temperature, calculate the Cp for the molecule.
{
    return CpA + CpB * T + CpC * pow(T, 2.0) + CpD * pow(T, 3.0);
}
