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

molecule::molecule()
{

}

molecule::molecule(const molecule &amolecule)
{
    abreviation = amolecule.abreviation;
    name = amolecule.name;
    molarmass = amolecule.molarmass; //kg / mole
    dynamicviscosity = amolecule.dynamicviscosity; //Pa·s  We need the dynamic viscosity here since it is used to calculate the
            //Renoulds number in liquid flow applications.  Later temperature dependance will be added.  Standard conditions for now.
    density = amolecule.density;          //kg/m3.  Also used in Re nr calc.  Later temp dependance will be added.  Standard conditions for now.
            //public double defaultmolefraction; //fraction of this molecule that is the default fraction for each stream.
    Tc = amolecule.Tc; //Kelvin; Critical temperature.
    Pc = amolecule.Pc; //Pascal; Critical pressure.
    omega = amolecule.omega; //unitless; Acentric factor.
    CpA = amolecule.CpA; //coeficient in equation to calculate Cp is a function of T
    CpB = amolecule.CpB; //coeficient in equation to calculate Cp is a function of T
    CpC = amolecule.CpC; //coeficient in equation to calculate Cp is a function of T
    CpD = amolecule.CpD; //coeficient in equation to calculate Cp is a function of T
}

molecule &molecule::operator=(const molecule &amolecule)
{
    molecule mol(amolecule);
    return mol;
}

double molecule::calcCp(double T) //given the Temperature, calculate the Cp for the molecule.
{
    return CpA + CpB * T + CpC * pow(T, 2.0) + CpD * pow(T, 3.0);
}
