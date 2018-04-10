#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>


class molecule
{
public:
    std::string abreviation;
    std::string name;
    double molarmass; //kg / mole
    double dynamicviscosity; //Pa·s  We need the dynamic viscosity here since it is used to calculate the
            //Renoulds number in liquid flow applications.  Later temperature dependance will be added.  Standard conditions for now.
    double density;          //kg/m3.  Also used in Re nr calc.  Later temp dependance will be added.  Standard conditions for now.
            //public double defaultmolefraction; //fraction of this molecule that is the default fraction for each stream.
    double Tc; //Kelvin; Critical temperature.
    double Pc; //Pascal; Critical pressure.
    double omega; //unitless; Acentric factor.
    double CpA; //coeficient in equation to calculate Cp is a function of T
    double CpB; //coeficient in equation to calculate Cp is a function of T
    double CpC; //coeficient in equation to calculate Cp is a function of T
    double CpD; //coeficient in equation to calculate Cp is a function of T


    molecule(std::string anabreviation, std::string aname, double amolarmass, double adynamicviscosity, double adensity,
             double aTc = 500, double aPc = 50*100000,
      double aomega = 0.3, double aCpA = -4.224, double aCpB = 0.3063, double aCpC = -1.586e-04, double aCpD = 3.215e-08);

    double calcCp(double T);
};

#endif // MOLECULE_H
