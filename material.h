#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>
#include <math.h>
#include <algorithm>
#include "controlvar.h"
#include "matrix.h"
#include "complex.h"
#include "global.h"
#include "utilities.h"

class material
{
public:
    controlvar *P = nullptr; //Pa

    controlvar *V = nullptr; //m3.  This will be total stage volume in the case of a distillation column.
    controlvar vmolar; //m3/mole.  Total volume per mole for the total material all phases.
    controlvar vmolarL; //m3/mole.  Total volume per mole for the total liquid phase.
    controlvar vmolarV; //m3/mole.  Total volume per mole for the total vapour phase.
    controlvar *density = nullptr; //kg/m3

    controlvar *T = nullptr; //Kelvin

    controlvar *mass = nullptr; //kg
    controlvar *n = nullptr; //amount of moles of the material (all phases) for all components.

    controlvar f; //fraction.  Vapour molar fraction.

    controlvar U; //Joule. Internal energy.

    controlvar relativehumidity; //% ; Using relative humidity at the moment since that is the data that we have on site.


    std::vector<double> x; //Molar fraction of the total liquid of the material per component.
    std::vector<double> y; //Molar fraction of the total vapour of the material per component.
    double hL; //Joule/mol.  Molar enthalpy for the liquid phase, e.g. on a stage in a distillation column.
    double hV; //Joule/mol.  Molar enthalpy for the vapour phase, e.g. on a stage in a distillation column.
    double ML; //moles.  Total molar liquid hold-up.
    double MV; //moles.  Total molar vapour hold-up.
    std::vector<double> z; //fraction.  Molar fraction per component of the total moles in the material.

    controlvar umolar; //Joule/mol.  Molar internal energy.
    std::vector<double> umolarL; //Joule/mol.  Molar internal energy for the total liquid phase.
    std::vector<double> umolarV; //Joule/mol.  Molar internal energy for the total vapour phase.
    double totalumolar; //last calculation of the total vmolar;
    std::vector<double> umolarideal; //Joule/mol.  Molar internal energy for the ideal gas case.
    std::vector<double> fugacityL; //fugacity per component for the liquid phase.
    std::vector<double> fugacityV; //fugacity per component for the vapour phase.
    std::vector<double> aT; //kg⋅m^5 /( s^2 ⋅mol^2 ); figure used in Peng Robinson equation.
    std::vector<double> ac; // kg⋅m^5 / (s^2 ⋅mol^2 ) ; figure used in Peng Robinson equation.
    std::vector<double> alpha; //part of the equation for aT, and used in other equations in the Peng Robinson collection.
    std::vector<double> b;  // m^3/mol ; figure used in Peng Robinson equation.
    std::vector<double> K;  //dimensionless - depends on omega.
    std::vector<double> A;  //used in equation for compressibility factor.
    std::vector<double> B;  //used in equation for compressibility factor.
    std::vector<std::vector<double>> Z;  //Compressibility factor. More than one value should be present based on the amount of phases in the material
    std::vector<double> discriminantZ; //The discriminant of the cubig PR equation for Z;
    std::vector<double> Tr; //reduced temperature.
    std::vector<double> acomp; //Cubic equation for compressibility factor coeficients.
    std::vector<double> bcomp; //Cubic equation for compressibility factor coeficients.
    std::vector<double> ccomp; //Cubic equation for compressibility factor coeficients.
    std::vector<double> dcomp; //Cubic equation for compressibility factor coeficients.
    double bs; //Used in scarlet cubic solver. bs for bscarlet, This is for another cubic solver: http://home.scarlet.be/~ping1339/cubic.htm
    double cs; //Used in scarlet cubic solver.
    double ds; //Used in scarlet cubic solver.
    double rs; // Used in scarlet cubic solver.
    double es; //Used in scarlet cubic solver.
    double fs; //Used in scarlet cubic solver.
    std::vector<complex> zs; //Used in scarlet cubic solver.
    std::vector<complex> us; //Used in scarlet cubic solver.
    std::vector<complex> ys; //Used in scarlet cubic solver.

    std::vector<double> delta0; //used to solve the cubic equation for Z.
    std::vector<double> delta1; //used to solve the cubic equation for Z.
    std::vector<double> C; //used to solve the cubic equation for Z.
    std::vector<complex> Cc; //used to solve the cubic equation for Z.
    std::vector<std::vector<double>> a; //used to solve the cubid equation for Z by making use of an alternative closed form formulat for them.
    std::vector<double> Qc; //used in alternative formula.
    std::vector<double> Rc; //used in alternative formula.
    std::vector<double> Dc; //used in alternative formula.
    std::vector<complex> Sc; //used in alternative formula.
    std::vector<complex> Tc; //used in alternative formula.
    std::vector<double> Cp; //Heat capacity of constant pressure for an ideal gas type scenario per component.

    double totalCp; //Heat capacity of the material in totallity, with all its components includded.
    std::vector<double> xvector;
    matrix jacobian; //The jacobian matrix for the system of equations as defined below in fflash().
    int uvflashsize; //The size of the k variables k equations problem to be solved for the UV flash.
    int origuvflashsize;

    std::vector<component> composition; //mass in kg for each component in the material
    double massofonemole; //kg/mole
    double volumeofonemole; //m^3/mole
    materialphase phase;

    material(double aV);
    material(const material &m); //copy constructor
    material(std::vector<component> &acomposition, double aTemp, double aV, double aP, double af);
    ~material();

    material &operator=(const material &m);

    void copyfrom(const material &materialcopyfrom);
    void copycompositiontothismat(material &amat);
    static void copycompositiontothiscomposition(std::vector<component> &compositioncopyto, std::vector<component> &compositioncopyfrom);
    void init(std::vector<component> &acomposition, double aTemp, double aV, double aP, double af);
    void setxybasedonf();
    void PTfVflash(double aTemp, double aV, double aP, double af);
    void nullZ(int j);
    void calcmass();
    component &match(std::string aname);
    void copycompositiontothisobject(material &m);
    void mapvarstox();
    void mapxtovars();
    complex &calcxk(complex &u, int j);
    void calcZ(int j);
    void calcfugacity(int j);
    double calcumolarwithZ(int j, int Zi);
    void calcumolarpercomponent();
    void calctotalCp();
    double calcumolar();
    double calcvmolaranddensity();
    double dfdx(int fi, int xi);
    double fflash(int i);
    void calcjacobian();
    void limitxvector();
    void calccompz();
    void uvflash();
    void zero();
    void update(int simi, bool historise);
};

#endif // MATERIAL_H
