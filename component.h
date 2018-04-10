#ifndef COMPONENT_H
#define COMPONENT_H

#include "molecule.h"

class component
{
public:
    molecule *m;
    double n; //Amount of moles of this molecule/compound that is in this material.
    double molefraction; //fraction of total moles in the material that is from this component.
    double massfraction; //mass fraction of this component of the material in a particular unitop or stream.
            //public double molefraction;

    component(molecule *am = nullptr, double amolefraction = 0.0, double an = 0.0);
    ~component();
    void copytothisobject(component *c);
};

#endif // COMPONENT_H
