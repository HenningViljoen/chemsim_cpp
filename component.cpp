#include "component.h"

component::component()
{
    molefraction = 0.0;
    n = 0.0;
    m = nullptr;
}

component::component(molecule *am, double amolefraction, double an)
{
    molefraction = amolefraction;
    n = an;
    m = am;
    //massfraction = amassfraction;
}

component::component(const component &c)
{
    copytothisobject(c);
}

component::~component()
{
    //if (molecule != nullptr) {delete molecule;}  I do not want to delete the molecule since it is allocated memory once in the fluidpackage
        //and not again

}

component &component::operator=(const component &c)
{
    component comp(c);
    return comp;
}

void component::copytothisobject(const component &c)
{
    m = c.m; //Do not make a copy of the c.m object, just point to the address of c.m, since c.m should never change between unit ops - it is part
    // of the fluid package.
    molefraction = c.molefraction;
    n = c.n;
    massfraction = c.massfraction;
    //molefraction = c.molefraction;
}
