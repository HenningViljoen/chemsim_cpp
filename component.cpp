#include "component.h"

component::component(molecule *am, double amolefraction, double an)
{
    molefraction = amolefraction;
    n = an;
    m = am;
    //massfraction = amassfraction;
}

component::~component()
{
    //no deletes needed as there were no news.
}

void component::copytothisobject(component *c)
{
    m = c->m; //Do not make a copy of the c.m object, just point to the address of c.m, since c.m should never change between unit ops - it is part
    // of the fluid package.
    molefraction = c->molefraction;
    n = c->n;
    massfraction = c->massfraction;
    //molefraction = c.molefraction;
}
