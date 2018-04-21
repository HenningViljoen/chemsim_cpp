#include "controlvar.h"



controlvar::controlvar(double av, bool aisbool)
{
    v = av;
    isbool = aisbool;
    //simvector = new double[global.SimIterations];
    excelsource = nullptr;
    datasource = Simulation;
}

controlvar::~controlvar()
{
    if (excelsource != nullptr) {delete excelsource;}
}

controlvar::controlvar(const controlvar &copyfrom)
{
     v = copyfrom.v;
     excelsource = copyfrom.excelsource;
     datasource = copyfrom.datasource;
}

controlvar &controlvar::operator=(const controlvar &cv)
{
    controlvar cont(cv);
    return cont;
}

void controlvar::copyfrom(controlvar *copyfrom)
{
     v = copyfrom->v;
     excelsource = copyfrom->excelsource;
     datasource = copyfrom->datasource;
}

std::string controlvar::ToString(std::string format)
{
     return std::to_string(v);
}



//Operator functions which are not part of the class:

controlvar &operator-(controlvar &c, double d)
{
    controlvar cvar(c.v - d);
    return cvar;
}

controlvar &operator-(controlvar &c)
{
    controlvar cvar(-c.v);
    return cvar;
}

bool operator<(controlvar &c, double d)
{
    return c.v < d;
}

bool operator>(controlvar &c, double d)
{
    return c.v > d;
}

controlvar &operator/(controlvar &c, double d)
{
    controlvar cvar(c.v / d);
    return cvar;
}

controlvar &operator*(controlvar &c, double d)
{
    controlvar cvar(c.v * d);
    return cvar;
}

controlvar &operator+(controlvar &c, double d)
{
    controlvar cvar(c.v + d);
    return cvar;
}

controlvar &operator*(double d, controlvar &c)
{
    controlvar cvar(d * c.v);
    return cvar;
}

controlvar &operator/(double d, controlvar &c)
{
    controlvar cvar(d / c.v);
    return cvar;
}

controlvar &operator*(controlvar &c1, controlvar &c2)
{
    controlvar cvar(c1.v * c2.v);
    return cvar;
}

controlvar &operator/(controlvar &c1, controlvar &c2)
{
    controlvar cvar(c1.v/c2.v);
    return cvar;
}

controlvar &operator+(controlvar &c1, controlvar &c2)
{
    controlvar cvar(c1.v + c2.v);
    return cvar;
}

controlvar &operator-(controlvar &v1, controlvar &v2)
{
    controlvar cvar(v1.v - v2.v);
    return cvar;
}

controlvar &Pow(controlvar &c, double y)
{
    controlvar cvar(pow(c.v, y));
    return cvar;
}

















