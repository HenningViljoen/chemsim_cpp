#ifndef CONTROLVAR_H
#define CONTROLVAR_H

#include <math.h>
#include <vector>
#include "global.h"
#include "exceldataset.h"

class controlvar
{
public:
    double v; //The variable in this object that will have the PV or OP to be controlled in the simulation.
    bool isbool; //Is this a boolean variable?  Could then be part of a hybrid system.
    std::vector<double> simvector;

    exceldataset *excelsource; // for the case that data will be drawn in from an Excel file.
    datasourceforvar datasource; //The source of data for the variable in the model.

    controlvar(double av = 0, bool aisbool = false);
    controlvar(const controlvar &copyfrom);
    ~controlvar();

    controlvar &operator=(const controlvar &cv);
    void copyfrom(controlvar *copyfrom);
    std::string ToString(std::string format = "");


};

controlvar &operator-(controlvar &c, double d);
controlvar &operator-(controlvar &c);
bool operator<(controlvar &c, double d);
bool operator>(controlvar &c, double d);
controlvar &operator/(controlvar &c, double d);
controlvar &operator*(controlvar &c, double d);
controlvar &operator+(controlvar &c, double d);
controlvar &operator*(double d, controlvar &c);
controlvar &operator/(double d, controlvar &c);
controlvar &operator *(controlvar &c1, controlvar &c2);
controlvar &operator/(controlvar &c1, controlvar &c2);
controlvar &operator+(controlvar &c1, controlvar &c2);
controlvar &operator-(controlvar &v1, controlvar &v2);
controlvar &Pow(controlvar &c, double y);

#endif // CONTROLVAR_H
