#ifndef VALVE_H
#define VALVE_H

#include "unitop.h"
#include "drawwidget.h"


class valve : public unitop
{
public:
    //variables
    controlvar *deltapressure = nullptr; //Pa
    //public double[] deltapressuresimvector; //Pa; history sim vector.
    double Cv; //Valve coeficient
    controlvar *op = nullptr; //Fraction : valve opening as a fraction
    //public double[] opsimvector;

    double deltapressurenew; //Pa
    double ddeltapressuredt; //Pa/s

    valve(int anr, double ax, double ay, double aCv, double aop);
    valve(baseclass *baseclasscopyfrom);
    ~valve();

    void draw(QImage *graphics);
    void initvalve(int anr, double ax, double ay, double aCv, double aop);
    void copyfrom(baseclass *baseclasscopyfrom);
    controlvar *selectedproperty(int selection);
    void ddt(int simi);
    void update(int simi, bool historise);
    void sizevalvefromactualvolumeflow();
    bool mouseover(double x, double y);
    void updateinoutpointlocations();
    void setproperties(simulation *asim);
};

#endif // VALVE_H
