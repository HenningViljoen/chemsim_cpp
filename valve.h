#ifndef VALVE_H
#define VALVE_H

#include "unitop.h"
#include "drawwidget.h"

class valve : public unitop
{
public:
    valve(int anr, double ax, double ay, double aCv, double aop);
    //void valve()

    void draw();
};

#endif // VALVE_H
