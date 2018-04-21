#ifndef UNITOP_H
#define UNITOP_H

#include "baseprocessclass.h"
#include "stream.h"

class unitop : public baseprocessclass
{
public:
    std::vector<baseprocessclass *> inflow;
    std::vector<baseprocessclass *> outflow;
    std::vector<point> inpoint;   //the points where the stream(s) are coming into the unitop
    std::vector<point> outpoint; //the points where the stream(s) are going out of the unitop.
    int nin, nout; //amount of streams in, and amount of streams out.

    unitop(int anr, double ax, double ay, int anin, int anout);
    unitop(unitop *unitopcopyfrom);
    ~unitop();
    void initunitop(int anin, int anout);
    void copyfrom(baseclass *baseclasscopyfrom);
    virtual void initinflow();
    virtual void initoutflow();
    virtual void initinpoint();
    virtual void initoutpoint();
    virtual void updateinoutpointlocations();
    void update(int i, bool historise);
    void draw(QImage *graphics);
};

#endif // UNITOP_H
