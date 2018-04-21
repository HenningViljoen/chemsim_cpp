#ifndef BASECLASS_H
#define BASECLASS_H

#include <string>
#include <QObject>
#include <QPainter>
#include <QPoint>
#include <QPolygon>
#include <vector>

#include "point.h"
#include "global.h"
#include "controlvar.h"
//#include "simulation.h"

class simulation;

class baseclass
{
public:
    int nr; //Index of this object in the collection that it is part of.
    std::string name; //Name of the piece of equipment / stream that is under scope.
            //public bool on; //Is the equipment/unitop/stream/controller, on, or off.
    objecttypes objecttype; //What kind of unit op is a particular one.
    bool highlighted;
    point location;
    std::vector<std::string> controlproperties;
    std::vector<std::string> controlpropthisclass;
    int nrcontrolpropinherited;
    std::vector<point> points; //Points used for the item on the PFD.  Expressed in meters.

    baseclass(int anr, double ax, double ay);
    ~baseclass();

    void copyfrom(baseclass *baseclasscopyfrom);
    virtual controlvar *selectedproperty(int selection);
    virtual void update(int i, bool historise);

    virtual void updatepoint(int i, double x, double y);
    virtual bool mouseover(double x, double y); //This function will indicate whether the mouse is over a particular unitop or stream at any moment in time.
    virtual void setproperties(simulation *asim); //Method that will be inherited and that will set the properties of the applicable object in a window
    //public virtual void showtrenddetail(simulation asim, List<Form> detailtrendslist) //Virtual method that will set up the trend detail window for the
    virtual void draw(QImage *graphics);
};

#endif // BASECLASS_H
