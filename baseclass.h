#ifndef BASECLASS_H
#define BASECLASS_H

#include <string>
#include <QObject>

#include "point.h"

class baseclass
{
public:
    int nr; //Index of this object in the collection that it is part of.
    std::string name; //Name of the piece of equipment / stream that is under scope.
            //public bool on; //Is the equipment/unitop/stream/controller, on, or off.
    //public objecttypes objecttype; //What kind of unit op is a particular one.
    bool highlighted;
    point *location;
    point *points; //Points used for the item on the PFD.  Expressed in meters.

    baseclass(int anr, double ax, double ay);
    ~baseclass();
    virtual void draw(QImage *graphics);
};

#endif // BASECLASS_H
