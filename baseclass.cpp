#include "baseclass.h"

baseclass::baseclass(int anr, double ax, double ay) :
    nr {anr},
    name {""},
    location(ax, ay),
    controlproperties {},
    controlpropthisclass {},
    nrcontrolpropinherited {0},
    highlighted {false}

{
    //on = true;
}

baseclass::~baseclass()
{

}

void baseclass::copyfrom(baseclass *baseclasscopyfrom)
{
    nr = baseclasscopyfrom->nr;
    name = baseclasscopyfrom->name;
    //on = baseclasscopyfrom.on;
    objecttype = baseclasscopyfrom->objecttype;
    highlighted = baseclasscopyfrom->highlighted;
    location.copyfrom(baseclasscopyfrom->location);
    controlproperties = baseclasscopyfrom->controlproperties;
    controlpropthisclass = baseclasscopyfrom->controlpropthisclass;
    nrcontrolpropinherited = baseclasscopyfrom->nrcontrolpropinherited;
}

controlvar *baseclass::selectedproperty(int selection)
{
    return nullptr;
}

void baseclass::update(int i, bool historise) //index for where in the simvectors the update is to be stored,
                                                            //boolean for whether or not history should be kept for the simulation
                                                            //at this time.
{
}

void baseclass::updatepoint(int i, double x, double y)
{
}

bool baseclass::mouseover(double x, double y) //This function will indicate whether the mouse is over a particular unitop or stream at any moment in time.
{
     return false;
}

void baseclass::setproperties(simulation *asim) //Method that will be inherited and that will set the properties of the applicable object in a window
{
}

//void baseclass::showtrenddetail(simulation *asim, List<Form> detailtrendslist) //Virtual method that will set up the trend detail window for the
//        //                                                     //applicable object.
//        {
//        }


void baseclass::draw(QImage *graphics)
{
}
