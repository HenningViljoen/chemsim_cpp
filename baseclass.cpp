#include "baseclass.h"

baseclass::baseclass(int anr, double ax, double ay)
{
    location = new point(ax, ay);
}

baseclass::~baseclass()
{
    delete location;
}

void baseclass::draw(QImage *graphics)
{

}
