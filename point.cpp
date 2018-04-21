#include "point.h"

point::point()
{
    setxy(0.0, 0.0);
    highlighted = false;
}

point::point(double ax, double ay)
{
    setxy(ax, ay);
    highlighted = false;
}

point::point(const point &pointcopyfrom)
{
    copyfrom(pointcopyfrom);
}

point &point::operator=(const point &p)
{
    point pnt(p);
    return pnt;
}

void point::copyfrom(const point &pointcopyfrom)
{
    setxy(pointcopyfrom.x, pointcopyfrom.y);
    highlighted = pointcopyfrom.highlighted;
}

void point::setxy(double ax, double ay)
{
    x = ax;
    y = ay;
}
