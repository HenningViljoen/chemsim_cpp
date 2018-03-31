#include "point.h"

point::point(double ax, double ay)
{
    setxy(ax, ay);
    highlighted = false;
}

point::point(point* pointcopyfrom)
{
    setxy(pointcopyfrom->x, pointcopyfrom->y);
    highlighted = pointcopyfrom->highlighted;

}

void point::setxy(double ax, double ay)
{
    x = ax;
    y = ay;
}
