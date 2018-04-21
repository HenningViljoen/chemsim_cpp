#include "utilities.h"

utilities::utilities()
{

}

double dndt2fps(double adndt)
{
    return adndt * (global::R * global::Ts) / global::Ps;
}

double utilities::calcdirection(double deltay, double deltax)
{
     double bias;
     if (deltax == 0) { deltax = 0.0001; }
     if (deltax < 0)
     {
          bias = M_PI;
     }
     else
     {
          bias = 0;
     }
     return bias + atan(deltay / deltax);
}

double utilities::distance(point &p1, point &p2)
{
     return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

double utilities::distance(double deltax, double deltay)
{
     return sqrt(pow(deltax, 2) + pow(deltay, 2));
}

