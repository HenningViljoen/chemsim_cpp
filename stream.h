#ifndef STREAM_H
#define STREAM_H

#include "baseprocessclass.h"

class stream : public baseprocessclass
{
public:
    double direction; //radians
    double distance; //m
    std::vector<point> inbetweenpoints; //points that are between the main in and out points, that will enable straight lines for the stream.
    std::vector<point> displaypoints; //The locations of the T, P and flow properties that are being displayed for the stream.

    stream(int anr, double p0x, double p0y, double p1x, double p1y);
    stream(stream *streamcopyfrom);
    void streaminit(double p0x, double p0y, double p1x, double p1y);
    void copyfrom(baseclass *baseclasscopyfrom);
    void updatedirection();
    void updatemassflowsource(std::string afilename);
    void updatepoint(int i, double x, double y);
    void update(int simi, bool historise);
};

#endif // STREAM_H
