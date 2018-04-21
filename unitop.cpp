#include "unitop.h"

unitop::unitop(int anr, double ax, double ay, int anin, int anout) : baseprocessclass{anr, ax, ay}
{
    initunitop(anin, anout);
}

unitop::unitop(unitop *unitopcopyfrom) : baseprocessclass{unitopcopyfrom->nr, unitopcopyfrom->location.x, unitopcopyfrom->location.y}
{
    initunitop(unitopcopyfrom->nin, unitopcopyfrom->nout);
    copyfrom(unitopcopyfrom);
}

unitop::~unitop()
{
    for (int i = 0; i < inflow.size(); i++)
    {
        if (inflow[i] != nullptr) {delete inflow[i];}
    }
    for (int i = 0; i < outflow.size(); i++)
    {
        if (outflow[i] != nullptr) {delete outflow[i];}
    }
}

void unitop::initunitop(int anin, int anout)
{
    nin = anin;
    nout = anout;
    initinflow();
    initoutflow();
    initinpoint();
    initoutpoint();
}

void unitop::copyfrom(baseclass *baseclasscopyfrom)
{
    unitop *unitopcopyfrom = (unitop *)baseclasscopyfrom;

    baseprocessclass::copyfrom(baseclasscopyfrom);

    //public baseprocessclass[] inflow;
    //public baseprocessclass[] outflow;
    //public point[] inpoint;   //the points where the stream(s) are coming into the unitop
    //public point[] outpoint; //the points where the stream(s) are going out of the unitop.
    //public int nin, nout; //amount of streams in, and amount of streams out.

    for (int i = 0; i < inflow.size(); i++) { inflow[i]->copyfrom(unitopcopyfrom->inflow[i]); }
    for (int i = 0; i < outflow.size(); i++) { outflow[i]->copyfrom(unitopcopyfrom->outflow[i]); }
    for (int i = 0; i < inpoint.size(); i++) { inpoint[i].copyfrom(unitopcopyfrom->inpoint[i]); }
    for (int i = 0; i < outpoint.size(); i++) { outpoint[i].copyfrom(unitopcopyfrom->outpoint[i]); }
    nin = unitopcopyfrom->nin;
    nout = unitopcopyfrom->nout;
}

void unitop::initinflow()
{
     inflow.resize(nin);
     for (int i = 0; i < nin; i++)
     {
          inflow[i] = new stream(0, location.x, location.y, location.x, location.y);
     }
}

void unitop::initoutflow()
{
     outflow.resize(nout);
     for (int i = 0; i < nout; i++)
     {
          outflow[i] = new stream(0, location.x, location.y, location.x, location.y);
     }
}

void unitop::initinpoint()
{
     inpoint.assign(nin, point(location.x, location.y));  //To be changed by derived classes.
}

void unitop::initoutpoint()
{
     outpoint.assign(nout, point(location.x, location.y));   //To be changed by derived classes.
}

void unitop::updateinoutpointlocations()
{
     for (int i = 0; i < nin; i++)
     {
          if (inflow[i] != nullptr) { inflow[i]->points[1].copyfrom(inpoint[i]); }
     }
     for (int i = 0; i < nout; i++)
     {
          if (outflow[i] != nullptr) { outflow[i]->points[0].copyfrom(outpoint[i]); }
     }
}

void unitop::update(int i, bool historise)
{
     baseprocessclass::update(i, historise);
}

void unitop::draw(QImage *graphics)
{
    QPainter painter(graphics);
    QBrush fillbrush;
    fillbrush.setColor(Qt::red);
    fillbrush.setStyle(Qt::SolidPattern);

    for (int i = 0; i < nin; i++)
    {
        QPolygon poly;
        poly << QPoint(global::OriginX + round(global::GScale*(inpoint[i].x)),
                            global::OriginY + round(global::GScale*(inpoint[i].y)))
             << QPoint(global::OriginX + round(global::GScale*(inpoint[i].x + global::InOutPointWidth)),
                            global::OriginY + round(global::GScale*(inpoint[i].y - global::InOutPointHeight)))
             << QPoint(global::OriginX + round(global::GScale*(inpoint[i].x + global::InOutPointWidth)),
                            global::OriginY + round(global::GScale*(inpoint[i].y + global::InOutPointHeight)));

        painter.drawPolygon(poly);
    }
}










