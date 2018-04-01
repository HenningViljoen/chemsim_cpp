#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
//#include <QWidget>
#include <QObject>
#include "unitop.h"

class simulation
{
public:
    std::vector<unitop *> unitops;

    simulation();
    ~simulation();

    void drawnetwork(QImage *graphics);
};

#endif // SIMULATION_H
