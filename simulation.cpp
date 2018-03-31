#include "simulation.h"
#include "valve.h"

simulation::simulation()
{

}

simulation::~simulation()
{
    for (int i = 0; i < unitops.size(); i++) {
        delete (unitops[i]);
    }
    unitops.clear();
}

void simulation::drawnetwork(QImage G)
{

}
