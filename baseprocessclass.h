#ifndef BASEPROCESSCLASS_H
#define BASEPROCESSCLASS_H

#include "baseclass.h"
#include "material.h"

class baseprocessclass : public baseclass
{
public:
    bool hasmaterial; //True if the stream/unitop/baseclass has material inside that can flow/pump or be processed.

    material mat;

    controlvar *actualvolumeflow = nullptr; //m3/s  non-standard conditions.
    controlvar *standardvolumeflow = nullptr; //m3/s  standard conditions (later to be descriminated between gases and liquids)
    controlvar *massflow = nullptr; //kg/second
    controlvar *molarflow = nullptr; //molar flow per second

    baseprocessclass(int anr, double ax, double ay);
    ~baseprocessclass();

    void copyfrom(baseclass *baseclasscopyfrom);
    controlvar *selectedproperty(int selection);
    void calcactualvolumeflowfrommassflow();
    void calcstandardflowfrommoleflow();
    void calcmolarflowfrommassflow();
    void calcmassflowfrommolarflow();
    void update(int i, bool historise);
};

#endif // BASEPROCESSCLASS_H
