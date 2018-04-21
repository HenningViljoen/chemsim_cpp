#include "baseprocessclass.h"

baseprocessclass::baseprocessclass(int anr, double ax, double ay)
    : baseclass {anr, ax, ay},
      mat(global::fluidpackage, global::baseprocessclassInitTemperature, global::baseprocessclassInitVolume,
                      global::baseprocessclassInitPressure, 0),
      hasmaterial {true}
{
    actualvolumeflow = new controlvar();
    standardvolumeflow = new controlvar();
    molarflow = new controlvar();
    massflow = new controlvar();

    //mat = new material(global.baseprocessclassInitVolume);
    //public material(string componentname, double aTemp, double aV, double aP, double af) //second constructor

    massflow->v = global::baseprocessclassInitMassFlow;

    calcactualvolumeflowfrommassflow();
    calcmolarflowfrommassflow();
    calcstandardflowfrommoleflow();

    //pressuresimvector = new double[global.SimIterations];

    controlpropthisclass.clear();
    controlpropthisclass = {"pressure",
                                                                            "volume",
                                                                            "density",
                                                                            "temperature",
                                                                            "mass",
                                                                            "n",
                                                                            "actualvolumeflow",
                                                                            "standardvolumeflow",
                                                                            "massflow",
                                                                            "molarflow"};
    for (int i = 0; i < controlpropthisclass.size(); i++)
    {
        controlproperties.push_back(controlpropthisclass[i]);
    }
}

baseprocessclass::~baseprocessclass()
{
    if (actualvolumeflow != nullptr) {delete actualvolumeflow;} //m3/s  non-standard conditions.
    if (standardvolumeflow != nullptr) {delete standardvolumeflow;} //m3/s  standard conditions (later to be descriminated between gases and liquids)
    if (massflow != nullptr) {delete massflow;} //kg/second
    if (molarflow != nullptr) {delete molarflow;} //molar flow per second
}

void baseprocessclass::copyfrom(baseclass *baseclasscopyfrom)
{
     baseprocessclass *baseprocessclasscopyfrom = (baseprocessclass *)baseclasscopyfrom;

     baseclass::copyfrom(baseclasscopyfrom);

     hasmaterial = baseprocessclasscopyfrom->hasmaterial; //True if the stream/unitop/baseclass has material inside that can flow/pump or be processed.

     mat.copyfrom(baseprocessclasscopyfrom->mat);

     actualvolumeflow->v = baseprocessclasscopyfrom->actualvolumeflow->v; //m3/s  non-standard conditions.
     standardvolumeflow->v = baseprocessclasscopyfrom->standardvolumeflow->v; //m3/s  standard conditions (later to be descriminated between gases and liquids)
     massflow->copyfrom(baseprocessclasscopyfrom->massflow); //kg/second  //At this stage we use copyfrom for this one as we need to copy the excel
                                                                  //source in particular as well for the mass flow.
     molarflow->v = baseprocessclasscopyfrom->molarflow->v; //molar flow per second
}

controlvar *baseprocessclass::selectedproperty(int selection)
{
     switch (selection)
     {
          case 0:
               return mat.P;
          case 1:
               return mat.V;
          case 2:
               return mat.density;
          case 3:
               return mat.T;
          case 4:
                return mat.mass;
          case 5:
               return mat.n;
          case 6:
               return actualvolumeflow;
          case 7:
               return standardvolumeflow;
          case 8:
               return massflow;
          case 9:
               return molarflow;
          default:
               return nullptr;

     }
}

        //public void calcmassflowfromactualvolflow()
        //{
        //    massflow = actualvolumeflow * density;
        //}

void baseprocessclass::calcactualvolumeflowfrommassflow()
{
     actualvolumeflow->v = massflow->v / (mat.density->v + 0.001);
}

void baseprocessclass::calcstandardflowfrommoleflow()
{
     standardvolumeflow->v = dndt2fps(molarflow->v);
}

void baseprocessclass::calcmolarflowfrommassflow()
{
     molarflow->v = 0;
     molarflow->v = massflow->v / (mat.massofonemole + 0.001);
}

        //public void calcmassflowfromstandardflow()
        //{
        //    molarflow = utilities.fps2dndt(standardvolumeflow);
        //    massflow = molarflow * mat.massofonemole;
        //}

void baseprocessclass::calcmassflowfrommolarflow()
{
     massflow->v = 0;
     for (int i = 0; i < mat.composition.size(); i++)
     {
          massflow->v += mat.composition[i].n / mat.n->v * molarflow->v *
                    mat.composition[i].m->molarmass;
     }
}

void baseprocessclass::update(int i, bool historise) //index for where in the simvectors the update is to be stored.
{
     mat.update(i, historise);
}

        //public virtual void updatepoint(int i, double x, double y)
        //{
        //}

        //public virtual bool mouseover(double x, double y) //This function will indicate whether the mouse is over a particular unitop or stream at any moment in time.
        //{
        //    return false;
        //}

        //public virtual void setproperties(simulation asim) //Method that will be inherited and that will set the properties of the applicable object in a window
        //{

        //}

        //public virtual void draw(Graphics G)
        //{
        //}


