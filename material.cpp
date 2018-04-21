#include "material.h"

material::material(double aV)
{
    P = new controlvar();
    V = new controlvar();
    //vmolar = new controlvar();
    //vmolarL = new controlvar();
    //vmolarV = new controlvar();
    density = new controlvar();
    T = new controlvar();
    mass = new controlvar();
    n = new controlvar();
    //f = new controlvar();
    //U = new controlvar();
    //umolar = new controlvar();
    //relativehumidity = new controlvar();

    double oldn = 0.0;
    for (int i = 0; i < global::fluidpackage.size(); i++)
    {
        if (global::fluidpackage[i].n != 0)
        {
             composition.push_back(global::fluidpackage[i]);
             oldn += composition[composition.size() - 1].n;
        }
    }

    x.assign(composition.size(),0.0);
    y.assign(composition.size(),0.0);
    z.assign(composition.size(),0.0);
    fugacityL.assign(composition.size(),0.0);
    fugacityV.assign(composition.size(),0.0);
    K.resize(composition.size());
    Z.resize(composition.size());
    for (int i = 0; i < composition.size(); i++)
    {
        Z[i].assign(6,0.0);
        nullZ(i);
    }

    discriminantZ.assign(composition.size(),0.0);
    Tr.assign(composition.size(),0.0);
    ac.assign(composition.size(),0.0);
    aT.assign(composition.size(),0.0);
    alpha.assign(composition.size(),0.0);
    b.assign(composition.size(),0.0);
    A.assign(composition.size(),0.0);
    B.assign(composition.size(),0.0);
    delta0.assign(composition.size(),0.0);
    delta1.assign(composition.size(),0.0);
    C.assign(composition.size(),0.0);
    Cc.resize(composition.size());

    Cp.assign(composition.size(),0.0);

    umolarL.assign(composition.size(),0.0);
    umolarV.assign(composition.size(),0.0);
    umolarideal.assign(composition.size(),0.0);

    acomp.assign(composition.size(),0.0); //Cubic equation for compressibility factor coeficients.
    bcomp.assign(composition.size(),0.0); //Cubic equation for compressibility factor coeficients.
    ccomp.assign(composition.size(),0.0); //Cubic equation for compressibility factor coeficients.
    dcomp.assign(composition.size(),0.0); //Cubic equation for compressibility factor coeficients.

    a.resize(composition.size());
    for (int i = 0; i < composition.size(); i++) { a[i].assign(3, 0.0); }

    Qc.assign(composition.size(),0.0); //used in alternative formula.
    Rc.assign(composition.size(),0.0); //used in alternative formula.
    Dc.assign(composition.size(),0.0); //used in alternative formula.
    Sc.resize(composition.size()); //used in alternative formula.
    Tc.resize(composition.size()); //used in alternative formula.

    zs.resize(6); //Used in scarlet cubic solver.
    us.resize(2); //Used in scarlet cubic solver.
    ys.resize(6); //Used in scarlet cubic solver.

    uvflashsize = 2 * composition.size() + 3;
    xvector.assign(uvflashsize, 0.0);
    for (int i = 0; i < composition.size(); i++)
    {
        x[i] = 1.0 / (composition.size()); //Add 0.5 to get starting values that can converge.
        y[i] = 1.0 / (composition.size());
    }
    if (composition.size() == 1) { uvflashsize = 3; } //x and y does not need to be part of the model anymore.
    origuvflashsize = uvflashsize;

    jacobian.initmatrix(uvflashsize, uvflashsize);

    V->v = aV;
    P->v = global::baseprocessclassInitPressure;
    T->v = global::baseprocessclassInitTemperature;
    f.v = global::fdefault;

    mapvarstox(); //Variables are allocated according to the Flatby article sequence of equations.
    n->v = V->v / calcvmolaranddensity();

    for (int i = 0; i < composition.size(); i++)
    {
         composition[i].n = composition[i].n / oldn * n->v;
    }

    calcmass();

    U.v = n->v * calcumolar();

    calccompz();
    uvflash(); //  This is commented out now, but needs to be put back in soon again.

    density->simvector.resize(global::SimVectorLength);
}

material::material(std::vector<component> &acomposition, double aTemp, double aV, double aP, double af) //second constructor
{
    P = new controlvar(); //Pa
    V = new controlvar(); //m3.  This will be total stage volume in the case of a distillation column.
    //vmolar = new controlvar(); //m3/mole.  Total volume per mole for the total material all phases.
    //vmolarL = new controlvar(); //m3/mole.  Total volume per mole for the total liquid phase.
    //vmolarV = new controlvar(); //m3/mole.  Total volume per mole for the total vapour phase.
    density = new controlvar(); //kg/m3
    T = new controlvar(); //Kelvin
    mass = new controlvar(); //kg
    n = new controlvar(); //amount of moles of the material (all phases) for all components.
    //f = new controlvar(); //fraction.  Vapour molar fraction.
    //U = new controlvar(); //Joule. Internal energy.
    //umolar = new controlvar();
    //relativehumidity = new controlvar();
    //composition = new List<component>();

    init(acomposition, aTemp, aV, aP, af);

            //P.simvector = new double[global.SimIterations];
            //T.simvector = new double[global.SimIterations];
            //f.simvector = new double[global.SimIterations];
            //n.simvector = new double[global.SimIterations];
            //U.simvector = new double[global.SimIterations];
    density->simvector.resize(global::SimVectorLength);
}

material::material(const material &m) //copy constructor
{
    copyfrom(m);
}

material::~material()
{
    if (P != nullptr) {delete P;}

    if (V != nullptr) {delete V;}

    if (density != nullptr) {delete density;}

    if (T != nullptr) {delete T;}

    if (mass != nullptr) {delete mass;}
    if (n != nullptr) {delete n;}
}

material &material::operator=(const material &m)
{
    material mat(m);
    return mat;
}

void material::copyfrom(const material &materialcopyfrom)
{
    P->v = materialcopyfrom.P->v; //Pa
    V->v = materialcopyfrom.V->v; //m3.  This will be total stage volume in the case of a distillation column.
    vmolar.v = materialcopyfrom.vmolar.v; //m3/mole.  Total volume per mole for the total material all phases.
    vmolarL.v = materialcopyfrom.vmolarL.v; //m3/mole.  Total volume per mole for the total liquid phase.
    vmolarV.v = materialcopyfrom.vmolarV.v; //m3/mole.  Total volume per mole for the total vapour phase.
    density->v = materialcopyfrom.density->v; //kg/m3
    T->copyfrom(materialcopyfrom.T); //Kelvin
    mass->v = materialcopyfrom.mass->v; //kg
    n->v = materialcopyfrom.n->v; //amount of moles of the material (all phases) for all components.
    f.v = materialcopyfrom.f.v; //fraction.  Vapour molar fraction.
    U.v = materialcopyfrom.U.v; //Joule. Internal energy.
    relativehumidity.v = materialcopyfrom.relativehumidity.v;

    x = materialcopyfrom.x; //Molar fraction of the total liquid of the material per component.
    y = materialcopyfrom.y; //Molar fraction of the total vapour of the material per component.
    hL = materialcopyfrom.hL; //Joule/mol.  Molar enthalpy for the liquid phase, e.g. on a stage in a distillation column.
    hV = materialcopyfrom.hV; //Joule/mol.  Molar enthalpy for the vapour phase, e.g. on a stage in a distillation column.
    ML = materialcopyfrom.ML; //moles.  Total molar liquid hold-up.
    MV = materialcopyfrom.MV; //moles.  Total molar vapour hold-up.
    z = materialcopyfrom.z; //fraction.  Molar fraction per component of the total moles in the material.

    umolar.v = materialcopyfrom.umolar.v; //Joule/mol.  Molar internal energy.
    umolarL = materialcopyfrom.umolarL; //Joule/mol.  Molar internal energy for the total liquid phase.
    umolarV = materialcopyfrom.umolarV; //Joule/mol.  Molar internal energy for the total vapour phase.
    double totalumolar = materialcopyfrom.totalumolar; //last calculation of the total vmolar;
    umolarideal = materialcopyfrom.umolarideal; //Joule/mol.  Molar internal energy for the ideal gas case.
    fugacityL = materialcopyfrom.fugacityL; //fugacity per component for the liquid phase.
    fugacityV = materialcopyfrom.fugacityV; //fugacity per component for the vapour phase.
    aT = materialcopyfrom.aT; //kg⋅m^5 /( s^2 ⋅mol^2 ); figure used in Peng Robinson equation.
    ac = materialcopyfrom.ac; // kg⋅m^5 / (s^2 ⋅mol^2 ) ; figure used in Peng Robinson equation.
    alpha = materialcopyfrom.alpha; //part of the equation for aT, and used in other equations in the Peng Robinson collection.
    b = materialcopyfrom.b;  // m^3/mol ; figure used in Peng Robinson equation.
    K = materialcopyfrom.K;  //dimensionless - depends on omega.
    A = materialcopyfrom.A;  //used in equation for compressibility factor.
    B = materialcopyfrom.B;  //used in equation for compressibility factor.
    Z = materialcopyfrom.Z;  //Compressibility factor. More than one value should be present based on the amount of phases in the material
    discriminantZ = materialcopyfrom.discriminantZ; //The discriminant of the cubig PR equation for Z;
    Tr = materialcopyfrom.Tr; //reduced temperature.
    acomp = materialcopyfrom.acomp; //Cubic equation for compressibility factor coeficients.
    bcomp = materialcopyfrom.bcomp; //Cubic equation for compressibility factor coeficients.
    ccomp = materialcopyfrom.ccomp; //Cubic equation for compressibility factor coeficients.
    dcomp = materialcopyfrom.dcomp; //Cubic equation for compressibility factor coeficients.
    bs = materialcopyfrom.bs; //Used in scarlet cubic solver. bs for bscarlet, This is for another cubic solver: http://home.scarlet.be/~ping1339/cubic.htm
    cs = materialcopyfrom.cs; //Used in scarlet cubic solver.
    ds = materialcopyfrom.ds; //Used in scarlet cubic solver.
    rs = materialcopyfrom.rs; // Used in scarlet cubic solver.
    es = materialcopyfrom.es; //Used in scarlet cubic solver.
    fs = materialcopyfrom.fs; //Used in scarlet cubic solver.
            //public complex[] zs; //Used in scarlet cubic solver. I DO NOT THINK WE NEED TO COPY THESE VALUES AS THEY ARE INITED ON EACH ITERATION.
            //public complex[] us; //Used in scarlet cubic solver.
            //public complex[] ys; //Used in scarlet cubic solver.


    delta0 = materialcopyfrom.delta0; //used to solve the cubic equation for Z.
    delta1 = materialcopyfrom.delta1; //used to solve the cubic equation for Z.
    C = materialcopyfrom.C; //used to solve the cubic equation for Z.
    Cc = materialcopyfrom.Cc; //used to solve the cubic equation for Z.
    a = materialcopyfrom.a; //used to solve the cubid equation for Z by making use of an alternative closed form formulat for them.
    Qc = materialcopyfrom.Qc; //used in alternative formula.
    Rc = materialcopyfrom.Rc; //used in alternative formula.
    Dc = materialcopyfrom.Dc; //used in alternative formula.
    Sc = materialcopyfrom.Sc; //used in alternative formula.
    Tc = materialcopyfrom.Tc; //used in alternative formula.
    Cp = materialcopyfrom.Cp; //Heat capacity of constant pressure for an ideal gas type scenario per component.

    totalCp = materialcopyfrom.totalCp; //Heat capacity of the material in totallity, with all its components includded.
    xvector = materialcopyfrom.xvector;
    //public matrix jacobian; //The jacobian matrix for the system of equations as defined below in fflash().
    uvflashsize = materialcopyfrom.uvflashsize; //The size of the k variables k equations problem to be solved for the UV flash.
    origuvflashsize = materialcopyfrom.origuvflashsize;

    //copycompositiontothismat(materialcopyfrom); THIS HAS TO BE UNCOMMENTED WHEN THE CODE PORT IS DONE.
    massofonemole = materialcopyfrom.massofonemole; //kg/mole
    volumeofonemole = materialcopyfrom.volumeofonemole; //m^3/mole
    phase = materialcopyfrom.phase;

            //componentindex = materialcopyfrom.componentindex;
}



void material::copycompositiontothismat(material &amat)
{
    composition.clear();
    for (int i = 0; i < amat.composition.size(); i++)
    {
        composition.push_back(component());

        composition[composition.size() - 1].copytothisobject(amat.composition[i]);
        //componentindex = i; //this will fall away later when there is more than one component per stream.
        composition[composition.size() - 1].molefraction = amat.composition[i].molefraction;
        //composition[composition.Count - 1].n = amat.composition[i].n / amat.n.v * n.v;  //this one will need to be tested might not be right.
     }
}

void material::copycompositiontothiscomposition(std::vector<component> &compositioncopyto, std::vector<component> &compositioncopyfrom)
{
     compositioncopyto.clear();
     for (int i = 0; i < compositioncopyfrom.size(); i++)
     {
         compositioncopyto.push_back(component());

         compositioncopyto[compositioncopyto.size() - 1].copytothisobject(compositioncopyfrom[i]);

         compositioncopyto[compositioncopyto.size() - 1].molefraction = compositioncopyfrom[i].molefraction;
     }
}

void material::init(std::vector<component> &acomposition, double aTemp, double aV, double aP, double af)
{
     //phase = global.MaterialInitPhase;

     composition.clear();
     for (int i = 0; i < acomposition.size(); i++)
     {
         //if (global.fluidpackage[i].m.name == componentname)
         //{
         composition.push_back(component());

         composition[composition.size() - 1].copytothisobject(acomposition[i]);
         //componentindex = i; //Should not be needed anymore soon.
         //}
     }
     x.assign(composition.size(), 0.0);
     y.assign(composition.size(), 0.0);
     z.assign(composition.size(), 0.0);
     fugacityL.assign(composition.size(), 0.0);
     fugacityV.assign(composition.size(), 0.0);
     K.assign(composition.size(), 0.0);

     Z.resize(composition.size());
     for (int i = 0; i < composition.size(); i++)
     {
          Z[i].assign(6, 0.0);
          nullZ(i);
     }

     discriminantZ.assign(composition.size(), 0.0);
     Tr.assign(composition.size(), 0.0);
     ac.assign(composition.size(), 0.0);
     aT.assign(composition.size(), 0.0);
     alpha.assign(composition.size(), 0.0);
     b.assign(composition.size(), 0.0);
     A.assign(composition.size(), 0.0);
     B.assign(composition.size(), 0.0);
     delta0.assign(composition.size(), 0.0);
     delta1.assign(composition.size(), 0.0);
     C.assign(composition.size(), 0.0);
     Cc.assign(composition.size(), complex());
            //for (int i = 0; i < composition.Count; i++)
            //{
            //    Cc[i] = new complex(0, 0);
            //}

     Cp.assign(composition.size(), 0.0);
     totalCp = 0;

     umolarL.assign(composition.size(), 0.0);
     umolarV.assign(composition.size(), 0.0);
     umolarideal.assign(composition.size(), 0.0);

     acomp.assign(composition.size(), 0.0); //Cubic equation for compressibility factor coeficients.
     bcomp.assign(composition.size(), 0.0); //Cubic equation for compressibility factor coeficients.
     ccomp.assign(composition.size(), 0.0); //Cubic equation for compressibility factor coeficients.
     dcomp.assign(composition.size(), 0.0); //Cubic equation for compressibility factor coeficients.

     a.resize(composition.size());
     for (int i = 0; i < composition.size(); i++) { a[i].assign(3, 0.0); }

     Qc.assign(composition.size(), 0.0); //used in alternative formula.
     Rc.assign(composition.size(), 0.0); //used in alternative formula.
     Dc.assign(composition.size(), 0.0); //used in alternative formula.
     Sc.assign(composition.size(), complex()); //used in alternative formula.
     Tc.assign(composition.size(), complex()); //used in alternative formula.

     zs.assign(6, complex()); //Used in scarlet cubic solver.
     us.assign(2, complex()); //Used in scarlet cubic solver.
     ys.assign(6, complex()); //Used in scarlet cubic solver.

     uvflashsize = 2 * composition.size() + 3;
     xvector.assign(uvflashsize, 0.0);
     for (int i = 0; i < composition.size(); i++)
     {
          x[i] = 1.0 / (composition.size()); //Add 0.5 to get starting values that can converge.
          y[i] = 1.0 / (composition.size());
     }
     if (composition.size() == 1) { uvflashsize = 3; } //x and y does not need to be part of the model anymore.
     origuvflashsize = uvflashsize;

     jacobian.initmatrix(uvflashsize, uvflashsize);

     PTfVflash(aTemp, aV, aP, af);
     //calccompz(); Already done in the previous flash ptfv

}

void material::setxybasedonf() //If f = 0 or f = 1, then the values of x and y can be preset based on z and full flashing is not needed
{
     for (int i = 0; i < composition.size(); i++)
     {
          if (f.v == 0)
          {
               x[i] = z[i];
               y[i] = 0;
          }
          else if (f.v == 1)
          {
               x[i] = 0;
               y[i] = z[i];
          }
     }
}

void material::PTfVflash(double aTemp, double aV, double aP, double af)
{
     //composition.Clear();
     //for (int i = 0; i < acomposition.Count; i++)
     //{
     //    //if (global.fluidpackage[i].m.name == componentname)
     //    //{
     //    composition.Add(new component());

     //    composition[composition.Count - 1].copytothisobject(acomposition[i]);
     //    //componentindex = i;
     //    //}
     //}

     P->v = aP;
     T->v = aTemp;
     f.v = af;
     V->v = aV;
     calccompz();
     setxybasedonf();
     mapvarstox(); // Variables are allocated according to the Flatby article sequence of equations.
     //U.v = 10000000.0;
     //uvflash();
     vmolar.v = calcvmolaranddensity();
     n->v = (V->v / vmolar.v);
     //composition[0].n = n.v;  //ASSUMING ONLY ONE COMPONENT AT THE MOMENT.  Should not be needed anymore.

     //calcn(); //Should not be needed anymore.
     calcmass();
     //massofonemole = mass.v / n.v; Already caculated in calcvmolarandensity;

     umolar.v = calcumolar();
     U.v = n->v * umolar.v;
     //calctotalCp(); THIS HAS TO BE UNCOMMENTED WHEN THE CODE PORT IS DONE.

}

void material::nullZ(int j)
{
     for (int i = 0; i < 3; i++) { Z[j][i] = global::ZNotDefined; }
}

//public void calcn() //This method could not be needed anymore.
        //{
        //    n.v = 0;
        //    for (int i = 0; i < composition.Count; i++)
        //    {
        //        n.v += composition[i].n;
        //    }
        //}

        //public void calcmassanddensity()
        //{
        //    calcn();
        //    massofonemole = 0;
        //    volumeofonemole = 0;
        //    double masscontribution, volumecontribution, moleof1kg = 0;
        //    for (int i = 0; i < composition.Count; i++)
        //    {
        //        moleof1kg += composition[i].massfraction / composition[i].m.molarmass;
        //    }
        //    for (int i = 0; i < composition.Count; i++)
        //    {
        //        composition[i].molefraction = composition[i].massfraction / composition[i].m.molarmass / moleof1kg;
        //        masscontribution = composition[i].molefraction * composition[i].m.molarmass;
        //        volumecontribution = masscontribution / (composition[i].m.density + 0.00001); //Adding a small number to negate
        //        //very singularities.
        //        massofonemole += masscontribution;
        //        volumeofonemole += volumecontribution;
        //    }
        //    density.v = massofonemole / volumeofonemole;
        //}

void material::calcmass()
{
     mass->v = 0.0;
     for (int i = 0; i < composition.size(); i++)
     {
          mass->v += n->v*composition[i].molefraction * composition[i].m->molarmass;
     }
}

component &material::match(std::string aname)
{
     component c;
     for (int i = 0; i < composition.size(); i++)
     {
           if (composition[i].m->name == aname) { c = composition[i]; }
     }
     return c;
}

//public double totalmolefraction()
        //{
        //    double total = 0;
        //    for (int i = 0; i < composition.Count; i++) { total += composition[i].molefraction; } //This should now be adding up to 100%, if it
        //    //does not add up, then a scaling will need to
        //    //be done in order to have the total to be 100%
        //    return total;
        //}

void material::copycompositiontothisobject(material &m)
{
    for (int i = 0; i < composition.size(); i++)
    {
         composition[i].copytothisobject(m.composition[i]);
    }
}

void material::mapvarstox() //Variables are allocated according to the Flatby article sequence of equations.
{
    xvector[0] = T->v;
    xvector[1] = P->v;
    xvector[2] = f.v;
    for (int i = 0; i < composition.size(); i++)
    {
         xvector[3 + i] = x[i];
         xvector[3 + composition.size() + i] = y[i];
    }
}

void material::mapxtovars() //Variables are allocated according to the Flatby article sequence of equations.
{
    T->v = xvector[0];
    P->v = xvector[1];
    f.v = xvector[2];

    for (int i = 0; i < composition.size(); i++)
    {
         x[i] = xvector[3 + i];
         y[i] = xvector[3 + composition.size() + i];
    }
}

complex &material::calcxk(complex &u, int j)
{
    return -1.0 / (3.0 * acomp[j]) * (bcomp[j] + u * Cc[j] + delta0[j] / (u * Cc[j]));
}

bool complexcomparer(complex &x, complex &y)
{
    if (abs(x.b) < global::ZeroImaginary && abs(y.b) >= global::ZeroImaginary)
    {
        return false;
    }
    else if (abs(x.b) >= global::ZeroImaginary && abs(y.b) < global::ZeroImaginary)
    {
         return true;
    }
    else if (abs(x.b) < global::ZeroImaginary && abs(y.b) < global::ZeroImaginary)
    {
         if (x.a < y.a) { return false; }
         else if (x.a > y.a) { return true; }
         else { return true; }
    }
    else
    {
         return true;
    }
}

void material::calcZ(int j) //Calculate the roots of the compressibility equation for the particular component.
{
    if (composition[j].m->omega <= 0.49)
    {
        K[j] = 0.37464 + 1.54226 * composition[j].m->omega - 0.26992 * pow(composition[j].m->omega, 2.0);
    }
    else
    {
        K[j] = 0.379642 + 1.48503 * composition[j].m->omega - 0.164423 * pow(composition[j].m->omega, 2.0) +
                0.016666 * pow(composition[j].m->omega, 3.0);
    }
    Tr[j] = xvector[0] / composition[j].m->Tc; //T / composition[j].m.Tc;
    ac[j] = 0.45724 * pow(global::R, 2.0) * pow(composition[j].m->Tc, 2.0) / composition[j].m->Pc;
    alpha[j] = pow(1 + K[j] * (1 - sqrt(Tr[j])), 2.0);
    aT[j] = alpha[j] * ac[j];
    b[j] = 0.07780 * global::R * composition[j].m->Tc / composition[j].m->Pc;
    A[j] = aT[j] * xvector[1] / (pow(global::R, 2.0) * pow(xvector[0], 2.0)); //xvector[1]: P
    B[j] = b[j] * xvector[1] / (global::R * xvector[0]);

    acomp[j] = 1.0;
    bcomp[j] = (B[j] - 1.0);
    ccomp[j] = A[j] - 2.0 * B[j] - 3.0 * pow(B[j], 2.0);
    dcomp[j] = pow(B[j], 3.0) + pow(B[j], 2.0) - A[j] * B[j];

    //Delta = 18abcd -4b^3d + b^2c^2 - 4ac^3 - 27a^2d^2
    discriminantZ[j] = 18.0 * acomp[j] * bcomp[j] * ccomp[j] * dcomp[j] - 4.0 * pow(bcomp[j], 3.0) * dcomp[j] + pow(bcomp[j], 2.0)
            * pow(ccomp[j], 2.0) -
                4.0 * acomp[j] * pow(ccomp[j], 3.0) - 27.0 * pow(acomp[j], 2.0) * pow(dcomp[j], 2.0);

    //Delta_0 = b^2-3 a c
    //delta0[j] = Math.Pow(bcomp[j], 2.0) - 3.0 * acomp[j] * ccomp[j];

    //Delta_1 = 2 b^3-9 a b c+27 a^2 d
    //delta1[j] = 2.0 * Math.Pow(bcomp[j], 3.0) - 9.0 * acomp[j] * bcomp[j] * ccomp[j] + 27.0 * Math.Pow(acomp[j], 2.0) * dcomp[j];

    //C[j] = Math.Pow((delta1[j] + Math.Sqrt(Math.Pow(delta1[j], 2.0) - 4.0 * Math.Pow(delta0[j], 3.0))) / 2.0, 1.0 / 3.0);

    //Cc[j].a = Cc[j].b = 0.0;
    //complex Coper = new complex(Math.Pow(delta1[j], 2.0) - 4.0 * Math.Pow(delta0[j], 3.0), 0.0);
    //if (Coper < 0) { Cc[j].b = Math.Sqrt(-Coper); }
    //else { Cc[j].a = Math.Sqrt(Coper); }
    //Cc[j] = complex.Pow((delta1[j] - complex.Pow(Coper,0.5))/2.0, 1.0/3.0);

    nullZ(j);
    std::vector<complex> xk(3);
    //complex[] u = new complex[] {new complex(1,00), new complex(-0.5,Math.Pow(3,0.5)/2.0), new complex(-0.5, -Math.Pow(3,0.5)/2.0)};

    //Now an alternative way of getting the roots will be tried.  Wolfram algorithm.
    a[j][0] = dcomp[j] / acomp[j];
    a[j][1] = ccomp[j] / acomp[j];
    a[j][2] = bcomp[j] / acomp[j];
    //Qc[j] = (3.0 * a[j][1] - Math.Pow(a[j][2], 2.0)) / 9.0;
    //Rc[j] = (9.0 * a[j][2] * a[j][1] - 27.0 * a[j][0] - 2.0 * Math.Pow(a[j][2], 3.0)) / 54.0;
    //Dc[j] = Math.Pow(Qc[j], 3.0) + Math.Pow(Rc[j], 3.0);
    //Sc[j] = complex.pow(Rc[j] + complex.pow(Dc[j], 0.5), 1.0/3.0);
    //Tc[j] = complex.pow(Rc[j] - complex.pow(Dc[j], 0.5), 1.0/3.0);
    //xk[0] = -1.0 / 3.0 * a[j][2] + Sc[j] + Tc[j];
    //xk[1] = -1.0 / 3.0 * a[j][2] - 1.0 / 2.0 * (Sc[j] + Tc[j]) + 1.0 / 2.0 * complex.I * Math.Pow(3.0, 0.5) * (Sc[j] - Tc[j]);
    //xk[2] = -1.0 / 3.0 * a[j][2] - 1.0 / 2.0 * (Sc[j] + Tc[j]) - 1.0 / 2.0 * complex.I * Math.Pow(3.0, 0.5) * (Sc[j] - Tc[j]);

    //This is still another algorithm for the solution.  The scarlet algorithm.  http://home.scarlet.be/~ping1339/cubic.htm

    std::vector<complex> xs; //Used in scarlet cubic solver.

    bs = a[j][2];
    cs = a[j][1];
    ds = a[j][0];
    rs = -bs / 3.0;

    //bs = 8.0/15.0;
    //cs = -7.0/45.0;
    //ds = -2.0/45.0;
    //rs = -bs / 3.0;

    //y^3  + (3 r + b) y^2  + (3 r^2  + 2 r b + c) y + r^3  + r^2  b + r c + d = 0
    //y^3  + e y + f = 0
    es = 3.0 * pow(rs, 2.0) + 2 * rs * bs + cs;
    fs = pow(rs, 3.0) + pow(rs, 2.0) * bs + rs * cs + ds;

    double bqz = fs; //The b term of the quadratic equation in z
    double cqz = -pow(es, 3.0) / 27.0;

    us[0] = (-bqz + complex::power(pow(bqz, 2.0) - 4 * cqz, 0.5)) / (2.0);
    us[1] = (-bqz - complex::power(pow(bqz, 2.0) - 4 * cqz, 0.5)) / (2.0);

    zs[0] = complex::power(us[0], 1.0 / 3.0, 0);
    zs[1] = complex::power(us[0], 1.0 / 3.0, 1);
    zs[2] = complex::power(us[0], 1.0 / 3.0, 2);
    zs[3] = complex::power(us[1], 1.0 / 3.0, 0);
    zs[4] = complex::power(us[1], 1.0 / 3.0, 1);
    zs[5] = complex::power(us[1], 1.0 / 3.0, 2);

    double temp = -es / 3.0;

    ys[0] = zs[0] - es / (3.0 * zs[0]);
    ys[1] = zs[1] - es / (3.0 * zs[1]);
    ys[2] = zs[2] - es / (3.0 * zs[2]);
    ys[3] = zs[3] - es / (3.0 * zs[3]);
    ys[4] = zs[4] - es / (3.0 * zs[4]);
    ys[5] = zs[5] - es / (3.0 * zs[5]);

    xs.push_back(ys[0] - bs / 3.0);
    xs.push_back(ys[1] - bs / 3.0);
    xs.push_back(ys[2] - bs / 3.0);
    xs.push_back(ys[3] - bs / 3.0);
    xs.push_back(ys[4] - bs / 3.0);
    xs.push_back(ys[5] - bs / 3.0);

    //complexcomparer complexcompare = new complexcomparer();
    std::sort(xs.begin(), xs.end(), complexcomparer);

    for (int k = 1; k < xs.size(); k++)
    {
        if (round_precision(xs[k].a, 6) == round_precision(xs[k - 1].a, 6)) { xs.erase(xs.begin() + k); }
    }

    //int nbelow05 = 0; //nr of roots that are below 0.5.
    //for (int k = 0; k < xs.Count; k++)
    //{
    //    if (xs[k].a < 0.5 && Math.Abs(xs[k].b) < global.ZeroImaginary) { nbelow05++; }
    //}
    if (xs.size() == 3) { xs.erase(xs.begin() + 1); } //Assume that the higher root value should be taken out.  To be verified.

    //End of Scarlet algorithm.

    int Zsvalid = 0;
    int Zsnotvalid = 0;
    int Zthatsvalid = 0;
    int Zthatsnotvalid = 0;
    int badZi = xs.size() - 1;
    int goodZi = 0;
    int rightiforZ;

    for (int k = 0; k < xs.size(); k++)
    {
        //xk[k] = calcxk(u[k], j);
        if (abs(xs[k].b) < global::ZeroImaginary && xs[k].a >= 0)
        {
            Z[j][goodZi] = xs[k].a;
            Zsvalid++;
            Zthatsvalid = goodZi;
            rightiforZ = goodZi;
            //for (int l = goodZi - 1; l >= 0; l--)
            //{
            //    if (Z[j][l] > Z[j][goodZi]) {rightiforZ = l;}
            //}
            //for (int m = goodZi; m > rightiforZ; m--)
            //{
            //    Z[j][m] = Z[j][m - 1];
            //}
            //Z[j][rightiforZ] = xk[k].a;

            goodZi++;
         }
         else
         {
             Z[j][badZi--] = global::ZNotDefined;
             Zsnotvalid++;
         }
    }
    if (Zsvalid > 1)
    {
         if (f.v == 0 || xvector[2] == 0)
         {
               //f.v = 0.1;
               //xvector[2] = f.v;
         }
         else if (f.v == 1 || xvector[2] == 1)
         {
               //f.v = 0.9;
               //xvector[2] = f.v;
         }
    }
    else if (Zsvalid == 1)
    {
         if (Z[j][Zthatsvalid] < 0.5 && composition.size() == 1)
         {
               f.v = 0;
               xvector[2] = f.v;
         } //only liquid
         else if (Z[j][Zthatsvalid] > 0.5 && composition.size() == 1)
         {
               f.v = 1;
               xvector[2] = f.v;
         } //only vapour
    }

    int dummy = 0;
    for (int k = 0; k < 3; k++)
    {
         if (Z[j][k] == global::ZNotDefined) { Z[j][k] = Z[j][Zthatsvalid]; }
    }



    //if (discriminantZ[j] < 0) //Should be the case if there is only one phase of the molecule present
    //{
    //    Z[j][0] = -1 / (3 * acomp[j]) * (bcomp[j] + C[j] + delta0[j] / C[j]);
    //    Z[j][1] = Z[j][0];
    //}
    //else
    //{
    //    Z[j][0] = (9*acomp[j]*dcomp[j] - bcomp[j]*ccomp[j])/(2*delta0[j]); //This should be the liquid root.
    //    Z[j][1] = (4*acomp[j]*bcomp[j]*ccomp[j] - 9*Math.Pow(acomp[j],2)*dcomp[j] - Math.Pow(bcomp[j],3))/(acomp[j]*delta0[j]); //This should be the gas root.
    //}
}

void material::calcfugacity(int j) //assume that Z has already been calculated for both roots.  Z[j][0] is the liquid root, Z[j][1] is the gas root.
{
    std::vector<double> logtheta = {0.0, 0.0};
    calcZ(j);
    for (int i = 0; i < 2; i++)
    {
        if (Z[j][i] != global::ZNotDefined)
        {
             logtheta[i] = (Z[j][i] - 1.0) - log(Z[j][i] - B[j]) -
                        A[j] / (2.0 * sqrt(2.0) * B[j]) * log((Z[j][i] + (sqrt(2.0) + 1.0) * B[j]) / (Z[j][i] - (sqrt(2.0) - 1.0) * B[j]));
        }
    }
    fugacityL[j] = exp(logtheta[0]) * xvector[1];
    fugacityV[j] = exp(logtheta[1]) * xvector[1];
    //if (Z[j][0] == global.ZNotDefined) {fugacityL[j] = fugacityV[j];}
    //else if (Z[j][1] == global.ZNotDefined) { fugacityV[j] = fugacityL[j]; }
}

double material::calcumolarwithZ(int j, int Zi) //calcs the molar internal energy for the given root of the Z equation, and for component j.
{
     return -A[j] / (B[j] * sqrt(8.0)) * (1.0 + K[j] * sqrt(Tr[j]) / sqrt(alpha[j])) *
          log((Z[j][Zi] + (1.0 + sqrt(2.0)) * B[j]) / (Z[j][Zi] + (1.0 - sqrt(2.0)) * B[j])) * global::R * xvector[0] +    //xvector[0]: T
          umolarideal[j];
}

void material::calcumolarpercomponent()
{
    for (int j = 0; j < composition.size(); j++)
    {
        calcZ(j);
        Cp[j] = composition[j].m->calcCp(xvector[0]);
        umolarideal[j] = (Cp[j] - global::R) * xvector[0];
        //umolarideal[j] = 0.0;
        umolarL[j] = calcumolarwithZ(j, 0);
        umolarV[j] = calcumolarwithZ(j, 1);
     }
}

void material::calctotalCp()
{
    calccompz();
    totalCp = 0;
    for (int i = 0; i < composition.size(); i++)
    {
         totalCp += Cp[i]*z[i];
    }
}

double material::calcumolar()
{
    //calcumolarpercomponent(); THIS HAS TO BE UNCOMMENTED WHEN THE CODE PORT IS DONE.

    double totalumolarL = 0.0;
    double totalumolarV = 0.0;
    for (int k = 0; k < composition.size(); k++)
    {
         totalumolarL += umolarL[k] * xvector[3 + k]*composition[k].molefraction;
         totalumolarV += umolarV[k] * xvector[3 + composition.size() + k] * composition[k].molefraction;
    }
    totalumolar = (1 - xvector[2]) * totalumolarL + xvector[2] * totalumolarV; //f: xvector[2]
    return totalumolar;
}

double material::calcvmolaranddensity()
{
    double totalvmolar = 0.0;
    double totalvmolarL = 0.0;
    double totalvmolarV = 0.0;
    massofonemole = 0.0;
    for (int k = 0; k < composition.size(); k++)
    {
        calcZ(k);
        totalvmolarL += b[k] * Z[k][0] / B[k] * xvector[3 + k] * composition[k].molefraction;
        totalvmolarV += b[k] * Z[k][1] / B[k] * xvector[3 + composition.size() + k] * composition[k].molefraction;
        //totalvmolarL += Z[k][0] * global.R * xvector[0] / (xvector[1] + 0.0000001) * xvector[3 + k];
        //totalvmolarV += Z[k][1] * global.R * xvector[0] / (xvector[1] + 0.0000001) * xvector[3 + composition.Count + k];
        massofonemole += composition[k].m->molarmass * composition[k].molefraction;
    }
    totalvmolar = (1 - xvector[2]) * totalvmolarL + xvector[2] * totalvmolarV;
    density->v = massofonemole / totalvmolar;
    return totalvmolar;
}

double material::dfdx(int fi, int xi)
{
    double y0 = fflash(fi);
    double x0 = xvector[xi];
    //xvector[xi] *= global.epsilonfrac;
    xvector[xi] += global::epsilonadd;
    double den = xvector[xi] - x0;
    double y1 = fflash(fi);
    xvector[xi] = x0;
    return (y1 - y0) / den;
}

double material::fflash(int i) //function that will define the functions to be solved for the UV flash
{
    double ffreturn = 0.0; //the value that will be returned.
    if (i == 0) //internal molar energy equation.
    {
        ffreturn = umolar.v - calcumolar(); //f: xvector[2]
        return ffreturn;
    }
    else if (i == 1) //molar volume equation
    {
        //vmolar = calcvmolar(); //We are nulling this equation just for this particular case now.
        ffreturn = vmolar.v - calcvmolaranddensity();
        return ffreturn;
    }
    else if (i >= 2 && i < composition.size() + 2) //fugacity equation
    {
        calcfugacity(i - 2);
        if (xvector[2] == 0.0) { ffreturn = 0.0; }
        else
        {
            ffreturn = xvector[3 + composition.size() + i - 2] * fugacityV[i - 2] -
                        xvector[3 + i - 2] * fugacityL[i - 2];
        } //y[j]*fl[j] - x[j]*fv[j] }
        return ffreturn;
    }
    else if (i >= composition.size() + 2 && i < 2 * composition.size() + 2)  //(1 - f)*x[j] + f*y[j] - z[j]
    {
        ffreturn = (1 - xvector[2]) * xvector[3 + (i - composition.size() - 2)] +
                    xvector[2] * xvector[3 + composition.size() + (i - composition.size() - 2)] - z[i - composition.size() - 2];
        return ffreturn;
    }
    else if (i == 2 * composition.size() + 2) //composition sum
    {
        double sum = 0;
        for (int k = 0; k < composition.size(); k++)
        {
            sum += xvector[3 + composition.size() + k] - xvector[3 + k];
        }
        return sum;
    }
    else { return 0; }
}

void material::calcjacobian()
{
    for (int r = 0; r < uvflashsize; r++)
    {
         for (int c = 0; c < uvflashsize; c++)
         {
              jacobian.m[r][c] = dfdx(r, c);
         }
    }
}

void material::limitxvector()
{
    if (xvector[0] < 0.0) { xvector[0] = 0.1; } //T
    if (xvector[1] < 0.0) { xvector[1] = 0.1; } //P
    if (xvector[2] < 0.0) { xvector[2] = 0.0; } //f
    else if (xvector[2] > 1.0) { xvector[2] = 1.0; }
}

void material::calccompz() //method to calculate the molar fraction of each component in the material
{
    for (int i = 0; i < composition.size(); i++)
    {
        z[i] = composition[i].molefraction;
    }
}

void material::uvflash()
{
    matrix lmatrix(uvflashsize, uvflashsize);
    matrix umatrix(uvflashsize, uvflashsize);
    matrix deltax(uvflashsize, 1);
    matrix ymatrix(uvflashsize, 1);
    matrix bmatrix(uvflashsize, 1);

    mapvarstox();
    //calcn();  This calcn will need to come in later when more than once component can be better simulated by the simulation.
    umolar.v = U.v / n->v;
    vmolar.v = V->v / n->v;

    calccompz();

    for (int i = 0; i < global::NMaterialIterations; i++)
    {
        if (xvector[2] == 0.0 || xvector[2] == 1.0)
        {
            uvflashsize = 2;
        }
        else { uvflashsize = origuvflashsize; }

        jacobian.resize(uvflashsize, uvflashsize);
        lmatrix.resize(uvflashsize, uvflashsize);
        umatrix.resize(uvflashsize, uvflashsize);
        deltax.resize(uvflashsize, 1);
        ymatrix.resize(uvflashsize, 1);
        bmatrix.resize(uvflashsize, 1);

        calcjacobian();
        //jacobian.swoprowsinthematrix(1, 3);
        //jacobian.ludecomposition(lmatrix, umatrix);
        //jacobian.swoprowsinthematrix(2, 4);
        jacobian.ludecomposition(lmatrix, umatrix);
        matrix tempm = lmatrix * umatrix;
        for (int j = 0; j < uvflashsize; j++)
        {
            bmatrix.m[j][0] = -fflash(j);
        }

        matrix::solveLYequalsB(lmatrix, ymatrix, bmatrix);
        matrix tempm2 = lmatrix * ymatrix;
        matrix::solveUXequalsY(umatrix, deltax, ymatrix);
        matrix tempm3 = umatrix * deltax;
        for (int j = 0; j < uvflashsize; j++)
        {
             xvector[j] += deltax.m[j][0];
        }
        limitxvector();
        f.v = xvector[2];
    }
    mapxtovars();
}

//public void addtothisobject(material m)
        //{
        //    for (int i = 0; i < composition.Length; i++)
        //    {
        //        composition[i].massfraction += m.composition[i].massfraction;
        //    }
        //}

void material::zero()
{
    for (int i = 0; i < composition.size(); i++)
    {
        composition[i].n = 0;
    }
}

void material::update(int simi, bool historise)
{
    if (historise && (simi % global::SimVectorUpdatePeriod == 0))
    {
        if (T->simvector.size() != 0) { T->simvector[simi/global::SimVectorUpdatePeriod] = T->v; }


        if (P->simvector.size() != 0) { P->simvector[simi / global::SimVectorUpdatePeriod] = P->v; }


        if (f.simvector.size() != 0) { f.simvector[simi / global::SimVectorUpdatePeriod] = f.v; }


        if (n->simvector.size() != 0) { n->simvector[simi / global::SimVectorUpdatePeriod] = n->v; }


        if (U.simvector.size() != 0) { U.simvector[simi / global::SimVectorUpdatePeriod] = U.v; }


        if (density->simvector.size() != 0) { density->simvector[simi / global::SimVectorUpdatePeriod] = density->v; }

    }
}













