#include "valve.h"
#include <iostream>
#include "global.h"

valve::valve(int anr, double ax, double ay, double aCv, double aop) : unitop{anr, ax, ay, 1, 1}
{
    int h = global::TimerInterval;
    std::cout << h;
    initvalve(anr, ax, ay, aCv, aop);
}

valve::valve(baseclass *baseclasscopyfrom)
    : unitop{0, 0, 0, 1, 1} //these numbes do not matter much since they will be sorted anyway by the copymethods down the
                                    //hierarchy
{
    initvalve(0, 0, 0, 0, 0);
    copyfrom(baseclasscopyfrom);
}

valve::~valve()
{
    if (deltapressure != nullptr) {delete deltapressure;}
    if (op != nullptr) {delete op;}
}

void valve::initvalve(int anr, double ax, double ay, double aCv, double aop)
{
    objecttype = Valve;

    deltapressure = new controlvar();
    op = new controlvar();

    name = std::to_string(nr) + " " + global::objecttypes_strings[objecttype];

    controlpropthisclass.clear();
    std::vector<std::string> newrange = {"deltapressure",
                                         "op"};
    controlpropthisclass.insert(std::end(controlpropthisclass), std::begin(newrange), std::end(newrange));
    nrcontrolpropinherited = controlproperties.size();
    controlproperties.insert(std::end(controlproperties), std::begin(controlpropthisclass), std::end(controlpropthisclass));

    Cv = aCv;
    op->v = aop;
    deltapressure->v = 0;
    //deltapressuresimvector = new double[global.SimVectorLength];
    //opsimvector = new double[global.SimVectorLength];

    actualvolumeflow->v = 0;

    deltapressurenew = 0; //Pa
    ddeltapressuredt = 0;

    updateinoutpointlocations();

    update(0, false);
}

void valve::copyfrom(baseclass *baseclasscopyfrom)
{
    valve *valvecopyfrom = (valve *)baseclasscopyfrom;

    unitop::copyfrom(baseclasscopyfrom);

    deltapressure->v = valvecopyfrom->deltapressure->v; //Pa
    Cv = valvecopyfrom->Cv; //Valve coeficient
    op->v = valvecopyfrom->op->v; //Fraction : valve opening as a fraction

    deltapressurenew = valvecopyfrom->deltapressurenew; //Pa
    ddeltapressuredt = valvecopyfrom->ddeltapressuredt;
}

controlvar *valve::selectedproperty(int selection)
{
    if (selection >= nrcontrolpropinherited)
    {
         switch (selection - nrcontrolpropinherited)
         {
              case 0:
                  return deltapressure;
              case 1:
                  return op;
              default:
                  return nullptr;
         }
    }
    else { return baseprocessclass::selectedproperty(selection); };
}

void valve::ddt(int simi)
{
    ddeltapressuredt = -1 / global::ValveHydraulicTau * deltapressure->v + 1 / global::ValveHydraulicTau * deltapressurenew;
}

void valve::update(int simi, bool historise)
{
    if (inflow[0] != nullptr && outflow[0] != nullptr)
    {
         mat.copycompositiontothisobject(inflow[0]->mat);
         mat.density->v = inflow[0]->mat.density->v;
         mat.T->v = inflow[0]->mat.T->v;
         massflow->v = inflow[0]->massflow->v;
         actualvolumeflow->v = massflow->v/mat.density->v;
         deltapressurenew = pow(actualvolumeflow->v / (Cv * pow(global::ValveEqualPercR, op->v - 1)), 2);

         ddt(simi);

         deltapressure->v += ddeltapressuredt * global::SampleT;

         inflow[0]->mat.P->v = outflow[0]->mat.P->v + deltapressure->v;

         calcmolarflowfrommassflow();
         calcstandardflowfrommoleflow();

         outflow[0]->mat.copycompositiontothisobject(mat);
         outflow[0]->massflow->v = massflow->v;
         outflow[0]->mat.density->v = mat.density->v;
         outflow[0]->mat.T->v = mat.T->v;
    }

    //if (op.v > 1) { op.v = 1; }
    if (op->v < 0.00) { op->v = 0.00; }

    if (historise && (simi % global::SimVectorUpdatePeriod == 0))
    {
         deltapressure->simvector.push_back(deltapressure->v);

         op->simvector.push_back(op->v);

         actualvolumeflow->simvector.push_back(actualvolumeflow->v);

         standardvolumeflow->simvector.push_back(standardvolumeflow->v);

         massflow->simvector.push_back(massflow->v);

         molarflow->simvector.push_back(molarflow->v);
    }
}

        //public void sizevalvefromstandardflow()
        //{
        //    calcmassflowfromstandardflow();
        //    calcactualvolumeflowfrommassflow();
        //    if (Math.Abs(deltapressure) > 0) { Cv = actualvolumeflow / (op * Math.Sqrt(Math.Abs(deltapressure))); }
        //}

void valve::sizevalvefromactualvolumeflow()
{
    if (abs(deltapressure->v) > 0) { Cv = actualvolumeflow->v / (op->v * sqrt(abs(deltapressure->v))); }
}

bool valve::mouseover(double x, double y)
{
    return (utilities::distance(x - location.x, y - location.y) <= global::ValveLength);
}

void valve::updateinoutpointlocations()
{
    inpoint[0].x = location.x - global::ValveLength / 2 - global::InOutPointWidth;
    inpoint[0].y = location.y;
    outpoint[0].x = location.x + global::ValveLength / 2 + global::InOutPointWidth;
    outpoint[0].y = location.y;

    unitop::updateinoutpointlocations();
}

void valve::setproperties(simulation *asim)
{/* THIS METHOD NEEDS TO WORK
    update(asim->simi, false);
    valveproperties valveprop = new valveproperties(this, asim);
    valveprop.Show(); */
}

void valve::draw(QImage *graphics)
{
    updateinoutpointlocations();

    QPolygon poly;
    poly << QPoint(global::OriginX + round(global::GScale*(location.x - global::ValveLength/2)),
                  global::OriginY + (global::GScale*(location.y - global::ValveWidth/2)))
         << QPoint(global::OriginX + round(global::GScale*(location.x + global::ValveLength / 2)),
                   global::OriginY + round(global::GScale*(location.y + global::ValveWidth/2)))
         << QPoint(global::OriginX + round(global::GScale*(location.x + global::ValveLength / 2)),
                   global::OriginY + round(global::GScale*(location.y - global::ValveWidth/2)))
         << QPoint(global::OriginX + round(global::GScale*(location.x - global::ValveLength / 2)),
                   global::OriginY + round(global::GScale*(location.y + global::ValveWidth/2)));

    QPainter painter(graphics);
    painter.drawPolygon(poly);
    if (highlighted)
    {
        QBrush fillbrush;
        fillbrush.setColor(Qt::red);
        fillbrush.setStyle(Qt::SolidPattern);

        QPainterPath path;
        path.addPolygon(poly);

        painter.fillPath(path,fillbrush);
    }
    unitop::draw(graphics);

    //QRgb value = QColor(Qt::black).rgb();
    //graphics->setPixel(location->x, location->y, value);
    //painter.drawLine(global::OriginX + static_cast<int>(global::GScale*location.x),
      //               global::OriginY + static_cast<int>(global::GScale*location.y), 400, 20);
    //painter.setRenderHint(QPainter::Antialiasing, true);
    //painter.setPen(QPen(Qt::black, 12, Qt::DashDotLine, Qt::RoundCap));
    //painter.setBrush(QBrush(Qt::green, Qt::SolidPattern));
    //painter.drawEllipse(80, 80, 400, 240);
}
