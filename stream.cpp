#include "stream.h"


stream::stream(int anr, double p0x, double p0y, double p1x, double p1y)
            : baseprocessclass(anr, (p0x + p1x) / 2, (p0y + p1y) / 2)
{
     streaminit(p0x, p0y, p1x, p1y);
}

stream::stream(stream *streamcopyfrom)
            : baseprocessclass(streamcopyfrom->nr, (streamcopyfrom->points[0].x + streamcopyfrom->points[1].x) / 2,
                (streamcopyfrom->points[0].y + streamcopyfrom->points[1].y) / 2)
{
     streaminit(streamcopyfrom->points[0].x, streamcopyfrom->points[0].y, streamcopyfrom->points[1].x,
                streamcopyfrom->points[1].y);
     copyfrom(streamcopyfrom);
}

void stream::streaminit(double p0x, double p0y, double p1x, double p1y)
{
     objecttype = StreamObjectType;
     name = std::to_string(nr) + " " + global::objecttypes_strings[objecttype];

     points.assign(2, point());
     points[0].setxy(p0x, p0y);
     points[1].setxy(p1x, p1y);

     //updatedirection(); NEED TO COMMENT OUT AFTER CODE PORT.
     //inbetweenpoints = new List<point>(0);
     displaypoints.assign(global::StreamNrPropDisplay, point());
     for (int i = 0; i < global::StreamNrPropDisplay; i++)
     {
          displaypoints[i].setxy(0, 0);
     }
     //update(0, false); NEED TO COMMENT OUT AFTER CODE PORT.
}

void stream::copyfrom(baseclass *baseclasscopyfrom)
{
     stream *streamcopyfrom = (stream *)baseclasscopyfrom;

     baseprocessclass::copyfrom(baseclasscopyfrom);

     direction = streamcopyfrom->direction; //radians
     distance = streamcopyfrom->distance; //m
}

void stream::updatedirection()
{
     direction = utilities::calcdirection(points[1].y - points[0].y, points[1].x - points[0].x);
     distance = utilities::distance(points[0], points[1]);
}

void stream::updatemassflowsource(std::string afilename)
{

}

void stream::updatepoint(int i, double x, double y)
{
     points[i].x = x;
     points[i].y = y;
     updatedirection();
}

void stream::update(int simi, bool historise)
{
     if (simi > 0)
     {
          // the reference flow is the mass flow
          calcactualvolumeflowfrommassflow();
          calcmolarflowfrommassflow();
          calcstandardflowfrommoleflow(); //should run after calcdndt since molar flow to be calculated first.
          if (massflow->datasource == Exceldata)
          {
               int i = (simi >= massflow->excelsource->data.size()) ? i = massflow->excelsource->data.size() - 1 : simi;
               massflow->v = massflow->excelsource->data[i];
          }

          if (mat.T->datasource == Exceldata) //This part should actually move to material.update().
          {
               int i = (simi >= mat.T->excelsource->data.size()) ? i = mat.T->excelsource->data.size() - 1 : simi;
               mat.T->v = mat.T->excelsource->data[i];
          }
          if (mat.relativehumidity.datasource == Exceldata) //This part should actually move to material.update().
          {
               int i = (simi >= mat.relativehumidity.excelsource->data.size()) ? i = mat.relativehumidity.excelsource->data.size() - 1 : simi;
               mat.relativehumidity.v = mat.relativehumidity.excelsource->data[i];
          }
     }

     if (historise && (simi % global::SimVectorUpdatePeriod == 0))
     {

          if (actualvolumeflow->simvector.size() != 0) { actualvolumeflow->simvector[simi/global::SimVectorUpdatePeriod] = actualvolumeflow->v; }

          if (standardvolumeflow->simvector.size() != 0) { standardvolumeflow->simvector[simi / global::SimVectorUpdatePeriod] =
                        standardvolumeflow->v; }

          if (massflow->simvector.size() != 0) { massflow->simvector[simi / global::SimVectorUpdatePeriod] = massflow->v; }

          if (molarflow->simvector.size() != 0) { molarflow->simvector[simi / global::SimVectorUpdatePeriod] = molarflow->v; }
                //pressuresimvector[simi] = mat.P.v;
      }
      baseprocessclass::update(simi, historise);
}


















