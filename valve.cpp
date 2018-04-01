#include "valve.h"
#include <QPainter>
#include "global.h"

valve::valve(int anr, double ax, double ay, double aCv, double aop) : unitop{anr, ax, ay, 1, 1}
{
    int h = global::TimerInterval;
}

void valve::draw(QImage *graphics)
{
    //QRgb value = QColor(Qt::black).rgb();
    //graphics->setPixel(location->x, location->y, value);

    QPainter painter(graphics);
    painter.drawLine(80, 80, 400, 240);

    //painter.setRenderHint(QPainter::Antialiasing, true);
    //painter.setPen(QPen(Qt::black, 12, Qt::DashDotLine, Qt::RoundCap));
    //painter.setBrush(QBrush(Qt::green, Qt::SolidPattern));
    //painter.drawEllipse(80, 80, 400, 240);

}
