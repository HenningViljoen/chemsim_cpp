#include "drawwidget.h"

#include <QMouseEvent>
#include <QPaintEvent>
#include <QPainter>
#include <QRgb>

#include "global.h"

DrawWidget::DrawWidget(MainWindow *parent) : QWidget(parent)
{
    parent_window = parent;
    m_drawColor = QColor(Qt::black);
}

DrawWidget::~DrawWidget()
{

}

void DrawWidget::drawPixel(QPoint pt){
    QRgb value = m_drawColor.rgb();
    m_canvas.setPixel(pt.x(),pt.y(),value);
}

void DrawWidget::mousePressEvent(QMouseEvent *event){
    if (event->buttons() & Qt::LeftButton) {
        if (parent_window->drawwing_mode_entity == ValveMode) {
            parent_window->addnewvalve();
            //parent_window->sim->unitops
            drawPixel(event->pos());
            repaint();

        }
    }
}

void DrawWidget::mouseMoveEvent(QMouseEvent *event){

    if(event->buttons() & Qt::LeftButton){
        /*
        if (parent_window->drawwing_mode_entity == ValveMode) {
            drawPixel(event->pos());
            repaint();
        }
        */
    }
}

void DrawWidget::resizeEvent(QResizeEvent *event)
{
    m_canvas = QImage(width(),height(),QImage::Format_RGBA8888);
}

QColor DrawWidget::drawColor()
{
    return m_drawColor;
}

void DrawWidget::setDrawColor(QColor color)
{
    m_drawColor = color;
}

void DrawWidget::paintEvent(QPaintEvent *event)
{
    QWidget::paintEvent(event);
    QPainter painter(this);

    painter.drawPixmap(0,0,QPixmap::fromImage(m_canvas));
}
