#ifndef DRAWWIDGET_H
#define DRAWWIDGET_H

#include <QObject>
#include <QWidget>

#include "mainwindow.h"

class QPaintEvent;
class QMouseEvent;

class DrawWidget : public QWidget
{
    Q_OBJECT
public:
    MainWindow *parent_window;

    explicit DrawWidget(MainWindow *parent = nullptr);
    ~DrawWidget();

    void drawPixel(QPoint pt);

signals:

protected:
    void paintEvent(QPaintEvent *);
    void mousePressEvent(QMouseEvent *);
    void mouseMoveEvent(QMouseEvent *);
    void resizeEvent(QResizeEvent *);

private:
    QColor m_drawColor;
    QPixmap m_pixmap;
    QImage m_canvas;

public slots:
    QColor drawColor();
    void setDrawColor(QColor color);
};

#endif // DRAWWIDGET_H
