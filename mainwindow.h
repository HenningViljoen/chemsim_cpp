#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "global.h"
#include "simulation.h"


namespace Ui {
class MainWindow;
}

class MainWindowPrivate;
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    DrawModeEntity drawwing_mode_entity;
    simulation *sim;

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void addnewvalve(int x, int y);

public slots:
    void colorPickTriggered();
    void valvebuttonclicked();

private:
    MainWindowPrivate *d;

};
#endif // MAINWINDOW_H
