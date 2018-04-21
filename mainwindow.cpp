#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "drawwidget.h"
#include "valve.h"
#include <QColorDialog>
#include <QLayout>

#include <QToolBar>

class MainWindowPrivate {
public:
    MainWindowPrivate(MainWindow *owner);
    ~MainWindowPrivate();

    void pickColor();

private:
    DrawWidget *m_drawWidget;
    MainWindow *m_Owner;
    Ui::MainWindow *m_ui;
    //QToolBar* FileBar;
};

MainWindowPrivate::MainWindowPrivate(MainWindow *owner) :
    m_Owner(owner),
    m_drawWidget(new DrawWidget(owner)),
    m_ui(new Ui::MainWindow)
{
    m_ui->setupUi(m_Owner);
    m_ui->centralWidget->setLayout(new QVBoxLayout(m_Owner));
    m_ui->centralWidget->layout()->addWidget(m_drawWidget);

}

void MainWindowPrivate::pickColor()
{
    /*QColor color = QColorDialog::getColor(
                m_drawWidget->drawColor(),
                m_Owner,
                QString("Select a draw color"),
                QColorDialog::ShowAlphaChannel
                );
    m_drawWidget->setDrawColor(color);*/
}

MainWindowPrivate::~MainWindowPrivate(){
    delete m_ui;
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    d(new MainWindowPrivate(this))
{
    drawwing_mode_entity = DefaultDrawModeEntity;
    sim = new simulation();
}

MainWindow::~MainWindow()
{
    delete sim;
    delete d;
}

void MainWindow::addnewvalve(double x, double y)
{
    sim->unitops.push_back(new valve(0,x, y,1,0.5));
}

void MainWindow::colorPickTriggered()
{
    d->pickColor();
}

void MainWindow::valvebuttonclicked()
{
    drawwing_mode_entity = ValveMode;
}

