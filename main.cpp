#include "mainwindow.h"
#include <QApplication>

//#include "global.h"

int main(int argc, char *argv[])
{
    //global_class global = global_class();
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
