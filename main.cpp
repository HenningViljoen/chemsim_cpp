#include "mainwindow.h"
#include <QApplication>

#include "exceldataset.h"

int main(int argc, char *argv[])
{
    //exceldataset excel("/Users/johanneshenningviljoen/Dropbox/Projects-DB/ChemSim/Electronic_data/135-ES-004/170728/135TI0015.PV.csv");
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
