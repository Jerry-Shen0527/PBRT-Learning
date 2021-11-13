#include "Raytracing.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Raytracing w;
    w.show();
    return a.exec();
}
