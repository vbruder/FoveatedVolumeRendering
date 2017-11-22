#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

//    // Load an application style
//    QFile styleFile( "../RaycastLight/style.qss" );
//    styleFile.open( QFile::ReadOnly );
//    // Apply the loaded stylesheet
//    QString style = QLatin1String( styleFile.readAll() );
//    a.setStyleSheet( style );

    MainWindow w;
//    w.ensurePolished();
    w.show();

    return a.exec();
}
