#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "lbgstippling.h"

#include <QtWidgets>

class StippleViewer;

class MainWindow : public QMainWindow {

  public:
    MainWindow(QImage& density, LBGStippling::Params& params);

  private:
    StippleViewer* m_stippleViewer;
    QStatusBar* m_statusBar;
};

#endif // MAINWINDOW_H
