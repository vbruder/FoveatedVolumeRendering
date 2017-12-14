#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QSettings>
#include <QFileDialog>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QMimeData>
#include <QFutureWatcher>
#include <QProgressBar>
#include <QTimer>
#include <QGraphicsScene>
#include <QGraphicsLineItem>
#include <QLabel>
#include <QSettings>
#include <QVector4D>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

protected slots:
    void openVolumeFile();

    void addProgress();
    void finishedLoading();

    void loadTff();
    void saveTff();
    void saveRawTff();

    void chooseBackgroundColor();
protected:
    void dragEnterEvent(QDragEnterEvent *ev) Q_DECL_OVERRIDE;
    void dropEvent(QDropEvent *ev) Q_DECL_OVERRIDE;
    void closeEvent(QCloseEvent *event) Q_DECL_OVERRIDE;

private:

    void setVolumeData(const QString &fileName);
    bool readVolumeFile(const QString &fileName);

    void readTff(const QString &fileName);

    void readSettings();
    void writeSettings();

    // ----- Members -----
    Ui::MainWindow *ui;

    QSettings *_settings;
    QFutureWatcher<void> *_watcher;
    QProgressBar _progBar;
    QTimer _timer;
    QString _fileName;
    QLabel *_statusLabel;
};

#endif // MAINWINDOW_H
