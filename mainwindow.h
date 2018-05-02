#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QSettings>
#include <QFileDialog>
#include <QDragEnterEvent>
#include <QDropEvent>
/**
 * \file
 *
 * \author Valentin Bruder
 *
 * \copyright Copyright (C) 2018 Valentin Bruder
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

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
    void loadRawTff();
    void saveTff();
    void saveRawTff();

    void chooseBackgroundColor();
    void saveCamState();
    void loadCamState();
    void showAboutDialog();
protected:
    void dragEnterEvent(QDragEnterEvent *ev) Q_DECL_OVERRIDE;
    void dropEvent(QDropEvent *ev) Q_DECL_OVERRIDE;
    void closeEvent(QCloseEvent *event) Q_DECL_OVERRIDE;
    void keyPressEvent(QKeyEvent *event) Q_DECL_OVERRIDE;
    void showEvent(QShowEvent *event) Q_DECL_OVERRIDE;

private:

    void setVolumeData(const QString &fileName);
    bool readVolumeFile(const QString &fileName);

    void readTff(const QString &fileName);

    void readSettings();
    void writeSettings();

    void setStatusText();

    // ----- Members -----
    Ui::MainWindow *ui;

    QSettings *_settings;
    QFutureWatcher<void> *_watcher;
    QProgressBar _progBar;
    QTimer _timer;
    QString _fileName;
    QLabel _statusLabel;
};

#endif // MAINWINDOW_H
