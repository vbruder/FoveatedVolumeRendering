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

#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFileDialog>
#include <QMessageBox>
#include <QDir>
#include <QString>
#include <QUrl>
#include <QFuture>
#include <QtConcurrentRun>
#include <QThread>
#include <QColorDialog>
#include <QtGlobal>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QJsonObject>
#include <QJsonDocument>

/**
 * @brief MainWindow::MainWindow
 * @param parent
 */
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent)
  , ui(new Ui::MainWindow)
  , _fileName("No volume data loaded yet.")
{
    QCoreApplication::setOrganizationName("VISUS");
    QCoreApplication::setOrganizationDomain("www.visus.uni-stuttgart.de");
    QCoreApplication::setApplicationName("Lightweight Volume Raycaster");
    _settings = new QSettings();

    setAcceptDrops( true );
    ui->setupUi(this);
    ui->gbTimeSeries->setVisible(false);
    connect(ui->volumeRenderWidget, &VolumeRenderWidget::timeSeriesLoaded,
            ui->gbTimeSeries, &QGroupBox::setVisible);
    connect(ui->volumeRenderWidget, &VolumeRenderWidget::timeSeriesLoaded,
            ui->sbTimeStep, &QSpinBox::setMaximum);
    connect(ui->volumeRenderWidget, &VolumeRenderWidget::timeSeriesLoaded,
            ui->sldTimeStep, &QSlider::setMaximum);
    connect(ui->sldTimeStep, &QSlider::valueChanged,
            ui->volumeRenderWidget, &VolumeRenderWidget::setTimeStep);

    // menu bar actions
    connect(ui->actionOpen, &QAction::triggered, this, &MainWindow::openVolumeFile);
    connect(ui->actionSaveCpTff, &QAction::triggered, this, &MainWindow::saveTff);
    connect(ui->actionSaveRawTff_2, &QAction::triggered, this, &MainWindow::saveRawTff);
    connect(ui->actionLoadCpTff, &QAction::triggered, this, &MainWindow::loadTff);
    connect(ui->actionLoadRawTff, &QAction::triggered, this, &MainWindow::loadRawTff);
    connect(ui->actionScreenshot, &QAction::triggered,
            ui->volumeRenderWidget, &VolumeRenderWidget::saveFrame);
    connect(ui->actionRecord, &QAction::triggered,
            ui->volumeRenderWidget, &VolumeRenderWidget::toggleVideoRecording);
    connect(ui->actionGenerateLowResVo, &QAction::triggered,
            ui->volumeRenderWidget, &VolumeRenderWidget::generateLowResVolume);
    connect(ui->actionResetCam, &QAction::triggered,
            ui->volumeRenderWidget, &VolumeRenderWidget::resetCam);
    connect(ui->actionSaveState, &QAction::triggered, this, &MainWindow::saveCamState);
    connect(ui->actionLoadState, &QAction::triggered, this, &MainWindow::loadCamState);
    connect(ui->actionShowOverlay, &QAction::toggled,
            ui->volumeRenderWidget, &VolumeRenderWidget::setShowOverlay);
    connect(ui->actionSelectOpenCL, &QAction::triggered,
            ui->volumeRenderWidget, &VolumeRenderWidget::showSelectOpenCL);
    connect(ui->actionAbout, &QAction::triggered, this, &MainWindow::showAboutDialog);


    // future watcher for concurrent data loading
    _watcher = new QFutureWatcher<void>(this);
    connect(_watcher, &QFutureWatcher<void>::finished, this, &MainWindow::finishedLoading);
    // loading progress bar
    _progBar.setRange(0, 100);
    _progBar.setTextVisible(true);
    _progBar.setAlignment(Qt::AlignCenter);
    connect(&_timer, &QTimer::timeout, this, &MainWindow::addProgress);

    // connect settings UI
    connect(ui->dsbSamplingRate, 
            static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
            ui->volumeRenderWidget, &VolumeRenderWidget::updateSamplingRate);
    connect(ui->dsbImgSampling, 
            static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
            ui->volumeRenderWidget, &VolumeRenderWidget::setImageSamplingRate);
    connect(ui->cbIllum, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            ui->volumeRenderWidget, &VolumeRenderWidget::setIllumination);
    connect(ui->pbBgColor, &QPushButton::released, this, &MainWindow::chooseBackgroundColor);
    // check boxes
    connect(ui->chbLinear, &QCheckBox::toggled,
            ui->volumeRenderWidget, &VolumeRenderWidget::setLinearInterpolation);
    connect(ui->chbAmbientOcclusion, &QCheckBox::toggled,
            ui->volumeRenderWidget, &VolumeRenderWidget::setAmbientOcclusion);
    connect(ui->chbContours, &QCheckBox::toggled,
            ui->volumeRenderWidget, &VolumeRenderWidget::setContours);
    connect(ui->chbAerial, &QCheckBox::toggled,
            ui->volumeRenderWidget, &VolumeRenderWidget::setAerial);
    connect(ui->chbBox, &QCheckBox::toggled,
            ui->volumeRenderWidget, &VolumeRenderWidget::setDrawBox);
    connect(ui->chbOrtho, &QCheckBox::toggled,
            ui->volumeRenderWidget, &VolumeRenderWidget::setCamOrtho);
    // connect tff editor
    connect(ui->transferFunctionEditor->getEditor(), &TransferFunctionEditor::gradientStopsChanged,
            ui->volumeRenderWidget, &VolumeRenderWidget::updateTransferFunction);
    connect(ui->pbResetTff, &QPushButton::clicked,
            ui->transferFunctionEditor, &TransferFunctionWidget::resetTransferFunction);
    connect(ui->cbInterpolation, 
            static_cast<void (QComboBox::*)(const QString &)>(&QComboBox::currentIndexChanged),
            ui->volumeRenderWidget, &VolumeRenderWidget::setTffInterpolation);
    connect(ui->cbInterpolation, SIGNAL(currentIndexChanged(QString)),
            ui->transferFunctionEditor, SLOT(setInterpolation(QString)));
//    connect(ui->cbInterpolation, qOverload<const QString &>(&QComboBox::currentIndexChanged),
//            ui->transferFunctionEditor, &TransferFunctionEditor::setInterpolation);

    connect(ui->transferFunctionEditor->getEditor(), &TransferFunctionEditor::selectedPointChanged,
            ui->colorWheel, &colorwidgets::ColorWheel::setColor);
    connect(ui->colorWheel, &colorwidgets::ColorWheel::colorChanged,
            ui->transferFunctionEditor, &TransferFunctionWidget::setColorSelected);

    ui->statusBar->addPermanentWidget(&_statusLabel);
    connect(ui->volumeRenderWidget, &VolumeRenderWidget::frameSizeChanged,
            this, &MainWindow::setStatusText);

    // restore settings
    readSettings();
}


/**
 * @brief MainWindow::~MainWindow
 */
MainWindow::~MainWindow()
{
    delete _watcher;
    delete _settings;
    delete ui;
}


/**
 * @brief MainWindow::closeEvent
 * @param event
 */
void MainWindow::closeEvent(QCloseEvent *event)
{
    writeSettings();
    event->accept();
}


void MainWindow::keyPressEvent(QKeyEvent *event)
{
    float factor = 0.01f;
    if (event->modifiers() & Qt::ShiftModifier)
        factor = 0.001f;

    switch (event->key())
    {
        case Qt::Key_Up:
        case Qt::Key_W: ui->volumeRenderWidget->updateView(+0.000f,-factor); break;
        case Qt::Key_Left:
        case Qt::Key_A: ui->volumeRenderWidget->updateView(-factor, 0.000f); break;
        case Qt::Key_Down:
        case Qt::Key_S: ui->volumeRenderWidget->updateView(+0.000f,+factor); break;
        case Qt::Key_Right:
        case Qt::Key_D: ui->volumeRenderWidget->updateView(+factor, 0.000f); break;
        // TODO: zoom
    }
    event->accept();
}


/**
 * @brief MainWindow::showEvent
 * @param event
 */
void MainWindow::showEvent(QShowEvent *event)
{
    ui->transferFunctionEditor->resetTransferFunction();
    event->accept();
}


/**
 * @brief MainWindow::writeSettings
 */
void MainWindow::writeSettings()
{
    _settings->beginGroup("MainWindow");
    _settings->setValue("geometry", saveGeometry());
    _settings->setValue("windowState", saveState());
    _settings->endGroup();

    _settings->beginGroup("Settings");
    // TODO
    _settings->endGroup();
}


/**
 * @brief MainWindow::readSettings
 */
void MainWindow::readSettings()
{
    _settings->beginGroup("MainWindow");
    restoreGeometry(_settings->value("geometry").toByteArray());
    restoreState(_settings->value("windowState").toByteArray());
    _settings->endGroup();

    _settings->beginGroup("Settings");
    // todo
    _settings->endGroup();
}


/**
 * @brief MainWindow::setVolumeData
 * @param fileName
 */
void MainWindow::setVolumeData(const QString &fileName)
{
    ui->volumeRenderWidget->setVolumeData(fileName);
    ui->volumeRenderWidget->updateTransferFunction(
                ui->transferFunctionEditor->getEditor()->getGradientStops());
    ui->volumeRenderWidget->updateView();
}


/**
 * @brief MainWindow::loadCamState
 */
void MainWindow::loadCamState()
{
    QFileDialog dialog;
    QString defaultPath = _settings->value( "LastStateFile" ).toString();
    QString pickedFile = dialog.getOpenFileName(this, tr("Save State"),
                                                defaultPath, tr("JSON files (*.json)"));
    if (pickedFile.isEmpty())
        return;
    _settings->setValue( "LastStateFile", pickedFile );

    QFile loadFile(pickedFile);
    if (!loadFile.open(QIODevice::ReadOnly))
    {
        qWarning() << "Couldn't open state file" << pickedFile;
        return;
    }

    QByteArray saveData = loadFile.readAll();
    QJsonDocument loadDoc(QJsonDocument::fromJson(saveData));
    QJsonObject json = loadDoc.object();

    if (json.contains("imgResFactor") && json["imgResFactor"].isDouble())
            ui->dsbImgSampling->setValue(json["imgResFactor"].toDouble());
    if (json.contains("rayStepSize") && json["rayStepSize"].isDouble())
            ui->dsbSamplingRate->setValue(json["rayStepSize"].toDouble());

    if (json.contains("useLerp") && json["useLerp"].isBool())
            ui->chbLinear->setChecked(json["useLerp"].toBool());
    if (json.contains("useAO") && json["useAO"].isBool())
            ui->chbAmbientOcclusion->setChecked(json["useAO"].toBool());
    if (json.contains("showContours") && json["showContours"].isBool())
            ui->chbContours->setChecked(json["showContours"].toBool());
    if (json.contains("useAerial") && json["useAerial"].isBool())
            ui->chbAerial->setChecked(json["useAerial"].toBool());
    if (json.contains("showBox") && json["showBox"].isBool())
            ui->chbBox->setChecked(json["showBox"].toBool());
    if (json.contains("useOrtho") && json["useOrtho"].isBool())
            ui->chbOrtho->setChecked(json["useOrtho"].toBool());
    // camera paramters
    ui->volumeRenderWidget->read(json);
}

/**
 * @brief MainWindow::showAboutDialog
 */
void MainWindow::showAboutDialog()
{
    QMessageBox::about(this, "About Volume Raycaster",
"<b>Volume Raycaster 2018</b><br><br>\
Check out the \
<a href='https://bitbucket.org/theVall/basicvolumeraycaster/overview'>Bitbucket repository</a> \
for more information.<br><br>\
Copyright 2017-2018 Valentin Bruder. All rights reserved. <br><br>\
The program is provided AS IS with NO WARRANTY OF ANY KIND, \
INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS \
FOR A PARTICULAR PURPOSE.");
}

/**
 * @brief MainWindow::saveCamState
 */
void MainWindow::saveCamState()
{
    QFileDialog dialog;
    QString defaultPath = _settings->value( "LastStateFile" ).toString();
    QString pickedFile = dialog.getSaveFileName(this, tr("Save State"),
                                                defaultPath, tr("JSON files (*.json)"));
    if (pickedFile.isEmpty())
        return;
    _settings->setValue( "LastStateFile", pickedFile );

    QFile saveFile(pickedFile);
    if (!saveFile.open(QIODevice::WriteOnly))
    {
        qWarning() << "Couldn't open save file" << pickedFile;
        return;
    }

    QJsonObject stateObject;
    // resolution
    stateObject["imgResFactor"] = ui->dsbImgSampling->value();
    stateObject["rayStepSize"] = ui->dsbSamplingRate->value();
    // rendering flags
    stateObject["useLerp"] = ui->chbLinear->isChecked();
    stateObject["useAO"] = ui->chbAmbientOcclusion->isChecked();
    stateObject["showContours"] = ui->chbContours->isChecked();
    stateObject["useAerial"] = ui->chbAerial->isChecked();
    stateObject["showBox"] = ui->chbBox->isChecked();
    stateObject["useOrtho"] = ui->chbOrtho->isChecked();
    // camera parameters
    ui->volumeRenderWidget->write(stateObject);

    QJsonDocument saveDoc(stateObject);
    saveFile.write(saveDoc.toJson());
}


/**
 * @brief MainWindow::readVolumeFile
 * @param fileName
 * @return
 */
bool MainWindow::readVolumeFile(const QString &fileName)
{
    ui->volumeRenderWidget->setLoadingFinished(false);
    _progBar.setFormat("Loading volume file: " + fileName);
    _progBar.setValue(1);
    _progBar.show();
    ui->statusBar->addPermanentWidget(&_progBar, 2);
    ui->statusBar->updateGeometry();
    QApplication::processEvents();

    if (fileName.isEmpty())
    {
        _progBar.deleteLater();
        throw std::invalid_argument("Invalid volume data file name.");
    }
    else
    {
        _fileName = fileName;
        QFuture<void> future = QtConcurrent::run(this, &MainWindow::setVolumeData, fileName);
        _watcher->setFuture(future);
        _timer.start(100);
    }
    return true;
}


/**
 * @brief MainWindow::readTff
 * @param fileName
 * @return
 */
void MainWindow::readTff(const QString &fileName)
{
    if (fileName.isEmpty())
    {
        throw std::invalid_argument("Invalid trtansfer function file name.");
    }
    else
    {
        QFile file(fileName);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            throw std::invalid_argument("Could not open transfer function file "
                                        + fileName.toStdString());
        }
        else
        {
            QTextStream in(&file);
            QGradientStops stops;
            while (!in.atEnd())
            {
                QStringList line = in.readLine().split(QRegExp("\\s"));
                if (line.size() < 5)
                    continue;
                QGradientStop stop(line.at(0).toDouble(),
                                   QColor(line.at(1).toInt(), line.at(2).toInt(),
                                          line.at(3).toInt(), line.at(4).toInt()));
                stops.push_back(stop);
            }
            if (!stops.isEmpty())
            {
                ui->transferFunctionEditor->getEditor()->setGradientStops(stops);
                ui->transferFunctionEditor->getEditor()->pointsUpdated();
            }
            else
                qCritical() << "Empty transfer function file.";
            file.close();
        }
    }
}


/**
 * @brief MainWindow::saveTff
 * @param fileName
 */
void MainWindow::saveTff()
{
    QFileDialog dia;
    QString defaultPath = _settings->value( "LastTffFile" ).toString();
    QString pickedFile = dia.getSaveFileName(
                this, tr("Save Transfer Function"),
                defaultPath, tr("Transfer function files (*.tff)"));

    if (!pickedFile.isEmpty())
    {
        QFile file(pickedFile);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            throw std::invalid_argument("Could not open file "
                                        + pickedFile.toStdString());
        }
        else
        {
            QTextStream out(&file);
            const QGradientStops stops =
                    ui->transferFunctionEditor->getEditor()->getGradientStops();
            foreach (QGradientStop s, stops)
            {
                out << s.first << " " << s.second.red() << " " << s.second.green()
                    << " " << s.second.blue() << " " << s.second.alpha() << "\n";
            }
            file.close();
        }
    }
}

/**
 * @brief MainWindow::saveRawTff
 */
void MainWindow::saveRawTff()
{
    QFileDialog dia;
    QString defaultPath = _settings->value( "LastRawTffFile" ).toString();
    QString pickedFile = dia.getSaveFileName(
                this, tr("Save Transfer Function"),
                defaultPath, tr("Transfer function files (*.tff)"));

    if (!pickedFile.isEmpty())
    {
        QFile file(pickedFile);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            throw std::invalid_argument("Could not open file "
                                        + pickedFile.toStdString());
        }
        else
        {
            QTextStream out(&file);
            const QGradientStops stops =
                    ui->transferFunctionEditor->getEditor()->getGradientStops();
            const std::vector<unsigned char> tff = ui->volumeRenderWidget->getRawTransferFunction(stops);
            foreach (unsigned char c, tff)
            {
                out << (int)c << " ";
            }
            file.close();
        }
    }
}


/**
 * @brief MainWindow::saveRawTff
 */
void MainWindow::loadRawTff()
{
    QFileDialog dia;
    QString defaultPath = _settings->value( "LastRawTffFile" ).toString();
    QString pickedFile = dia.getOpenFileName(
                this, tr("Open Transfer Function"),
                defaultPath, tr("Transfer function files (*.tff)"));
    if (!pickedFile.isEmpty())
    {
        qDebug() << "Loading transfer funtion data defined in" << pickedFile;
        std::ifstream tff_file(pickedFile.toStdString(), std::ios::in);
        float value = 0;
        std::vector<unsigned char> values;

        // read lines from file and split on whitespace
        if (tff_file.is_open())
        {
            while (tff_file >> value)
            {
                values.push_back((char)value);
            }
            tff_file.close();
            ui->volumeRenderWidget->setRawTransferFunction(values);
        }
        else
        {
            qDebug() << "Could not open transfer function file " + pickedFile;
        }
        _settings->setValue( "LastRawTffFile", pickedFile );
    }
}


/**
 * @brief MainWindow::getStatus
 * @return
 */
void MainWindow::setStatusText()
{
    QString status = "No data loaded yet.";
    if (ui->volumeRenderWidget->hasData())
    {
        status = "File: ";
        status += _fileName;
        status += " | Volume: ";
        status += QString::number(ui->volumeRenderWidget->getVolumeResolution().x());
        status += "x";
        status += QString::number(ui->volumeRenderWidget->getVolumeResolution().y());
        status += "x";
        status += QString::number(ui->volumeRenderWidget->getVolumeResolution().z());
        status += " | Frame: ";
        status += QString::number(ui->volumeRenderWidget->size().width());
        status += "x";
        status += QString::number(ui->volumeRenderWidget->size().height());
        status += " ";
    }
    _statusLabel.setText(status);
}


/**
 * @brief MainWindow::finishedLoading
 */
void MainWindow::finishedLoading()
{
    _progBar.setValue(100);
    _progBar.hide();
    _timer.stop();
//    qDebug() << "finished.";
    this->setStatusText();
    ui->volumeRenderWidget->setLoadingFinished(true);
    ui->volumeRenderWidget->updateView();
}


/**
 * @brief MainWindow::addProgress
 */
void MainWindow::addProgress()
{
    if (_progBar.value() < _progBar.maximum() - 5)
    {
        _progBar.setValue(_progBar.value() + 1);
    }
}


/**
 * @brief MainWindow::loadTff
 */
void MainWindow::loadTff()
{
    QFileDialog dia;
    QString defaultPath = _settings->value( "LastTffFile" ).toString();
    QString pickedFile = dia.getOpenFileName(
                this, tr("Open Transfer Function"),
                defaultPath, tr("Transfer function files (*.tff)"));
    if (!pickedFile.isEmpty())
    {
        readTff(pickedFile);
        _settings->setValue( "LastTffFile", pickedFile );
    }
}


/**
 * @brief MainWindow::openVolumeFile
 */
void MainWindow::openVolumeFile()
{
    QFileDialog dia;
    QString defaultPath = _settings->value( "LastVolumeFile" ).toString();
    QString pickedFile = dia.getOpenFileName(
                this, tr("Open Volume Data"), defaultPath, tr("Volume data files (*.dat)"));
    if (!pickedFile.isEmpty())
    {
        if (!readVolumeFile(pickedFile))
        {
            QMessageBox msgBox;
            msgBox.setIcon( QMessageBox::Critical );
            msgBox.setText( "Error while trying to create OpenCL memory objects." );
            msgBox.exec();
        }
        else
        {
            _settings->setValue( "LastVolumeFile", pickedFile );
        }
    }
}


/**
 * @brief MainWindow::dragEnterEvent
 * @param ev
 */
void MainWindow::dragEnterEvent(QDragEnterEvent *ev)
{
    if (ev->mimeData()->hasUrls())
    {
        bool valid = false;
        foreach(QUrl url, ev->mimeData()->urls())
        {
            if (!url.fileName().isEmpty())
            {
                QString sFileEnding = url.fileName().split(".", QString::SkipEmptyParts).at(1);
                if (sFileEnding == "dat")
                {
                    valid = true;
                }
            }
        }
        if (valid)
        {
            ev->acceptProposedAction();
        }
    }
}


/**
 * @brief MainWindow::dropEvent
 * @param ev
 */
void MainWindow::dropEvent(QDropEvent *ev)
{
    foreach(QUrl url, ev->mimeData()->urls())
    {
        QString sFileEnding = url.fileName().split(".", QString::SkipEmptyParts).at(1);
        if (sFileEnding == "dat")
        {
            // extract path and remove leading '/'
            QString fileName = url.path(); //.remove( 0, 1 );
            qDebug() << "Loading volume data file" << fileName;
            readVolumeFile(fileName);
        }
    }
}


/**
 * @brief MainWindow::chooseBackgroundColor
 */
void MainWindow::chooseBackgroundColor()
{
    QColorDialog dia;
    QColor col = dia.getColor();
    if (col.isValid())
        ui->volumeRenderWidget->setBackgroundColor(col);
}
