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

/**
 * @brief MainWindow::MainWindow
 * @param parent
 */
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    QCoreApplication::setOrganizationName("VISUS");
    QCoreApplication::setOrganizationDomain("www.visus.uni-stuttgart.de");
    QCoreApplication::setApplicationName("Lightweight Volume Raycaster");
    _settings = new QSettings();

    setAcceptDrops( true );
    ui->setupUi(this);

    // menu bar actions
    connect(ui->actionOpen, &QAction::triggered, this, &MainWindow::openVolumeFile);
    connect(ui->actionSaveTff, &QAction::triggered, this, &MainWindow::saveTff);
    connect(ui->actionLoadTff, &QAction::triggered, this, &MainWindow::loadTff);

    // future watcher for concurrent data loading
    _watcher = new QFutureWatcher<void>(this);
    connect(_watcher, &QFutureWatcher<void>::finished, this, &MainWindow::finishedLoading);
    // loading progress bar
    _progBar.setRange(0, 100);
    _progBar.setTextVisible(true);
    _progBar.setAlignment(Qt::AlignCenter);
    connect(&_timer, &QTimer::timeout, this, &MainWindow::addProgress);

    // connect settings UI
    connect(ui->dsbSamplingRate, qOverload<double>(&QDoubleSpinBox::valueChanged),
            ui->volumeRenderWidget, &VolumeRenderWidget::updateSamplingRate);
    connect(ui->chbLinear, &QCheckBox::toggled,
            ui->volumeRenderWidget, &VolumeRenderWidget::setLinearInterpolation);
    connect(ui->chbIllum, &QCheckBox::toggled,
            ui->volumeRenderWidget, &VolumeRenderWidget::setIllumination);
    connect(ui->chbBox, &QCheckBox::toggled,
            ui->volumeRenderWidget, &VolumeRenderWidget::setDrawBox);
    connect(ui->chbOrtho, &QCheckBox::toggled,
            ui->volumeRenderWidget, &VolumeRenderWidget::setCamOrtho);
    connect(ui->pbBgColor, &QPushButton::released, this, &MainWindow::chooseBackgroundColor);
    // connect tff editor
    connect(ui->transferFunctionEditor->getEditor(), &TransferFunctionEditor::gradientStopsChanged,
            ui->volumeRenderWidget, &VolumeRenderWidget::updateTransferFunction);
    connect(ui->pbResetTff, &QPushButton::clicked,
            ui->transferFunctionEditor, &TransferFunctionWidget::resetTransferFunction);
    connect(ui->cbInterpolation, qOverload<const QString &>(&QComboBox::currentIndexChanged),
            ui->volumeRenderWidget, &VolumeRenderWidget::setTffInterpolation);
    connect(ui->cbInterpolation, SIGNAL(currentIndexChanged(QString)),
            ui->transferFunctionEditor, SLOT(setInterpolation(QString)));
//    connect(ui->cbInterpolation, qOverload<const QString &>(&QComboBox::currentIndexChanged),
//            ui->transferFunctionEditor, &TransferFunctionEditor::setInterpolation);

    _statusLabel = new QLabel("No data loaded yet.");
    ui->statusBar->addPermanentWidget(_statusLabel);

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


/**
 * @brief MainWindow::writeSettings
 */
void MainWindow::writeSettings()
{
    _settings->beginGroup("MainWindow");
    _settings->setValue("size", size());
    _settings->setValue("pos", pos());
    _settings->endGroup();

    _settings->beginGroup("Settings");
    // todo
    _settings->endGroup();
}


/**
 * @brief MainWindow::readSettings
 */
void MainWindow::readSettings()
{
    _settings->beginGroup("MainWindow");
    resize(_settings->value("size", QSize(400, 400)).toSize());
    move(_settings->value("pos", QPoint(200, 200)).toPoint());
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
}


/**
 * @brief MainWindow::readVolumeFile
 * @param fileName
 * @return
 */
bool MainWindow::readVolumeFile(const QString &fileName)
{
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
 * @brief MainWindow::finishedLoading
 */
void MainWindow::finishedLoading()
{
    _progBar.setValue(100);
    _progBar.hide();
    _timer.stop();
    qDebug() << "finished.";

    QString status;
    status += _fileName;
    status += " - Resolution: ";
    status += QString::number(ui->volumeRenderWidget->getVolumeResolution().x());
    status += " x ";
    status += QString::number(ui->volumeRenderWidget->getVolumeResolution().y());
    status += " x ";
    status += QString::number(ui->volumeRenderWidget->getVolumeResolution().z());
    _statusLabel->setText(status);

//    const QVector3D volRes = ui->volumeRenderWidget->getVolumeResolution();

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
