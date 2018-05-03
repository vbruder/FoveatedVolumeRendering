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

#include "volumerenderwidget.h"

#include <QPainter>
#include <QGradient>
#include <QImage>
#include <QCoreApplication>
#include <QScreen>
#include <QInputDialog>
#include <QJsonObject>
#include <QErrorMessage>
#include <QLoggingCategory>
#include <QMessageBox>

const static double Z_NEAR = 1.0;
const static double Z_FAR = 500.0;

static const char *pVsScreenQuadSource =
    "#version 330\n"
    "layout(location = 0) in highp vec3 vertex;\n"
    "out highp vec2 texCoord;\n"
    "uniform mat4 projMatrix;\n"
    "uniform mat4 mvMatrix;\n"
    "void main() {\n"
    "   texCoord = vec2(0.5f) + 0.5f * vertex.xy;\n"
    "   gl_Position = projMatrix * mvMatrix * vec4(vertex.xy, 1.0f, 1.0f);\n"
    "}\n";

static const char *pFsScreenQuadSource =
    "#version 330\n"
    "in highp vec2 texCoord;\n"
    "out highp vec4 fragColor;\n"
    "uniform highp sampler2D outTex;\n"
    "void main() {\n"
    "   fragColor = texture(outTex, texCoord);\n"
    "   fragColor.a = 1.0f;\n"
    "}\n";


/**
 * @brief VolumeRenderWidget::VolumeRenderWidget
 * @param parent
 */
VolumeRenderWidget::VolumeRenderWidget(QWidget *parent)
    : QOpenGLWidget(parent)
    , _tffRange(QPoint(0, 255))
    , _timestep(0)
    , _lastLocalCursorPos(QPoint(0,0))
    , _rotQuat(QQuaternion(1, 0, 0, 0))
    , _translation(QVector3D(0, 0, 2.0))
    , _noUpdate(true)
    , _loadingFinished(false)
    , _writeImage(false)
    , _recordVideo(false)
    , _imgCount(0)
    , _imgSamplingRate(1)
    , _useGL(true)
    , _showOverlay(true)
{
    this->setMouseTracking(true);
}


/**
 * @brief VolumeRenderWidget::~VolumeRenderWidget
 */
VolumeRenderWidget::~VolumeRenderWidget()
{
}


/**
 * @brief VolumeRenderWidget::paintOrientationAxis
 */
void VolumeRenderWidget::paintOrientationAxis(QPainter &p)
{
    QMatrix4x4 proj;
    proj.perspective(53.14f, 1.0f, 0.1f, 1.0);
    QMatrix4x4 viewProj(proj * _coordViewMX);
    QVector4D x =         viewProj * QVector4D( 20,  0,  0, 0);
    QVector4D xArrLeft =  viewProj * QVector4D( 16, -2,  0, 0);
    QVector4D xArrRight = viewProj * QVector4D( 16, +2,  0, 0);
    QVector4D y =         viewProj * QVector4D(  0, 20,  0, 0);
    QVector4D yArrLeft =  viewProj * QVector4D( -2, 16,  0, 0);
    QVector4D yArrRight = viewProj * QVector4D( +2, 16,  0, 0);
    QVector4D z =         viewProj * QVector4D(  0,  0, 20, 0);
    QVector4D zArrLeft =  viewProj * QVector4D( -2,  0, 16, 0);
    QVector4D zArrRight = viewProj * QVector4D( +2,  0, 16, 0);

    p.resetTransform();
    p.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
    p.translate(66, height() - 66);
    int textOffset = 5;
    // x axis
    p.setPen(Qt::red);
    p.drawLine(0, 0, x.x(), x.y());
    p.drawLine(xArrLeft.x(), xArrLeft.y(), x.x(), x.y());
    p.drawLine(xArrRight.x(), xArrRight.y(), x.x(), x.y());
    p.drawText(x.x() + textOffset, x.y() + textOffset, "x");
    // y axis
    p.setPen(Qt::green);
    p.drawLine(0, 0, y.x(), y.y());
    p.drawLine(yArrLeft.x(), yArrLeft.y(), y.x(), y.y());
    p.drawLine(yArrRight.x(), yArrRight.y(), y.x(), y.y());
    p.drawText(y.x() + textOffset, y.y() + textOffset, "y");
    // z axis
    p.setPen(Qt::blue);
    p.drawLine(0, 0, z.x(), z.y());
    p.drawLine(zArrLeft.x(), zArrLeft.y(), z.x(), z.y());
    p.drawLine(zArrRight.x(), zArrRight.y(), z.x(), z.y());
    p.drawText(z.x() + textOffset, z.y() + textOffset, "z");
}


/**
 * @brief paintFPS
 */
void VolumeRenderWidget::paintFPS(QPainter &p, const double fps, const double lastTime)
{
    p.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
    p.setPen(Qt::darkGreen);
    p.setFont(QFont("Helvetica", 11));
    QString s = "FPS: " + QString::number(fps);
    p.drawText(10, 20, s);
    s = "Last: " + QString::number(lastTime);
    p.drawText(10, 36, s);
}


/**
 * @brief VolumeRenderWidget::initializeGL
 */
void VolumeRenderWidget::initializeGL()
{
    connect( context(), &QOpenGLContext::aboutToBeDestroyed, this,
             &VolumeRenderWidget::cleanup );

    initializeOpenGLFunctions();
    makeCurrent();

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    _spScreenQuad.addShaderFromSourceCode(QOpenGLShader::Vertex, pVsScreenQuadSource);
    _spScreenQuad.addShaderFromSourceCode(QOpenGLShader::Fragment, pFsScreenQuadSource);
    _spScreenQuad.bindAttributeLocation("vertex", 0);
    _spScreenQuad.link();

    _spScreenQuad.bind();
    _screenQuadVao.create();
    _screenQuadVao.bind();

    const int numQuadVertices = 8;
    GLfloat quadVertices[numQuadVertices] =
    {
        -1.0f, -1.0f,
        +1.0f, -1.0f,
        -1.0f, +1.0f,
        +1.0f, +1.0f
    };

    // Setup vertex buffer object.
    _quadVbo.create();
    _quadVbo.bind();
    _quadVbo.allocate( quadVertices, numQuadVertices * sizeof( GLfloat ) );
    _quadVbo.release();
    // Store the vertex attribute bindings for the program.
    setupVertexAttribs();
    _viewMX.setToIdentity();
    _viewMX.translate( 0.0f, 0.0f, -1.0f );
    // set quad model matrix
    _modelMX.setToIdentity();
    _modelMX.rotate( 180.0f, 1, 0, 0 );
    _spScreenQuad.release();
    _screenQuadVao.release();

    try
    {
        _volumerender.initialize(true, false);
    }
    catch (std::invalid_argument e)
    {
        qCritical() << e.what();
    }
    catch (std::runtime_error e)
    {
        qCritical() << e.what() << "\nSwitching to CPU fallback mode.";
        _useGL = false;
        _volumerender.initialize(false, true);
    }
    catch (...)
    {
        qCritical() << "An unknown error occured initializing OpenCL/OpenGL.";
    }
}


/**
 * @brief VolumeRenderWidget::saveFrame
 */
void VolumeRenderWidget::saveFrame()
{
    _writeImage = true;
    update();
}

/**
 * @brief VolumeRenderWidget::toggleVideoRecording
 */
void VolumeRenderWidget::toggleVideoRecording()
{
    qInfo() << (_recordVideo ? "Stopped recording." : "Started recording.");

    _recordVideo = !_recordVideo;
    _writeImage = true;
    update();
}

/**
 * @brief VolumeRenderWidget::setTimeStep
 * @param timestep
 */
void VolumeRenderWidget::setTimeStep(int timestep)
{
    _timestep = timestep;
    update();
}

/**
 * @brief VolumeRenderWidget::setImageSamplingRate
 * @param samplingRate
 */
void VolumeRenderWidget::setImageSamplingRate(const double samplingRate)
{
    _imgSamplingRate = samplingRate;
    this->resizeGL(this->width(), this->height());
}

/**
 * @brief VolumeRenderWidget::paintGL
 */
void VolumeRenderWidget::paintGL()
{
    double fps = 0.0;
    if (this->_loadingFinished && _volumerender.hasData() && !_noUpdate)
    {
        // OpenCL raycast
        try
        {
            if (_useGL)
                _volumerender.runRaycast(floor(this->size().width() * _imgSamplingRate),
                                         floor(this->size().height()* _imgSamplingRate), _timestep);
            else
            {
                std::vector<float> d;
                _volumerender.runRaycastNoGL(floor(this->size().width() * _imgSamplingRate),
                                             floor(this->size().height()* _imgSamplingRate),
                                             _timestep, d);
                glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F,
                             floor(this->size().width() * _imgSamplingRate),
                             floor(this->size().height()* _imgSamplingRate),
                             0, GL_RGBA, GL_FLOAT,
                             d.data());
                glGenerateMipmap(GL_TEXTURE_2D);
                _volumerender.updateOutputImg(static_cast<size_t>(width()),
                                              static_cast<size_t>(height()), _outTexId);
            }
        }
        catch (std::runtime_error e)
        {
            qCritical() << e.what();
        }
        fps = calcFPS();
    }

    QPainter p(this);
    p.beginNativePainting();
    {
        // render the ray casting output
        // clear to white to avoid getting colored borders outside the quad
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

        // draw screen quad
        //
        _screenQuadVao.bind();
        _quadVbo.bind();
        glVertexAttribPointer( 0, 2, GL_FLOAT, GL_FALSE, 0, 0 );

        // render screen quad
        //
        _spScreenQuad.bind();
        _spScreenQuad.setUniformValue( _spScreenQuad.uniformLocation( "projMatrix" ),
                                          _screenQuadProjMX );
        _spScreenQuad.setUniformValue( _spScreenQuad.uniformLocation( "mvMatrix" ),
                                          _viewMX * _modelMX );

        _spScreenQuad.setUniformValue(_spScreenQuad.uniformLocation("outTex"), GL_TEXTURE0);
        glDrawArrays( GL_TRIANGLE_STRIP, 0, 4 );
        _screenQuadVao.release();
        _quadVbo.release();
        _spScreenQuad.release();

        glDisable(GL_CULL_FACE);
        glDisable(GL_DEPTH_TEST);

        if (_volumerender.hasData() && _writeImage)
        {
            QImage img = this->grabFramebuffer();
            QString number = QString("%1").arg(_imgCount++, 6, 10, QChar('0'));
            if (!_recordVideo)
            {
                QLoggingCategory category("screenshot");
                qCInfo(category, "Writing current frame img/frame_%s.png", number.toStdString().c_str());
                _writeImage = false;
            }
            if (!QDir("img").exists())
                QDir().mkdir("img");
            img.save("img/frame_" + number + ".png");
        }
    }
    p.endNativePainting();

    // render overlays
    if (_showOverlay)
    {
        paintFPS(p, fps, _volumerender.getLastExecTime());
        paintOrientationAxis(p);
    }

    // recover opengl texture
    p.beginNativePainting();
    {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, _outTexId);
    }
    p.endNativePainting();
    p.end();
}


/**
 * @brief VolumeRenderWidget::resizeGL
 * @param w
 * @param h
 */
void VolumeRenderWidget::resizeGL(int w, int h)
{
    _screenQuadProjMX.setToIdentity();
    _screenQuadProjMX.perspective(53.14f, 1.0f, Z_NEAR, Z_FAR);

    _overlayProjMX.setToIdentity();
    _overlayProjMX.perspective(53.14f, qreal(w)/qreal(h ? h : 1), Z_NEAR, Z_FAR);

    try
    {
        generateOutputTextures(floor(w*_imgSamplingRate), floor(h*_imgSamplingRate));
    }
    catch (std::runtime_error e)
    {
        qCritical() << "An error occured while generating output texture." << e.what();
    }

    emit frameSizeChanged(this->size());
}


/**
 * @brief VolumeRenderWidget::generateOutputTextures
 */
void VolumeRenderWidget::generateOutputTextures(int width, int height)
{
    glGenTextures(1, &_outTexId);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, _outTexId);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    if (!_volumerender.hasData())
    {
        QImage img(width, height, QImage::Format_RGBA8888);
        img.fill(Qt::white);
        QPainter p(&img);
        p.setFont(QFont("Helvetia", 12));
        p.drawText(width/2 - 110, height/2, "Drop your volume data file here.");
        p.end();
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8,
                     width, height, 0,
                     GL_RGBA, GL_UNSIGNED_BYTE,
                     img.bits());
    }
    else
    {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8,
                     width, height, 0,
                     GL_RGBA, GL_UNSIGNED_BYTE,
                     NULL);
    }
    glGenerateMipmap(GL_TEXTURE_2D);

    _volumerender.updateOutputImg(static_cast<size_t>(width), static_cast<size_t>(height),
                                      _outTexId);

    updateView(0, 0);
}

void VolumeRenderWidget::setShowOverlay(bool showOverlay)
{
    _showOverlay = showOverlay;
    updateView();
}

QQuaternion VolumeRenderWidget::getCamRotation() const
{
    return _rotQuat;
}

void VolumeRenderWidget::setCamRotation(const QQuaternion &rotQuat)
{
    _rotQuat = rotQuat;
}

QVector3D VolumeRenderWidget::getCamTranslation() const
{
    return _translation;
}

void VolumeRenderWidget::setCamTranslation(const QVector3D &translation)
{
    _translation = translation;
}

/**
 * @brief VolumeRenderWidget::showSelectOpenCL
 */
void VolumeRenderWidget::showSelectOpenCL()
{
    std::vector<std::string> names;
    try
    {
        names = _volumerender.getPlatformNames();
    }
    catch (std::runtime_error e)
    {
        qCritical() << "An error occured while trying to retrieve OpenCL platform information."
                    << e.what();
        return;
    }

    QStringList platforms;
    for (std::string &s : names)
        platforms.append(QString::fromStdString(s));

    bool ok;
    QString platform = QInputDialog::getItem(this, tr("Select platform"),
                                             tr("Select OpenCL platform:"),
                                             platforms, 0, false, &ok);
    if (ok && !platform.isEmpty())
    {
        cl_vendor vendor = VENDOR_ANY;
        QString type = "GPU";
        _useGL = false;

        // filter GPU only platforms
        if (platform.contains("NVIDIA"))
            vendor = VENDOR_NVIDIA;
        else
        {
            if (!platform.contains("Graphics"))
            {
                QStringList types;
                types << "GPU" << "CPU";
                type = QInputDialog::getItem(this, tr("Select type"),
                                             tr("Select device type:"), types, 0, false, &ok);
            }
            if (platform.contains("Advanced Micro Devices", Qt::CaseInsensitive))
                vendor = VENDOR_AMD;
            else if (platform.contains("Intel", Qt::CaseInsensitive))
                vendor = VENDOR_INTEL;

        }
        if (!type.isEmpty())
        {
            if (type == "GPU")
            {
                QMessageBox msgBox;
                msgBox.setText("Do you wish to try OpenGL context sharing using this platform?");
                msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
                msgBox.setDefaultButton(QMessageBox::Yes);
                int ret = msgBox.exec();
                _useGL = ret == QMessageBox::Yes;
            }

            int platformId = platforms.indexOf(platform);
            try
            {
                names = _volumerender.getDeviceNames(platformId, type.toStdString());
            }
            catch (std::runtime_error e)
            {
                qCritical() << "No capable device found on using the selected platform and type."
                            << e.what();
                return;
            }
            QStringList devices;
            for (std::string &s : names)
                devices.append(QString::fromStdString(s));

            QString device;
            if (devices.empty())
                throw std::runtime_error("No capable device found on the selected platform.");
            else if (devices.size() == 1)
                device = devices.front();
            else
                device = QInputDialog::getItem(this, tr("Select device"),
                                                   tr("Select OpenCL device:"),
                                                   devices, 0, false, &ok);
            if (!device.isEmpty())
            {
                try {
                    _volumerender.initialize(_useGL, type == "CPU", vendor);
                } catch (std::runtime_error e) {
                    qCritical() << e.what() << "\nSwitching to CPU fallback mode.";
                    _useGL = false;
                    _volumerender.initialize(false, true);
                }
                catch (...) {
                    qCritical() << "An unknown error occured initializing OpenCL/OpenGL.";
                }

                try {
                    generateOutputTextures(floor(width()*_imgSamplingRate),
                                           floor(height()*_imgSamplingRate));
                } catch (std::runtime_error e) {
                    qCritical() << "An error occured while generating output texture." << e.what();
                }
            }
        }
    }
}


/**
 * @brief VolumeRenderWidget::setupVertexAttribs
 */
void VolumeRenderWidget::setupVertexAttribs()
{
    // screen quad
    _quadVbo.bind();
    glEnableVertexAttribArray( 0 );
    glVertexAttribPointer( 0, 2, GL_FLOAT, GL_FALSE, 0, 0 );
    _quadVbo.release();
}


/**
 * @brief VolumeRenderWidget::setVolumeData
 * @param fileName
 */
void VolumeRenderWidget::setVolumeData(const QString &fileName)
{
    this->_noUpdate = true;
    int timesteps = 0;
    try
    {
        timesteps = _volumerender.loadVolumeData(fileName.toStdString());
    }
    catch (std::invalid_argument e)
    {
        qCritical() << e.what();
    }
    catch (std::runtime_error e)
    {
        qCritical() << e.what();
    }
    emit timeSeriesLoaded(timesteps - 1);

    _overlayModelMX.setToIdentity();
    QVector3D res = getVolumeResolution();
    _overlayModelMX.scale(res / qMax(res.x(), qMax(res.y(), res.z())));
    this->_noUpdate = false;
    update();
}


/**
 * @brief VolumeRenderWidget::hasData
 * @return
 */
bool VolumeRenderWidget::hasData()
{
    return _volumerender.hasData();
}


/**
 * @brief VolumeRenderWidget::getVolumeResolution
 * @return
 */
const QVector3D VolumeRenderWidget::getVolumeResolution()
{
    if (_volumerender.hasData() == false)
        return QVector3D();

    return QVector3D(_volumerender.getResolution().at(0),
                     _volumerender.getResolution().at(1),
                     _volumerender.getResolution().at(2));
}


/**
 * @brief VolumeRenderWidget::updateStepSize
 * @param stepSize
 */
void VolumeRenderWidget::updateSamplingRate(double samplingRate)
{
    _volumerender.updateSamplingRate(samplingRate);
    update();
}


/**
 * @brief VolumeRenderWidget::setInterpolation
 * @param method
 */
void VolumeRenderWidget::setTffInterpolation(const QString method)
{
    if (method.contains("Quad"))
        _tffInterpol = QEasingCurve::InOutQuad;
    else if (method.contains("Linear"))
        _tffInterpol = QEasingCurve::Linear;
}

/**
 * @brief VolumeRenderWidget::setRawTransferFunction
 * @param tff
 */
void VolumeRenderWidget::setRawTransferFunction(std::vector<unsigned char> tff)
{
    try
    {
        _volumerender.setTransferFunction(tff);
        update();
    }
    catch (std::runtime_error e)
    {
        qCritical() << e.what();
    }
}

/**
 * @brief VolumeRenderWidget::updateTransferFunction
 * @param stops
 */
void VolumeRenderWidget::updateTransferFunction(QGradientStops stops)
{
    const int tffSize = 256;
    const qreal granularity = 4096.0;
    std::vector<uchar> tff(tffSize*4);
    std::vector<unsigned int> prefixSum(tffSize);

    QPropertyAnimation interpolator;
    interpolator.setEasingCurve(_tffInterpol);
    interpolator.setDuration(granularity);
    foreach (QGradientStop stop, stops)
    {
        interpolator.setKeyValueAt(stop.first, stop.second);
    }
//    tff.at(0) = (uchar)0;
//    tff.at(1) = (uchar)0;
//    tff.at(2) = (uchar)0;
//    tff.at(3) = (uchar)0;
#pragma omp for
    for (int i = 0; i < tffSize; ++i)
    {
        interpolator.setCurrentTime((i/static_cast<double>(tffSize)) * granularity);
        tff.at(i*4 + 0) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().red()   - 3);
        tff.at(i*4 + 1) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().green() - 3);
        tff.at(i*4 + 2) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().blue()  - 3);
        tff.at(i*4 + 3) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().alpha() - 3);
        prefixSum.at(i) = tff.at(i*4 + 3);
    }
    try
    {
        _volumerender.setTransferFunction(tff);
        // TODO: replace with std::exclusicve_scan(std::par, ... )
        std::partial_sum(prefixSum.begin(), prefixSum.end(), prefixSum.begin());
        _volumerender.setTffPrefixSum(prefixSum);
    }
    catch (std::runtime_error e)
    {
        qCritical() << e.what();
    }

    update();
}

std::vector<unsigned char> VolumeRenderWidget::getRawTransferFunction(QGradientStops stops) const
{
    const size_t tffSize = 256;
    const qreal granularity = 4096.0;
    std::vector<uchar> tff(tffSize*4);

    QPropertyAnimation interpolator;
    interpolator.setEasingCurve(_tffInterpol);
    interpolator.setDuration(granularity);
    foreach (QGradientStop stop, stops)
    {
        interpolator.setKeyValueAt(stop.first, stop.second);
    }

    for (size_t i = 0; i < tffSize; ++i)
    {
        interpolator.setCurrentTime((i/static_cast<double>(tffSize)) * granularity);
        tff.at(i*4 + 0) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().red()   - 3);
        tff.at(i*4 + 1) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().green() - 3);
        tff.at(i*4 + 2) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().blue()  - 3);
        tff.at(i*4 + 3) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().alpha() - 3);
    }
    return tff;
}


/**
 * @brief VolumeRenderWidget::cleanup
 */
void VolumeRenderWidget::cleanup()
{
//    makeCurrent();
//    if (_quadVbo.isCreated())
//        _quadVbo.destroy();
}



/**
 * @brief VolumeRenderWidget::mousePressEvent
 * @param event
 */
void VolumeRenderWidget::mousePressEvent(QMouseEvent *event)
{
    // nothing yet
    _lastLocalCursorPos = event->pos();
}


/**
 * @brief VolumeRenderWidget::mouseReleaseEvent
 * @param event
 */
void VolumeRenderWidget::mouseReleaseEvent(QMouseEvent *event)
{
    event->accept();
}

/**
 * @brief VolumeRenderWidget::resetCam
 */
void VolumeRenderWidget::resetCam()
{
    _rotQuat = QQuaternion(1, 0, 0, 0);
    _translation = QVector3D(0, 0, 2.0);
    updateView();
}

/**
 * @brief update camera view
 * @param dx
 * @param dy
 */
void VolumeRenderWidget::updateView(float dx, float dy)
{
    QVector3D rotAxis = QVector3D(dy, dx, 0.0f).normalized();
    double angle = QVector2D(dx, dy).length()*500;
    _rotQuat = _rotQuat * QQuaternion::fromAxisAndAngle(rotAxis, -angle);

    QMatrix4x4 viewMat;
//    viewMat.scale(QVector3D(1,1,1./_zScale));
    viewMat.rotate(_rotQuat);

    _coordViewMX.setToIdentity();
    _coordViewMX.scale(1, -1, 1);
    _coordViewMX.translate(_translation * -1.0);
    _coordViewMX *= QMatrix4x4(_rotQuat.toRotationMatrix().transposed());

    viewMat.translate(_translation);
    viewMat.scale(_translation.z());

    std::array<float, 16> viewArray;
    for (size_t i = 0; i < viewArray.size(); ++i)
    {
        viewArray.at(i) = viewMat.transposed().constData()[i];
    }

    try
    {
        _volumerender.updateView(viewArray);
    }
    catch (std::runtime_error e)
    {
        qCritical() << e.what();
    }
    update();
}


/**
 * @brief VolumeRenderWidget::mouseMoveEvent
 * @param event
 */
void VolumeRenderWidget::mouseMoveEvent(QMouseEvent *event)
{
    float dx = (float)(event->pos().x() - _lastLocalCursorPos.x()) / width();
    float dy = (float)(event->pos().y() - _lastLocalCursorPos.y()) / height();

    // rotate object
    if (event->buttons() & Qt::LeftButton)
    {
        if (event->modifiers() & Qt::ShiftModifier)
        {
            dx *= 0.1f;
            dy *= 0.1f;
        }
        updateView(dx, dy);
    }
    // translate object
    if (event->buttons() & Qt::MiddleButton)
    {
        int sensitivity = 6;
        if (event->modifiers() & Qt::ShiftModifier)
            sensitivity = 1;

        _translation.setX(_translation.x() - dx*sensitivity);
        _translation.setY(_translation.y() + dy*sensitivity);
        updateView();
    }

    _lastLocalCursorPos = event->pos();
    event->accept();
}


/**
 * @brief VolumeRenderWidget::wheelEvent
 * @param event
 */
void VolumeRenderWidget::wheelEvent(QWheelEvent *event)
{
    double t = 1600.0;
    if (event->modifiers() & Qt::ShiftModifier)
        t *= 6.0;

    // limit translation to origin, otherwise camera setup breaks (flips)
    _translation.setZ(qMax(0.01, _translation.z() - event->angleDelta().y() / t));
    updateView();
    event->accept();
}


/**
 * @brief VolumeRenderWidget::mouseDoubleClickEvent
 * @param event
 */
void VolumeRenderWidget::mouseDoubleClickEvent(QMouseEvent *event)
{
    // nothing yet
    event->accept();
}



/**
 * @brief VolumeRenderWidget::keyReleaseEvent
 * @param event
 */
void VolumeRenderWidget::keyReleaseEvent(QKeyEvent *event)
{
    // nothing yet
    event->accept();
}


/**
 * @brief VolumeRenderWidget::getLoadingFinished
 * @return
 */
bool VolumeRenderWidget::getLoadingFinished() const
{
    return _loadingFinished;
}


/**
 * @brief VolumeRenderWidget::setLoadingFinished
 * @param loadingFinished
 */
void VolumeRenderWidget::setLoadingFinished(bool loadingFinished)
{
    this->setTimeStep(0);
    _loadingFinished = loadingFinished;
}


/**
 * @brief VolumeRenderWidget::setCamOrtho
 * @param camOrtho
 */
void VolumeRenderWidget::setCamOrtho(bool camOrtho)
{
    _volumerender.setCamOrtho(camOrtho);
    _overlayProjMX.setToIdentity();
    if (camOrtho)
        _overlayProjMX.ortho(QRect(0, 0, width(), height()));
    else
        _overlayProjMX.perspective(53.14f, qreal(width())/qreal(height() ? height() : 1), Z_NEAR, Z_FAR);
    this->updateView();
}


/**
 * @brief VolumeRenderWidget::setIllumination
 * @param illum
 */
void VolumeRenderWidget::setIllumination(int illum)
{
    _volumerender.setIllumination(illum);
    this->updateView();
}


/**
 * @brief VolumeRenderWidget::setAmbientOcclusion
 * @param ao
 */
void VolumeRenderWidget::setAmbientOcclusion(bool ao)
{
    _volumerender.setAmbientOcclusion(ao);
    this->updateView();
}

/**
 * @brief VolumeRenderWidget::setLinearInterpolation
 * @param linear
 */
void VolumeRenderWidget::setLinearInterpolation(bool linear)
{
    _volumerender.setLinearInterpolation(linear);
    this->updateView();
}

/**
 * @brief VolumeRenderWidget::setContours
 * @param contours
 */
void VolumeRenderWidget::setContours(bool contours)
{
    _volumerender.setContours(contours);
    this->updateView();
}

/**
 * @brief VolumeRenderWidget::setAerial
 * @param aerial
 */
void VolumeRenderWidget::setAerial(bool aerial)
{
    _volumerender.setAerial(aerial);
    this->updateView();
}


/**
 * @brief VolumeRenderWidget::setDrawBox
 * @param box
 */
void VolumeRenderWidget::setDrawBox(bool box)
{
    _volumerender.setBoundingBox(box);
    this->updateView();
}


/**
 * @brief VolumeRenderWidget::setBackgroundColor
 * @param col
 */
void VolumeRenderWidget::setBackgroundColor(const QColor col)
{
    std::array<float, 4> color = {{(float)col.redF(), (float)col.greenF(),
                                   (float)col.blueF(), (float)col.alphaF()}};
    _volumerender.setBackground(color);
    this->updateView();
}


/**
 * @brief VolumeRenderWidget::calcFPS
 * @return
 */
double VolumeRenderWidget::calcFPS()
{
    _times.push_back(_volumerender.getLastExecTime());
    if (_times.length() > 60)
        _times.pop_front();

    double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < _times.size(); ++i)
        sum += _times.at(i);

    double fps = 1.0/(sum/_times.length());

//    qDebug() << fps;
    return fps;
}


/**
 * @brief VolumeRenderWidget::generateLowResVolume
 * @param factor
 */
void VolumeRenderWidget::generateLowResVolume()
{
    bool ok;
    int factor = QInputDialog::getInt(this, tr("Factor"),
                                         tr("Select downsampling factor:"), 2, 2, 64, 1, &ok);
    if (ok)
    {
        try
        {
            std::string name = _volumerender.volumeDownsampling(_timestep, factor);
            QLoggingCategory category("volumeDownSampling");
            qCInfo(category, "Successfully created down-sampled volume data set: '%s.raw'",
                   name.c_str());
        }
        catch (std::invalid_argument e)
        {
            qCritical() << e.what();
        }
        catch (std::runtime_error e)
        {
            qCritical() << e.what();
        }
    }
}

/**
 * @brief VolumeRenderWidget::read
 * @param json
 */
void VolumeRenderWidget::read(const QJsonObject &json)
{
    if (json.contains("camRotation"))
    {
        QStringList sl = json["camRotation"].toVariant().toString().split(' ');
        if (sl.length() >= 4)
        {
            _rotQuat.setScalar(sl.at(0).toFloat());
            _rotQuat.setX(sl.at(1).toFloat());
            _rotQuat.setY(sl.at(2).toFloat());
            _rotQuat.setZ(sl.at(3).toFloat());
        }
    }
    if (json.contains("camTranslation"))
    {
        QStringList sl = json["camTranslation"].toVariant().toString().split(' ');
        if (sl.length() >= 3)
        {
            _translation.setX(sl.at(0).toFloat());
            _translation.setY(sl.at(1).toFloat());
            _translation.setZ(sl.at(2).toFloat());
        }
    }

    updateView();
}

/**
 * @brief VolumeRenderWidget::write
 * @param json
 */
void VolumeRenderWidget::write(QJsonObject &json) const
{
    QString sTmp = QString::number(_rotQuat.scalar()) + " " + QString::number(_rotQuat.x())
                   + " " + QString::number(_rotQuat.y()) + " " + QString::number(_rotQuat.z());
    json["camRotation"] = sTmp;
    sTmp = QString::number(_translation.x()) + " " + QString::number(_translation.y())
            + " " + QString::number(_translation.z());
    json["camTranslation"] = sTmp;
}
