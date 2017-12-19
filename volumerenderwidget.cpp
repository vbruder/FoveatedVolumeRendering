#include "volumerenderwidget.h"

#include <QPainter>
#include <QGradient>
#include <QImage>
#include <QCoreApplication>
#include <QScreen>

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
    , _imgCount(0)
    , _imgSamplingRate(1)
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

    _volumerender.initialize();
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
        _volumerender.runRaycast(floor(this->size().width() * _imgSamplingRate),
                                 floor(this->size().height()* _imgSamplingRate), _timestep);
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
            QString number = QString("%1").arg(_imgCount++, 5, 10, QChar('0'));
            img.save("frame_" + number + ".png");
            _writeImage = false;
            qInfo() << "Frame written to" << "frame_" + number + ".png";
        }
    }
    p.endNativePainting();

    // render overlays
    paintFPS(p, fps, _volumerender.getLastExecTime());
    paintOrientationAxis(p);

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
    if (w < h)
        qWarning() << "A frame width smaller than height might negatively impact performance!";

    _screenQuadProjMX.setToIdentity();
    _screenQuadProjMX.perspective(53.14f, 1.0f, Z_NEAR, Z_FAR);

    _overlayProjMX.setToIdentity();
    _overlayProjMX.perspective(53.14f, qreal(w)/qreal(h ? h : 1), Z_NEAR, Z_FAR);

    generateOutputTextures(floor(w*_imgSamplingRate), floor(h*_imgSamplingRate));
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

    _volumerender.updateOutputImg(static_cast<size_t>(width), static_cast<size_t>(height), _outTexId);
    updateView(0, 0);
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
    int timesteps = _volumerender.loadVolumeData(fileName.toStdString());
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
 * @brief VolumeRenderWidget::updateTransferFunction
 * @param stops
 */
void VolumeRenderWidget::updateTransferFunction(QGradientStops stops)
{
    const int tffSize = 256;
    const qreal granularity = 4096.0;
    std::vector<uchar> tff(tffSize*4);
    std::vector<ushort> prefixSum(tffSize);

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
    _volumerender.setTransferFunction(tff);

    // TODO: replace with std::exclusicve_scan(std::par, ... )
    std::partial_sum(prefixSum.begin(), prefixSum.end(), prefixSum.begin());
    _volumerender.setTffPrefixSum(prefixSum);

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
    _volumerender.updateView(viewArray);
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
    double t = 800.0;
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
        _overlayProjMX.ortho(QRect(0, 0, width(), height())); // TODO: fix
    else
        _overlayProjMX.perspective(53.14f, qreal(width())/qreal(height() ? height() : 1), Z_NEAR, Z_FAR);
    this->updateView();
}


/**
 * @brief VolumeRenderWidget::resetCam
 */
void VolumeRenderWidget::resetCam()
{
    _rotQuat = QQuaternion(1, 0, 0, 0);
    QVector3D rotAxis = QVector3D(0.f, 0.f, 1.0f).normalized();
    _rotQuat = _rotQuat * QQuaternion::fromAxisAndAngle(rotAxis, -180.0);
    rotAxis = QVector3D(1.f, 0.f, 0.0f).normalized();
    _rotQuat = _rotQuat * QQuaternion::fromAxisAndAngle(rotAxis, -180.0);

    _translation = QVector3D(0, 0, 2.0);
    this->updateView();
}


/**
 * @brief VolumeRenderWidget::setIllumination
 * @param illum
 */
void VolumeRenderWidget::setIllumination(bool illum)
{
    _volumerender.setIllumination(illum);
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
