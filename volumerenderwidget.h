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

#include <QObject>
#include <QWidget>
#include <QDir>
#include <QPointer>
#include <QQuaternion>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShader>
#include <QOpenGLShaderProgram>
#include <QMatrix4x4>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QOpenGLTexture>
#include <QPropertyAnimation>
#include <qopenglfunctions_4_3_core.h>
#include <QPainter>

#include "volumerendercl.h"

class VolumeRenderWidget : public QOpenGLWidget, protected QOpenGLFunctions_4_3_Core
{
    Q_OBJECT

public:
    explicit VolumeRenderWidget(QWidget *parent = 0);
    virtual ~VolumeRenderWidget();

    void setupVertexAttribs();

    void setVolumeData(const QString &fileName);

    bool hasData();

    const QVector3D getVolumeResolution();

    void mousePressEvent(QMouseEvent *event) Q_DECL_OVERRIDE;

    void mouseReleaseEvent(QMouseEvent *event) Q_DECL_OVERRIDE;

    void mouseMoveEvent(QMouseEvent *event) Q_DECL_OVERRIDE;

    void wheelEvent(QWheelEvent *event) Q_DECL_OVERRIDE;

    void mouseDoubleClickEvent(QMouseEvent *event) Q_DECL_OVERRIDE;

    void keyReleaseEvent(QKeyEvent *event) Q_DECL_OVERRIDE;

    void updateView(float dx = 0, float dy = 0);


    bool getLoadingFinished() const;
    void setLoadingFinished(bool loadingFinished);

    QVector3D getCamTranslation() const;
    void setCamTranslation(const QVector3D &translation);

    QQuaternion getCamRotation() const;
    void setCamRotation(const QQuaternion &rotQuat);


public slots:
    void cleanup();
    void resetCam();

    void updateSamplingRate(double samplingRate);
    void updateTransferFunction(QGradientStops stops);
    std::vector<unsigned char> getRawTransferFunction(QGradientStops stops) const;
    void setRawTransferFunction(std::vector<unsigned char> tff);

#undef Bool
    void setTffInterpolation(const QString method);
    void setCamOrtho(bool camOrtho);
    void setIllumination(int illum);
    void setLinearInterpolation(bool linear);
    void setContours(bool contours);
    void setAerial(bool aerial);
    void setDrawBox(bool box);
    void setBackgroundColor(const QColor col);
    void setImageSamplingRate(const double samplingRate);
    void setShowOverlay(bool showOverlay);

    void saveFrame();
    void toggleVideoRecording();
    void setTimeStep(int timestep);
    void setAmbientOcclusion(bool ao);

    void generateLowResVolume();

    void read(const QJsonObject &json);
    void write(QJsonObject &json) const;

    void showSelectOpenCL();
signals:
    void fpsChanged(double);
    void frameSizeChanged(QSize);
    void timeSeriesLoaded(int);

protected:
    // Qt specific QOpenGLWidget methods
    void initializeGL() Q_DECL_OVERRIDE;
    void paintGL() Q_DECL_OVERRIDE;
    void resizeGL(int w, int h) Q_DECL_OVERRIDE;

private:
    void paintOrientationAxis(QPainter &p);
    void paintFPS(QPainter &p, const double fps, const double lastTime);
    double calcFPS();

    void generateOutputTextures(int width, int height);

    // -------Members--------
    //
    // OpenGL
    QOpenGLVertexArrayObject _screenQuadVao;
    QOpenGLShaderProgram _spScreenQuad;
    QOpenGLShaderProgram _spOverlaysGL;
    QOpenGLBuffer _quadVbo;
    GLuint _overlayFboId;
    GLuint _overlayTexId;

    QMatrix4x4 _screenQuadProjMX;
    QMatrix4x4 _viewMX;
    QMatrix4x4 _modelMX;
    QMatrix4x4 _coordViewMX;
    QMatrix4x4 _overlayProjMX;
    QMatrix4x4 _overlayModelMX;

    QPoint _tffRange;
    GLuint _outTexId;
    VolumeRenderCL _volumerender;
    QEasingCurve _tffInterpol;
    int _timestep;

    // global rendering flags
    QPoint _lastLocalCursorPos;
    QQuaternion _rotQuat;
    QVector3D _translation;

    bool _noUpdate;
    bool _loadingFinished;
    bool _writeImage;
    bool _recordVideo;
    qint64 _imgCount;
    QVector<double> _times;
    double _imgSamplingRate;       // image oversampling rate
    bool _useGL;
    bool _showOverlay;
};
