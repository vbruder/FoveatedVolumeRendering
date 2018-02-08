/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the demonstration applications of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef TRANSFERFUNCTIONWIDGET_H
#define TRANSFERFUNCTIONWIDGET_H

#include <QPointer>
#include <QGradientStop>
#include <QPainter>
#include <QWidget>

class HoverPoints;

/**
 * @brief The ShadeWidget class
 */
class ShadeWidget : public QWidget
{
    Q_OBJECT

public:
    enum ShadeType
    {
        RedShade,
        GreenShade,
        BlueShade,
        ARGBShade
    };

    ShadeWidget(ShadeType type, QWidget *parent);

    void setGradientStops(const QGradientStops &stops);

    void paintEvent(QPaintEvent *e) override;

    QSize sizeHint() const override { return QSize(150, 40); }
    QPolygonF points() const;
    QVector<QColor> colors() const;

    QSharedPointer<HoverPoints> hoverPoints() const { return _pHoverPoints; }

    uint colorAt(int x);

signals:
    void colorsChanged();
    void selectedPointChanged(const QColor color);

private:
    void generateShade();

    ShadeType _pShadeType;
    QImage _pShade;
    QSharedPointer<HoverPoints> _pHoverPoints;
    QLinearGradient _pAlphaGradient;
};


/**
 * @brief The TransferFunctionEditor class
 */
class TransferFunctionEditor : public QWidget
{
    Q_OBJECT
public:
    explicit TransferFunctionEditor(QWidget *parent = 0);

    void setGradientStops(const QGradientStops &stops);

    const QGradientStops getGradientStops() const;

    void resetPoints();

    void setInterpolation(const QString method);

    void setColorSelected(const QColor color);

public slots:
    void pointsUpdated();
    void selectedPointUpdated(const QColor color);

signals:
    void gradientStopsChanged(const QGradientStops &stops);
    void selectedPointChanged(const QColor color);

private:
    ShadeWidget *_pAlphaShade;

    QVector<ShadeWidget *> _shades;
    QGradientStops _stops;
};


/**
 * @brief The TransferFunctionWidget class
 */
class TransferFunctionWidget : public QWidget
{
    Q_OBJECT

public:
    TransferFunctionWidget(QWidget *parent);
    ~TransferFunctionWidget(){}

    QPointer<TransferFunctionEditor> getEditor();

public slots:
    void resetTransferFunction();

    void setInterpolation(QString method);

    void setColorSelected(const QColor color);

private:
    QPointer<TransferFunctionEditor> _pEditor;
};

#endif // TRANSFERFUNCTIONWIDGET_H
