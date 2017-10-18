#include "transferfunctionwidget.h"
#include "inc/painting/hoverpoints.h"

#include <algorithm>

/**
 * @brief ShadeWidget::ShadeWidget
 * @param type
 * @param parent
 */
ShadeWidget::ShadeWidget(ShadeType type, QWidget *parent)
    : QWidget(parent), _pShadeType(type), _pAlphaGradient(QLinearGradient(0, 0, 0, 0))
{
    // Checkers background
    if (_pShadeType == ARGBShade)
    {
        QPixmap pm(20, 20);
        QPainter pmp(&pm);
        pmp.fillRect(0, 0, 10, 10, Qt::white);
        pmp.fillRect(10, 10, 10, 10, Qt::white);
        pmp.fillRect(0, 10, 10, 10, Qt::lightGray);
        pmp.fillRect(10, 0, 10, 10, Qt::lightGray);
        pmp.end();
        QPalette pal = palette();
        pal.setBrush(backgroundRole(), QBrush(pm));
        setAutoFillBackground(true);
        setPalette(pal);
    }
    else
    {
        setAttribute(Qt::WA_NoBackground);
    }

    QPolygonF points;
    points << QPointF(0, sizeHint().height())
           << QPointF(sizeHint().width(), 0);

    _pHoverPoints = new HoverPoints(this, HoverPoints::CircleShape);
    _pHoverPoints->setConnectionType(HoverPoints::LineConnection);
    _pHoverPoints->setPoints(points);
    _pHoverPoints->setPointLock(0, HoverPoints::LockToLeft);
    _pHoverPoints->setPointLock(1, HoverPoints::LockToRight);
    _pHoverPoints->setSortType(HoverPoints::XSort);

    setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);

    connect(_pHoverPoints, SIGNAL(pointsChanged(QPolygonF)), this, SIGNAL(colorsChanged()));
    connect(_pHoverPoints, &HoverPoints::selectionChanged, this, &ShadeWidget::selectedPointChanged);
}

QPolygonF ShadeWidget::points() const
{
    return _pHoverPoints->points();
}

QVector<QColor> ShadeWidget::colors() const
{
    return _pHoverPoints->colors();
}

uint ShadeWidget::colorAt(int x)
{
    generateShade();

    QPolygonF pts = _pHoverPoints->points();
    for (int i=1; i < pts.size(); ++i)
    {
        if (pts.at(i-1).x() <= x && pts.at(i).x() >= x)
        {
            QLineF l(pts.at(i-1), pts.at(i));
            l.setLength(l.length() * ((x - l.x1()) / l.dx()));
            return _pShade.pixel(qRound(qMin(l.x2(), (qreal(_pShade.width() - 1)))),
                                 qRound(qMin(l.y2(), qreal(_pShade.height() - 1))));
        }
    }
    return 0;
}

void ShadeWidget::setGradientStops(const QGradientStops &stops)
{
    if (_pShadeType == ARGBShade)
    {
        _pAlphaGradient = QLinearGradient(0, 0, width(), 0);
        for (int i=0; i<stops.size(); ++i)
        {
            QColor c = stops.at(i).second;
            _pAlphaGradient.setColorAt(stops.at(i).first, c);// QColor(c.red(), c.green(), c.blue()));
        }
        _pShade = QImage();
        generateShade();
        update();
    }
}

void ShadeWidget::paintEvent(QPaintEvent *)
{
    generateShade();

    QPainter p(this);
    p.drawImage(0, 0, _pShade);

    p.setPen(QColor(146, 146, 146));
    p.drawRect(0, 0, width() - 1, height() - 1);
}

void ShadeWidget::generateShade()
{
    if (_pShade.isNull() || _pShade.size() != size())
    {
        if (_pShadeType == ARGBShade)
        {
            _pShade = QImage(size(), QImage::Format_ARGB32_Premultiplied);
            _pShade.fill(0);

            QPainter p(&_pShade);
            p.fillRect(rect(), _pAlphaGradient);

            p.setCompositionMode(QPainter::CompositionMode_DestinationIn);
            QLinearGradient fade(0, 0, 0, height());
            fade.setColorAt(0, QColor(0, 0, 0, 255));
            fade.setColorAt(1, QColor(0, 0, 0, 0));
            p.fillRect(rect(), fade);
        }
        else
        {
            _pShade = QImage(size(), QImage::Format_RGB32);
            QLinearGradient shade(0, 0, 0, height());
            shade.setColorAt(1, Qt::black);

            if (_pShadeType == RedShade)
                shade.setColorAt(0, Qt::red);
            else if (_pShadeType == GreenShade)
                shade.setColorAt(0, Qt::green);
            else
                shade.setColorAt(0, Qt::blue);

            QPainter p(&_pShade);
            p.fillRect(rect(), shade);
        }
    }
}


/**
 * @brief TransferFunctionWidget::TransferFunctionWidget
 * @param parent
 */
TransferFunctionEditor::TransferFunctionEditor(QWidget *parent) : QWidget(parent)
{
    QVBoxLayout *vbox = new QVBoxLayout(this);
    vbox->setSpacing(1);
    vbox->setMargin(1);

//    _pRedShade = new ShadeWidget(ShadeWidget::RedShade, this);
//    _shades.push_back(_pRedShade);
//    _pGreenShade = new ShadeWidget(ShadeWidget::GreenShade, this);
//    _shades.push_back(_pGreenShade);
//    _pBlueShade = new ShadeWidget(ShadeWidget::BlueShade, this);
//    _shades.push_back(_pBlueShade);
    _pAlphaShade = new ShadeWidget(ShadeWidget::ARGBShade, this);
    _shades.push_back(_pAlphaShade);

    foreach (ShadeWidget *s, _shades)
    {
        vbox->addWidget(s);
    }

//    connect(_pRedShade, SIGNAL(colorsChanged()), this, SLOT(pointsUpdated()));
//    connect(_pGreenShade, SIGNAL(colorsChanged()), this, SLOT(pointsUpdated()));
//    connect(_pBlueShade, SIGNAL(colorsChanged()), this, SLOT(pointsUpdated()));
    connect(_pAlphaShade, &ShadeWidget::colorsChanged, this, &TransferFunctionEditor::pointsUpdated);
    connect(_pAlphaShade, &ShadeWidget::selectedPointChanged,
            this, &TransferFunctionEditor::selectedPointUpdated);
}


inline static bool x_less_than(const QPointF &p1, const QPointF &p2)
{
    return p1.x() < p2.x();
}

void TransferFunctionEditor::resetPoints()
{
    QPolygonF pts;
    int h_off = this->width() / 10;
    int v_off = this->height() / 8;
    pts << QPointF(this->width() / 2, this->height() / 2)
        << QPointF(this->width() / 2 - h_off, this->height() / 2 - v_off);

    foreach (ShadeWidget *s, _shades)
    {
        s->hoverPoints()->setPoints(pts);
    }
}

void TransferFunctionEditor::pointsUpdated()
{
    qreal w = _pAlphaShade->width();
    qreal h = _pAlphaShade->height();
    QGradientStops stops;
    QPolygonF points;
    QVector<QColor> colors;

    foreach (ShadeWidget *s, _shades)
    {
        colors.append(s->colors());
        points += s->points();
    }
    std::sort(points.begin(), points.end(), x_less_than);

    for (int i = 0; i < points.size(); ++i)
    {
        qreal x = int(points.at(i).x());
        qreal y = int(points.at(i).y());
        if (i + 1 < points.size() && x == points.at(i + 1).x())
            continue;
//        QColor color((0x00ff0000 & _pAlphaShade->colorAt(int(x))) >> 16,
//                     (0x0000ff00 & _pAlphaShade->colorAt(int(x))) >> 8,
//                     (0x000000ff & _pAlphaShade->colorAt(int(x))),
//                     (0xff000000 & _pAlphaShade->colorAt(int(x))) >> 24);
//        color = QColor(_pAlphaShade->colorAt(int(x)));

        if (x / w > 1)
            return;

        QColor col = colors.at(i);
        col.setAlphaF(1.f - y/h);
        stops << QGradientStop(x / w, col);
    }

    _pAlphaShade->setGradientStops(stops);
    _stops = stops;
    emit gradientStopsChanged(stops);
}


void TransferFunctionEditor::selectedPointUpdated(const QColor color)
{
    emit selectedPointChanged(color);
}


static void set_shade_points(const QPolygonF &points, ShadeWidget *shade)
{
    shade->hoverPoints()->setPoints(points);
    shade->hoverPoints()->setPointLock(0, HoverPoints::LockToLeft);
    shade->hoverPoints()->setPointLock(points.size() - 1, HoverPoints::LockToRight);
    shade->update();
}

static void setShadePointsColored(ShadeWidget *shade, const QPolygonF &points,
                                  const QVector<QColor> &colors)
{
    shade->hoverPoints()->setColoredPoints(points, colors);
    shade->hoverPoints()->setPointLock(0, HoverPoints::LockToLeft);
    shade->hoverPoints()->setPointLock(points.size() - 1, HoverPoints::LockToRight);
    shade->update();
}

void TransferFunctionEditor::setGradientStops(const QGradientStops &stops)
{
    _stops = stops;
//    QPolygonF pts_red, pts_green, pts_blue, pts_alpha;
    QPolygonF points;
    QVector<QColor> colors;

//    qreal h_red = _pRedShade->height();
//    qreal h_green = _pGreenShade->height();
//    qreal h_blue = _pBlueShade->height();
    qreal h_alpha = _pAlphaShade->height();

    for (int i = 0; i < stops.size(); ++i)
    {
        qreal pos = stops.at(i).first;
        QRgb color = stops.at(i).second.rgba();
        points << QPointF(pos * _pAlphaShade->width(), h_alpha - qAlpha(color) * h_alpha / 255);
        colors.push_back(color);

//        pts_red << QPointF(pos * _pRedShade->width(), h_red - qRed(color) * h_red / 255);
//        pts_green << QPointF(pos * _pGreenShade->width(), h_green - qGreen(color) * h_green / 255);
//        pts_blue << QPointF(pos * _pBlueShade->width(), h_blue - qBlue(color) * h_blue / 255);
//        pts_alpha << QPointF(pos * _pAlphaShade->width(), h_alpha - qAlpha(color) * h_alpha / 255);
    }

//    set_shade_points(pts_red, _pRedShade);
//    set_shade_points(pts_green, _pGreenShade);
//    set_shade_points(pts_blue, _pBlueShade);
//    set_shade_points(pts_alpha, _pAlphaShade);
    setShadePointsColored(_pAlphaShade, points, colors);
}

const QGradientStops TransferFunctionEditor::getGradientStops() const
{
    return _stops;
}


/**
 * @brief TransferFunctionEditor::setInterpolation
 * @param method
 */
void TransferFunctionEditor::setInterpolation(const QString method)
{
    if (method.contains("Quad"))
    {
        foreach (ShadeWidget *s, _shades)
             s->hoverPoints()->setConnectionType(HoverPoints::CurveConnection);
    }
    else if (method.contains("Linear"))
    {
        foreach (ShadeWidget *s, _shades)
             s->hoverPoints()->setConnectionType(HoverPoints::LineConnection);
    }
    foreach (ShadeWidget *s, _shades)
         s->hoverPoints()->firePointChange();

    emit pointsUpdated();
}


/**
 * @brief TransferFunctionEditor::setColorSelected
 * @param color
 */
void TransferFunctionEditor::setColorSelected(const QColor color)
{
//    foreach (ShadeWidget *s, _shades)
    _pAlphaShade->hoverPoints()->setColorSelected(color);
    emit pointsUpdated();
}


//------------------------------------------------------
/**
 * @brief TransferFunctionWidget::TransferFunctionWidget
 * @param parent
 */
TransferFunctionWidget::TransferFunctionWidget(QWidget *parent) : QWidget(parent)
{
    _pEditor = new TransferFunctionEditor();
    QVBoxLayout *mainLayout = new QVBoxLayout(this);
    mainLayout->addWidget(_pEditor);
}

/**
 * @brief TransferFunctionWidget::resetTransferFunction
 */
void TransferFunctionWidget::resetTransferFunction()
{
    QGradientStops stops;
    stops << QGradientStop(0.00, QColor::fromRgba(0xff000000));
    stops << QGradientStop(1.00, QColor::fromRgba(0xffffffff));
    _pEditor->setGradientStops(stops);
    _pEditor->pointsUpdated();
}


/**
 * @brief TransferFunctionWidget::setInterpolation
 * @param method
 */
void TransferFunctionWidget::setInterpolation(QString method)
{
    _pEditor->setInterpolation(method);
    _pEditor->pointsUpdated();
}


/**
 * @brief TransferFunctionWidget::setColorSelected
 * @param color
 */
void TransferFunctionWidget::setColorSelected(const QColor color)
{
   _pEditor->setColorSelected(color);
}


/**
 * @brief TransferFunctionWidget::getEditor
 * @return
 */
QPointer<TransferFunctionEditor> TransferFunctionWidget::getEditor()
{
   return _pEditor;
}
