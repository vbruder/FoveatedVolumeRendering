#ifndef STIPPLEVIEWER_H
#define STIPPLEVIEWER_H

#include <QGraphicsView>

#include "lbgstippling.h"

class StippleViewer : public QGraphicsView {

    Q_OBJECT

  public:
    StippleViewer(const QImage& img, QWidget* parent);
    void stipple(const LBGStippling::Params params);
    QPixmap getImage();
    void setInputImage(const QImage& img);
    void saveImageSVG(const QString& path);
    void saveImagePDF(const QString& path);
    void displayPoints(const std::vector<Stipple>& stipples);
    void displayCells(const IndexMap& cells);

protected slots:
    void saveImage();

  signals:
    void finished();
    void inputImageChanged();
    void iterationStatus(int iteration, int numberPoints, int splits, int merges, float hysteresis);

  private:
    LBGStippling m_stippling;
    QImage m_image;
    uint m_count = 0;
};

#endif // STIPPLEVIEWER_H
