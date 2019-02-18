#include "lbgstippling.h"

#include <cassert>

#include <QElapsedTimer>

#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonValue>
#include <QtMath>
#include <QPainter>
#include <algorithm>

#include <nanoflann.hpp>

// And this is the "dataset to kd-tree" adaptor class:

struct QVectorAdaptor {
    const QVector<QVector2D>& obj; //!< A const ref to the data set origin

    QVectorAdaptor(const QVector<QVector2D>& obj_) : obj(obj_) {}

    inline const QVector<QVector2D>& derived() const { return obj; }

    inline size_t kdtree_get_point_count() const { return derived().size(); }

    inline float kdtree_get_pt(const size_t idx, const size_t dim) const {
        const auto& point = derived()[idx];
        if (dim == 0)
            return point.x();
        else
            return point.y();
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const {
        return false;
    }
};

namespace Random {
std::random_device rd;
std::mt19937 gen(42);
} // namespace Random

using Params = LBGStippling::Params;
using Status = LBGStippling::Status;

QVector<QVector2D> sites(const std::vector<Stipple>& stipples) {
    QVector<QVector2D> sites(static_cast<int>(stipples.size()));
    std::transform(stipples.begin(), stipples.end(), sites.begin(),
                   [](const auto& s) { return s.pos; });
    return sites;
}

std::vector<Stipple> randomStipples(size_t n, float size) {
    std::uniform_real_distribution<float> dis(0.01f, 0.99f);
    std::vector<Stipple> stipples(n);
    std::generate(stipples.begin(), stipples.end(), [&]() {
        return Stipple{QVector2D(dis(Random::gen), dis(Random::gen)), size, Qt::black};
    });
    return stipples;
}

template <class T>
inline T pow2(T x) {
    return x * x;
}

QVector2D jitter(QVector2D s) {
    using namespace Random;
    std::uniform_real_distribution<float> jitter_dis(-0.001f, 0.001f);
    return s += QVector2D(jitter_dis(gen), jitter_dis(gen));
}

float getSplitValueUpper(float pointDiameter, float hysteresis, size_t superSampling) {
    const float pointArea = M_PI * pow2(pointDiameter / 2.0f);
    return (1.0f + hysteresis / 2.0f) * pointArea * pow2(superSampling);
}

float getSplitValueLower(float pointDiameter, float hysteresis, size_t superSampling) {
    const float pointArea = M_PI * pow2(pointDiameter / 2.0f);
    return (1.0f - hysteresis / 2.0f) * pointArea * pow2(superSampling);
}

float stippleSize(const VoronoiCell& cell, const Params& params) {
    if (params.adaptivePointSize) {
        float norm;
        switch (params.mapping) {
        case LBGStippling::PointMappingFunction::LINEAR:
            norm = cell.sumDensity / cell.area;
            break;
        case LBGStippling::PointMappingFunction::SQUAREROOT:
            norm = std::sqrt(cell.sumDensity / cell.area);
            break;
        case LBGStippling::PointMappingFunction::EXPONENTIAL:
            norm = 2.0f - 2.0f / (1.0f + cell.sumDensity / cell.area);
            break;
        case LBGStippling::PointMappingFunction::SQUARE:
            norm = pow2(cell.sumDensity / cell.area);
        };
        return params.pointSizeMin * (1.0f - norm) + params.pointSizeMax * norm;
    } else {
        return params.initialPointSize;
    }
}

float currentHysteresis(size_t i, const Params& params) {
    return params.hysteresis + i * params.hysteresisDelta;
}

bool notFinished(const Status& status, const Params& params) {
    auto [iteration, size, splits, merges, hysteresis] = status;
    return !((splits == 0 && merges == 0) || (iteration == params.maxIterations));
}

LBGStippling::LBGStippling() {
    m_statusCallback = [](const Status&) {};
    m_stippleCallback = [](const std::vector<Stipple>&) {};
    m_cellCallback = [](const IndexMap&) {};
}

void LBGStippling::setStatusCallback(Report<Status> statusCB) {
    m_statusCallback = statusCB;
}

void LBGStippling::setStippleCallback(Report<std::vector<Stipple>> stippleCB) {
    m_stippleCallback = stippleCB;
}

void LBGStippling::setCellCallback(Report<IndexMap> cellCB) {
    m_cellCallback = cellCB;
}

/**
 * @brief encodeMorton
 * @see http://www.graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN
 * @param p point
 * @param width Index map width
 * @param height Index map height
 * @return morton encoded index
 */
unsigned int encodeMorton(QVector2D p, size_t width, size_t height)
{
    static const unsigned int B[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};
    static const unsigned int S[] = {1, 2, 4, 8};

    // x and y must initially be less than 65536.
    unsigned int x = static_cast<unsigned int>(p.x() * width); // Interleave lower 16 bits of x and y, so the bits of x
    unsigned int y = static_cast<unsigned int>(p.y() * height); // are in the even positions and bits from y in the odd;

    x = (x | (x << S[3])) & B[3];
    x = (x | (x << S[2])) & B[2];
    x = (x | (x << S[1])) & B[1];
    x = (x | (x << S[0])) & B[0];

    y = (y | (y << S[3])) & B[3];
    y = (y | (y << S[2])) & B[2];
    y = (y | (y << S[1])) & B[1];
    y = (y | (y << S[0])) & B[0];

    // return the resulting 32-bit Morton Number.
    return  x | (y << 1);
}

LBGStippling::Result LBGStippling::stipple(const QImage& density, const Params& params,
                                           const int batchCount, const int batchNo) const {
    QImage densityGray =
        density
            .scaledToWidth(static_cast<int>(params.superSamplingFactor) * density.width(),
                           Qt::SmoothTransformation)
            .convertToFormat(QImage::Format_Grayscale8);

    VoronoiDiagram voronoi(densityGray);

    std::vector<Stipple> stipples = randomStipples(params.initialPoints, params.initialPointSize);

    IndexMap indexMap;

    Status status = {0, 0, 1, 1, 0.5f};

    QVector<QVector2D> points;

    qDebug() << "LBG: Starting...";

    while (notFinished(status, params)) {
        status.splits = 0;
        status.merges = 0;

        points = sites(stipples);
        if (points.empty()) {
            // Stop if image is empty (all white).
            m_stippleCallback(stipples);
            break;
        }

        indexMap = voronoi.calculate(points);
        std::vector<VoronoiCell> cells = accumulateCells(indexMap, densityGray);

        assert(cells.size() == stipples.size());

        stipples.clear();

        float hysteresis = currentHysteresis(status.iteration, params);
        status.hysteresis = hysteresis;

        for (const auto& cell : cells) {
            float totalDensity = cell.sumDensity;
            float diameter = stippleSize(cell, params);

            if (totalDensity <
                    getSplitValueLower(diameter, hysteresis, params.superSamplingFactor) ||
                cell.area == 0.0f) {
                // cell too small - merge
                ++status.merges;
                continue;
            }

            if (totalDensity <
                getSplitValueUpper(diameter, hysteresis, params.superSamplingFactor)) {
                // cell size within acceptable range - keep
                stipples.push_back({cell.centroid, diameter, Qt::black});
                continue;
            }

            // cell too large - split
            float area = std::max(1.0f, cell.area);
            float circleRadius = std::sqrt(area / M_PI);
            QVector2D splitVector = QVector2D(0.5f * circleRadius, 0.0f);

            float a = cell.orientation;
            QVector2D splitVectorRotated =
                QVector2D(splitVector.x() * std::cos(a) - splitVector.y() * std::sin(a),
                          splitVector.y() * std::cos(a) + splitVector.x() * std::sin(a));

            splitVectorRotated.setX(splitVectorRotated.x() / densityGray.width());
            splitVectorRotated.setY(splitVectorRotated.y() / densityGray.height());

            QVector2D splitSeed1 = cell.centroid - splitVectorRotated;
            QVector2D splitSeed2 = cell.centroid + splitVectorRotated;

            // check boundaries
            splitSeed1.setX(std::max(0.0f, std::min(splitSeed1.x(), 1.0f)));
            splitSeed1.setY(std::max(0.0f, std::min(splitSeed1.y(), 1.0f)));

            splitSeed2.setX(std::max(0.0f, std::min(splitSeed2.x(), 1.0f)));
            splitSeed2.setY(std::max(0.0f, std::min(splitSeed2.y(), 1.0f)));

            stipples.push_back({jitter(splitSeed1), diameter, Qt::red});
            stipples.push_back({jitter(splitSeed2), diameter, Qt::red});

            ++status.splits;
        }
        status.size = stipples.size();
        m_statusCallback(status);
        m_cellCallback(indexMap);
        m_stippleCallback(stipples);

        ++status.iteration;
    }
    qDebug() << "LBG: Done";
    qDebug() << "Number of points: " << points.size();

    // stipples+points -> z ordering map
    QMap<unsigned int, unsigned int> mortonMap;
    for (unsigned int i = 0; i < points.size(); ++i)
    {
        unsigned int mc = encodeMorton(points.at(i), indexMap.width, indexMap.height);
        mortonMap.insert(mc, i);
    }

    // debug: draw morton ordering
    QImage dbgImg = QPixmap(2048,2048).toImage();
    dbgImg.fill(Qt::white);
    QPainter dbgPainter(&dbgImg);
    QVector<QPointF> dbgPoints;
    for (auto &a: mortonMap.values())
    {
        dbgPoints.push_back(points.at(a).toPointF()*2048.);
    }
    dbgPainter.drawPolyline(dbgPoints);
    dbgPainter.end();
    dbgImg.save(QString::number(points.size()) + "morton.png");

    // map: point id -> morton order id (compacted)
    QList<unsigned int> point2morton = mortonMap.values();

    const size_t batchSize = indexMap.height / batchCount;
    const size_t BucketCount = 16;
    std::vector<uint32_t> neighborIndexMap;
    std::vector<float> neighborWeightMap;

    neighborIndexMap.resize(indexMap.width * batchSize * BucketCount, 0);
    neighborWeightMap.resize(indexMap.width * batchSize * BucketCount, 0.0f);

#if true
    qDebug() << "Starting batch" << batchNo << "/" << batchCount << "with size" << batchSize
             << indexMap.width << "x" << indexMap.height;

    using namespace nanoflann;
    typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<float, QVectorAdaptor>, QVectorAdaptor, 2>
        my_kd_tree_t;

    QVectorAdaptor adaptor(points);
    my_kd_tree_t index(2, adaptor, KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    QElapsedTimer progressTimer;
    progressTimer.start();

//    QElapsedTimer perfTimer;
//    perfTimer.start();

    const size_t k = BucketCount; // + 4;
    std::vector<size_t> ret_indices(k);
    std::vector<float> out_dists_sqr(k);

    for (int y = batchNo*batchSize; y < batchNo*batchSize + batchSize; ++y)
    {
        int yOffset = (y - batchNo * batchSize);
        for (int x = 0; x < indexMap.width; ++x)
        {
            if (x == 1 && (y % 10 == 0))
            {
                float pointsTotal = static_cast<float>(indexMap.width * batchSize);
                float pointsProgress = static_cast<float>(yOffset * indexMap.width + x) / pointsTotal;
                auto elapsedTime = progressTimer.elapsed();
                auto remainingTime = ((1. - pointsProgress) / pointsProgress) * elapsedTime;
                qDebug() << "Natural Neighbor" << qSetRealNumberPrecision(5) << pointsProgress
                         << elapsedTime / 1000. << "sec" << remainingTime / 1000. / 60. << "min";
            }

            // Insert point at pixel to compute the modified voronoi diagram.
            QVector2D modifierPoint(static_cast<float>(x) / indexMap.width,
                                    static_cast<float>(y) / indexMap.height);
            QVector<QVector2D> pointsModified;

            float query_pt[2] = {modifierPoint.x(), modifierPoint.y()};
            nanoflann::KNNResultSet<float> resultSet(k);
            resultSet.init(&ret_indices[0], &out_dists_sqr[0]);
            bool yay = index.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams());
            assert(yay && "Ohne KNN gehts halt net");
            for (size_t i = 0; i < k; ++i) { pointsModified.append(points[ret_indices[i]]); }

            pointsModified.append(modifierPoint);
            IndexMap indexMapModified = voronoi.calculate(pointsModified);

            // Count intersection for each cell in the original index map.
            float intersectionSum = 0;
            QMap<uint32_t, float> intersectionSet;
            uint32_t modifierPointIndex = indexMapModified.get(x, y);

            int kernelHeight = 40;  //indexMapModified.height /4;
            int kernelWidth = 40;   //indexMapModified.width /4;
            bool overflow = false;
//#pragma omp parallel for
            for (int yy = std::max(yOffset - kernelHeight, 0); yy < std::min(y + kernelHeight, static_cast<int>(indexMapModified.height)); ++yy)
            {
                for (int xx = std::max(x - kernelWidth, 0); xx < std::min(x + kernelWidth, static_cast<int>(indexMapModified.width)); ++xx)
                {
                    if (indexMapModified.get(xx, yy) == modifierPointIndex)
                    {
                        if (intersectionSet.size() >= BucketCount)
                        {
                            qDebug() << "__Bucket overflow on pixel" << x << y;
                            overflow = true;
                        }
                        else
                        {
                            intersectionSum++;
                        }
                        // encode in morton order
                        auto originalIndex = point2morton.value(indexMap.get(xx, yy));
                        intersectionSet[originalIndex]++;
                    }
                }
            }            
//            qDebug() << "# pixels" << intersectionSum << "time " << perfTimer.restart();
//            assert(intersectionSet.size() < BucketCount && "Mehr geht halt net erstmal...");

            auto smap = intersectionSet.toStdMap();
            while (smap.size() > BucketCount)
            {
                auto it = min_element(smap.begin(), smap.end(),
                                      [](decltype(smap)::value_type & l, decltype(smap)::value_type& r) -> bool { return l.second < r.second; });
                smap.erase(it);
            }
            intersectionSet = QMap<uint32_t, float>(smap);

            // Normalize weights and copy to neighbor maps.
            size_t bucketIndex = 0;
            for (auto key : intersectionSet.keys())
            {
                intersectionSet[key] = intersectionSet[key] / intersectionSum;
                auto offset = (yOffset * indexMap.width + x) * BucketCount + bucketIndex;
                neighborIndexMap[offset] = key;
                neighborWeightMap[offset] = intersectionSet[key];
                bucketIndex++;
            }
        }
    }

    qDebug() << "Natural Neighbor: Done";
    qDebug() << "Batch no" << batchNo << " / " << batchCount;
    qDebug() << "Number of points: " << points.size() << params.pointSizeMin << params.pointSizeMax;
#endif
    return {stipples, indexMap, neighborWeightMap, neighborIndexMap, point2morton};
}

void LBGStippling::Params::saveParametersJSON(const QString& path) {
    QJsonObject json;
    json.insert("initialPoints", QJsonValue(initialPoints));
    json.insert("initialPointSize", QJsonValue(initialPointSize));
    json.insert("adaptivePointSize", QJsonValue(adaptivePointSize));
    json.insert("pointSizeMin", QJsonValue(pointSizeMin));
    json.insert("pointSizeMax", QJsonValue(pointSizeMax));
    QString sMapping;
    using pmf = LBGStippling::PointMappingFunction;
    switch (mapping) {
    case pmf::LINEAR:
        sMapping = "Linear";
        break;
    case pmf::SQUAREROOT:
        sMapping = "SquareRoot";
        break;
    case pmf::EXPONENTIAL:
        sMapping = "Exponential";
        break;
    case pmf::SQUARE:
        sMapping = "Square";
    };
    json.insert("pointSizeMapping", QJsonValue(sMapping));
    json.insert("superSamplingFactor", QJsonValue((int)superSamplingFactor));
    json.insert("maxIterations", QJsonValue((int)maxIterations));
    json.insert("hysteresis", QJsonValue(hysteresis));
    json.insert("hysteresisDelta", QJsonValue(hysteresisDelta));

    QJsonDocument jsonDoc(json);
    QFile jsonFile(path + ".json");
    jsonFile.open(QFile::WriteOnly);
    jsonFile.write(jsonDoc.toJson());
    jsonFile.close();
}

void LBGStippling::Params::loadParametersJSON(const QString& path) {
    QFile jsonFile(path);
    jsonFile.open(QFile::ReadOnly);
    QJsonDocument jsonDoc = QJsonDocument::fromJson(jsonFile.readAll());
    jsonFile.close();
    QJsonObject json = jsonDoc.object();

    initialPoints = json["initialPoints"].toInt();
    initialPointSize = json["initialPointSize"].toDouble();
    adaptivePointSize = json["adaptivePointSize"].toBool();
    pointSizeMin = json["pointSizeMin"].toDouble();
    pointSizeMax = json["pointSizeMax"].toDouble();
    QString sMapping = json["pointSizeMapping"].toString();
    using pmf = LBGStippling::PointMappingFunction;
    if (sMapping == "Linear")
        mapping = pmf::LINEAR;
    else if (sMapping == "SquareRoot")
        mapping = pmf::SQUAREROOT;
    else if (sMapping == "Exponential")
        mapping = pmf::EXPONENTIAL;
    else if (sMapping == "Square")
        mapping = pmf::SQUARE;
    superSamplingFactor = json["superSamplingFactor"].toInt();
    maxIterations = json["maxIterations"].toInt();
    hysteresis = json["hysteresis"].toDouble();
    hysteresisDelta = json["hysteresisDelta"].toDouble();
}
